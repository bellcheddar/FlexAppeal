"""Resolving a construct's motifs, core and reference ensembles.

Everything here needs a catalogue: KLIFS pocket positions for a kinase or
GPCRdb generic numbers for a receptor, the deposited structures behind them,
and `gemmi` to read those. None of it is per frame, and none of it can run on a
machine that has an OpenMM environment and nothing else.

That is why it is a module rather than part of the run. It was written inside
`tools/analyse_run.py`, where `flexappeal/refstates/cli.py` reached it by
loading the script from a path, and it is imported from both places now: the
Prepare tab calls it when a bundle is built, which is where the network and the
caches already are, and the bundle carries the few hundred kilobytes of motif
coordinates that come out so that scoring a frame needs nothing but numpy.
"""

from __future__ import annotations

import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

KLIFS_CACHE = ROOT / "data" / "cache" / "klifs"
GPCRDB_CACHE = ROOT / "data" / "cache" / "gpcrdb"

# The kinase motifs section 5.4 asks to be scored, by KLIFS region name. The
# pocket scaffold is excluded on purpose: it is the superposition reference, and
# scoring a motif against the frame it was superposed on would report zero.
KINASE_MOTIFS = ("g.l", "alphaC", "hinge", "c.l", "xDFG", "A-loop")

# Section 5.5's receptor motifs, by GPCRdb segment, plus the NPxxY stretch which
# is a span of generic numbers rather than a segment.
GPCR_MOTIFS = ("TM5", "TM6", "TM7", "H8", "NPxxY")

# The key the run's own structure takes in an intersection alongside the
# references. A name no PDB code can collide with, because it sits in the same
# dictionary as them.
RUN_CORE_KEY = "__run__"


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def _best_entries(entries: list, limit: int | None, *, key: str,
                  ascending: bool = False) -> list:
    """The `limit` best catalogue rows, or all of them when there is no limit.

    The worked examples need no limit: ABL1 yields 12 active and 45 inactive,
    beta2AR 72 and 28, and every committed result is computed over all of them.
    An arbitrary protein is a different matter. CDK2 yields 337 active, and
    since a state score is the distance to the *nearest* reference, the three
    hundredth structure adds a download and almost nothing else. So the reuse
    path caps it and the production path does not, which keeps the committed
    numbers exactly as they were.

    Sorted by the catalogue's own quality column where it has one, descending,
    and left in catalogue order where it does not, so the choice is the best
    structures rather than the first ones.
    """
    if not limit or len(entries) <= limit:
        return list(entries)

    def score(row) -> float:
        try:
            return float(row.get(key))
        except (TypeError, ValueError):
            return float("-inf")

    if any(row.get(key) is not None for row in entries):
        # A missing value sorts last either way, which is why the fallback is
        # -inf for a descending sort and +inf for an ascending one.
        if ascending:
            def rank(row) -> float:
                try:
                    return float(row.get(key))
                except (TypeError, ValueError):
                    return float("inf")
            return sorted(entries, key=rank)[:limit]
        return sorted(entries, key=score, reverse=True)[:limit]
    return list(entries)[:limit]


def kinase_setup(start_pdb: Path, reference_pdb: str,
                 protein: str = "ABL1", max_references: int | None = None):
    """Motif residues, core and reference ensembles for a kinase run.

    `protein` is the name KLIFS knows the kinase by, and it has to reach the
    lookup: this took the default and ignored its caller, so the reuse hook
    annotated every kinase against ABL1's reference ensembles and reported a
    confident state for it. That is section 1's G-REUSE goal failing silently,
    which is the worst way for it to fail, because the output looks the same.
    """
    from flexappeal.refstates import Catalogue, kinase_id, kinase_references
    from flexappeal.refstates import ensembles, pocket, refs, sse
    from flexappeal.refstates.traj import (core_agreement, intersect_core,
                                           scaffold_core)

    catalogue = Catalogue(KLIFS_CACHE)
    kid = kinase_id(catalogue, name=protein)
    ids = pocket.structure_ids(catalogue, kid)
    log(f"KLIFS {protein} is kinase {kid}, {len(ids)} structures")

    import gemmi

    structure = gemmi.read_structure(str(start_pdb))
    structure.setup_entities()
    chain = structure[0][0]

    # The reference entry whose pocket map resolves this construct's motif
    # residues has to be one of *this* kinase's structures. The default is an
    # ABL1 entry, so for any other kinase it is absent and this used to raise a
    # bare KeyError naming a PDB code with no hint that the code belonged to a
    # different protein.
    code = (reference_pdb or "").upper()
    if code and code in ids:
        chosen = code
    elif code:
        raise SystemExit(
            f"--reference-pdb {code} is not one of KLIFS's {len(ids)} "
            f"structures for {protein}. Pass one of that kinase's own entries, "
            f"for example {', '.join(sorted(ids)[:4])}, or leave it unset to "
            f"let one be chosen.")
    elif ids:
        # Deterministic rather than arbitrary, and said out loud: any structure
        # of this kinase resolves the pocket, but which one was used belongs in
        # the log.
        chosen = sorted(ids)[0]
        log(f"no --reference-pdb given; resolving the pocket from {chosen}")
    else:
        raise SystemExit(f"KLIFS has no structures for {protein}")

    own = pocket.pocket_map(catalogue, ids[chosen])
    last = max(r.seqid.num for r in chain)
    motif_residues: dict[str, list[int]] = {}
    for name in KINASE_MOTIFS:
        if name == "A-loop":
            residues = pocket.activation_loop(own, chain, last_residue=last)
        else:
            residues = own.region(name)
        present = [n for n in residues if any(r.seqid.num == n for r in chain)]
        if present:
            motif_residues[name] = present

    core = scaffold_core(own)
    ensemble_rows = kinase_references(catalogue, kid)

    # Each reference resolves its own residue numbers, which is what makes the
    # ABL1a and ABL1b entries comparable at all: 3PYY is 19 residues offset and
    # a shared residue list would read every motif at the wrong place in it.
    maps: dict[str, pocket.PocketMap] = {}

    def pocket_for(code: str):
        if code not in maps:
            maps[code] = pocket.pocket_map(catalogue, ids[code])
        return maps[code]

    rows: dict[str, list] = {}
    paths: dict[str, Path] = {}
    for which in ("active", "inactive"):
        available = _best_entries(ensemble_rows[which].entries, max_references,
                                  key="quality_score")
        entries = []
        for entry in available:
            code = str(entry["pdb"]).upper()
            if code not in ids:
                continue
            paths[code] = refs.fetch_mmcif(code, KLIFS_CACHE)
            entries.append({"pdb": code, "path": str(paths[code]),
                            "chain": entry.get("chain") or ""})
        rows[which] = entries
        total = len(ensemble_rows[which].entries)
        # Rows and structures are not the same number and the difference is
        # large. KLIFS lists one row per chain, and the ensemble is keyed on the
        # deposition, so ABL1's "12 active references" are 5 structures and its
        # 45 inactive are 14. The score rests on the structures, so both are
        # reported rather than the flattering one.
        unique = len({e["pdb"] for e in entries})
        best = (f", the best of {total} by quality"
                if max_references and total > len(entries) else "")
        log(f"{which}: {len(entries)} catalogue rows over {unique} "
            f"structures{best}")

    # The activation loop runs from the DFG aspartate to the APE motif, and APE
    # is found in each structure's own sequence, so a reference needs its own
    # chain to resolve it. Computed without one it falls back to KLIFS's `a.l`,
    # which is two residues, and every A-loop was then dropped for a length
    # mismatch against the trajectory's 26 without anything saying so.
    aloop_cache: dict[str, list[int]] = {}

    def aloop_for(code: str) -> list[int]:
        if code not in aloop_cache:
            import gemmi as _gemmi

            found = pocket_for(code)
            try:
                reference = _gemmi.read_structure(str(paths[code]))
                reference.setup_entities()
                aloop_cache[code] = pocket.activation_loop(
                    found, reference[0][0])
            except Exception:                                   # noqa: BLE001
                aloop_cache[code] = pocket.activation_loop(found)
        return aloop_cache[code]

    def motifs_for(code: str):
        """What this reference's own numbering calls each motif.

        Returned per code rather than shared, which is what makes the ABL1a and
        ABL1b entries comparable: 3PYY is offset by 19 residues and one shared
        list would read every motif at the wrong place in it.
        """
        found = pocket_for(code)
        out: dict[str, list[int]] = {}
        for name, residues in motif_residues.items():
            if name == "A-loop":
                # The one motif whose length is not fixed by the KLIFS
                # alignment. Both spans start at the DFG aspartate, so taking
                # the same count from that shared anchor is a real residue to
                # residue correspondence; trimming the tail is only correct
                # because the anchor is shared.
                theirs = aloop_for(code)[:len(residues)]
            else:
                theirs = found.region(name)
            if theirs and len(theirs) == len(residues[:len(theirs)]):
                out[name] = theirs
        return out

    # D025's deferred fix, built. The core is the KLIFS scaffold positions that
    # the run's own structure and every reference all resolve, taken over
    # positions rather than residue numbers because ABL1's references are split
    # between the ABL1a and ABL1b numbering and intersecting numbers across them
    # collapses to nothing. Resolved back per structure, so the i-th core
    # coordinate is the same place in the fold everywhere, which is what the
    # superposition was assuming all along.
    #
    # A no-op on a set that already resolves the same positions, which both
    # production ensembles do; it is CDK2, with references resolving 38 of 42,
    # that this is for.
    shared_core = intersect_core(
        {RUN_CORE_KEY: own,
         **{code: pocket_for(code) for code in sorted(paths)}})
    if shared_core["positions"]:
        kept = shared_core["residues_by"]
        narrowed = set(kept[RUN_CORE_KEY])
        if narrowed != set(core.residues):
            core.source += (f", intersected to {len(shared_core['positions'])} "
                            f"positions shared by {shared_core['n_structures']} "
                            f"structures")
            core.helices = {region: [n for n in residues if n in narrowed]
                            for region, residues in core.helices.items()}
            core.helices = {r: v for r, v in core.helices.items() if v}
            log(f"core intersected to {len(shared_core['positions'])} KLIFS "
                f"positions shared by {shared_core['n_structures']} structures, "
                f"from {len(core.residues)} on the start model")
        core.residues = list(kept[RUN_CORE_KEY])

    def core_for(code: str):
        kept = shared_core["residues_by"].get(code)
        if kept:
            return kept
        return scaffold_core(pocket_for(code)).residues

    # No coverage intersection here, unlike the receptor path, and the reason is
    # worth stating because trying it made things much worse. ABL1's references
    # are split between the ABL1a and ABL1b numbering, 19 residues apart, so
    # intersecting residue *numbers* across them compares 290 in one structure
    # with 309 in another and the intersection collapses to nothing: every motif
    # but one was dropped. The receptor's references all share one numbering,
    # which is what makes the same operation safe there.
    #
    # So each reference keeps its own numbering and `build` drops a reference
    # for a motif it does not fully resolve, keeping it for the others.
    built = ensembles.build(rows, motif_residues_for=motifs_for,
                            core_residues_for=core_for)

    # Section 5.6's SSE comparison against the reference ensembles, built per
    # motif because that is the only level at which the kinase references are
    # reliably aligned: each one's residues come from its own pocket map.
    sse_profiles = sse.motif_profiles(rows, motif_residues_for=motifs_for,
                                      motifs=sorted(motif_residues),
                                      cache_dir=GPCRDB_CACHE)
    log("  reference secondary structure: "
        + ", ".join(f"{state} {len(by_motif)} motifs"
                    for state, by_motif in sorted(sse_profiles.items())))

    alignment = ensembles.align_to_common_frame(built)
    log(f"aligned {alignment['n_aligned']} references onto "
        f"{alignment['target']}: core residual mean "
        f"{alignment['mean_core_rmsd_a']} A, max {alignment['max_core_rmsd_a']} A"
        + (f", {alignment['n_dropped']} dropped" if alignment['n_dropped'] else ""))

    cores = {code: motifs[ensembles.CORE_KEY]
             for state in built for code, motifs in built[state].by_code.items()}
    agreement = core_agreement(cores)
    log(f"core agreement across {agreement.get('n')} references: mean "
        f"{agreement.get('mean_rmsd_a')} A, worst {agreement.get('max_rmsd_a')} A "
        f"({agreement.get('worst')})")
    return motif_residues, core, built, own, {"rigidity": agreement,
                                              "alignment": alignment,
                                              "sse_profiles": sse_profiles}, \
        {"rows_by_state": rows, "motifs_for": motifs_for,
         "core_for": core_for}


def gpcr_setup(start_pdb: Path, protein: str = "adrb2_human",
               max_references: int | None = None):
    """Motif residues, core and reference ensembles for a receptor run.

    Every residue comes through GPCRdb generic numbers, so the segments are the
    catalogue's own and the beta2AR numbers are never typed. NPxxY is the one
    motif that is a span of generic numbers rather than a segment, because that
    is what it is: 7.49 to 7.53.

    `protein` is GPCRdb's entry name, `adrb2_human` for the worked example. Like
    the kinase side, this ignored its caller and every receptor was annotated
    against beta2AR.
    """
    import gemmi

    from flexappeal.refstates import (Catalogue, ensembles, receptor_references,
                                      receptor_residues, refs, sse)
    from flexappeal.refstates.gpcr import NPXXY, build_numbering
    from flexappeal.refstates.traj import core_agreement, gpcr_core

    catalogue = Catalogue(GPCRDB_CACHE)
    numbering = build_numbering(receptor_residues(catalogue, protein))

    structure = gemmi.read_structure(str(start_pdb))
    structure.setup_entities()
    present = {r.seqid.num for chain in structure[0] for r in chain}

    motif_residues: dict[str, list[int]] = {}
    for segment in GPCR_MOTIFS:
        if segment == "NPxxY":
            lo, hi = (int(float(NPXXY[0]) * 100) % 100,
                      int(float(NPXXY[1]) * 100) % 100)
            residues = [numbering.residue(f"7.{n}") for n in range(lo, hi + 1)]
            residues = [n for n in residues if n is not None and n in present]
        else:
            residues = [n for n in numbering.segment_residues(segment)
                        if n in present]
        if residues:
            motif_residues[segment] = residues

    core = gpcr_core(numbering, present)
    ensemble_rows = receptor_references(catalogue, protein)

    # One numbering for the whole receptor: unlike the kinase references, every
    # beta2AR entry here is checked by `ca_distance` against the catalogue's own
    # residue identities, so a renumbered construct is refused rather than
    # silently measured at the wrong residues.
    def motifs_for(code: str):
        return motif_residues

    def core_for(code: str):
        return core.residues

    rows = {}
    for which in ("active", "inactive"):
        entries = []
        # GPCRdb has no quality column, so a cap takes the best resolution.
        available = _best_entries(ensemble_rows[which].entries, max_references,
                                  key="resolution", ascending=True)
        # `receptor_references` returns Ensemble objects, not lists, and the
        # GPCRdb row names the entry `pdb_code` where KLIFS uses `pdb`.
        for entry in available:
            code = str(entry.get("pdb_code") or "").upper()
            if not code:
                continue
            try:
                path = refs.fetch_mmcif(code, GPCRDB_CACHE)
            except Exception as exc:                            # noqa: BLE001
                log(f"  {code}: not fetched ({str(exc)[:50]})")
                continue
            entries.append({"pdb": code, "path": str(path),
                            "chain": str(entry.get("preferred_chain") or "")[:1]})
        rows[which] = entries
        total = len(ensemble_rows[which].entries)
        unique = len({e["pdb"] for e in entries})
        best = (f", the best of {total} by resolution"
                if max_references and total > len(entries) else "")
        log(f"{which}: {len(entries)} catalogue rows over {unique} "
            f"structures{best}")

    # Restricted to what every reference resolves. Without this TM6, the
    # headline activation motif, went unscored: 38 residues through the most
    # mobile helix of the receptor and no structure resolves all of them.
    common, coverage = ensembles.common_motif_residues(
        rows, motif_residues_for=motifs_for)
    log(f"  motif coverage over {coverage['n_references']} references at "
        f"{coverage['coverage']:.0%}:")
    for name, counts in sorted(coverage["motifs"].items()):
        state = "dropped" if name in coverage["dropped"] else "kept"
        log(f"    {name}: {counts['kept']} of {counts['requested']} residues "
            f"across {counts['references_with_motif']} references that have it, {state}")
    motif_residues = {k: v for k, v in common.items()}

    def motifs_for(code: str):                      # noqa: F811
        return motif_residues

    built = ensembles.build(rows, motif_residues_for=motifs_for,
                            core_residues_for=core_for)
    sse_profiles = sse.motif_profiles(rows, motif_residues_for=motifs_for,
                                      motifs=sorted(motif_residues),
                                      cache_dir=GPCRDB_CACHE)
    log("  reference secondary structure: "
        + ", ".join(f"{state} {len(by_motif)} motifs"
                    for state, by_motif in sorted(sse_profiles.items())))
    alignment = ensembles.align_to_common_frame(built)
    log(f"aligned {alignment['n_aligned']} references onto "
        f"{alignment['target']}: core residual mean "
        f"{alignment['mean_core_rmsd_a']} A, max {alignment['max_core_rmsd_a']} A")

    cores = {code: motifs[ensembles.CORE_KEY]
             for state in built for code, motifs in built[state].by_code.items()}
    rigidity = core_agreement(cores)
    log(f"core agreement across {rigidity.get('n')} references: mean "
        f"{rigidity.get('mean_rmsd_a')} A, worst {rigidity.get('max_rmsd_a')} A")

    # Section 5.6's second map axis: NPxxY RMSD to the nearest inactive
    # reference. Without it the receptor map is a row of points with no y and
    # the panel renders empty while the endpoint returns 200.
    npxxy = ensembles.rmsd_to_nearest(built, "NPxxY", "inactive")
    if npxxy:
        log(f"  map axis: NPxxY RMSD to nearest inactive for "
            f"{len(npxxy)} references, "
            f"{min(npxxy.values()):.2f} to {max(npxxy.values()):.2f} A")
    return (motif_residues, core, built, numbering,
            {"rigidity": rigidity, "alignment": alignment,
             "sse_profiles": sse_profiles,
             "map_axis": {"motif": "NPxxY", "against": "inactive",
                          "values": npxxy}},
            {"rows_by_state": rows, "motifs_for": motifs_for,
             "core_for": core_for})


# Section 6.3: coordinates for the web, full resolution kept in data/md. The
# stride is the trajectory's own frame interval, so the viewer gets every frame
# that was written.
#
# It was 1.0 ns, and the comment above it said "a 50 ns run is then 500 frames".
# A 1 ns stride over 50 ns is fifty frames, not five hundred: the sentence
# justifying the number was out by a factor of ten, and nothing downstream
# recomputed it. Fifty frames is a frame per nanosecond, which is a jerky
# viewer by construction, and with playback pacing a fixed-length pass over
# however many frames it had, it also ran at one frame per second. Reported as
# "jerky and showing slow motion", which is exactly what those two faults
# multiply out to. 4bt.
#
# At 0.1 ns the kinase bundle is about 8 MB and the receptor 10 MB against
# section 6.3's 150 MB ceiling, so the budget the old number was protecting was
# never at risk.
VIEWER_STRIDE_NS = 0.1
