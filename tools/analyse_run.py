#!/usr/bin/env python3
"""Analyse one finished or in-flight production run.

Sections 5.4 to 5.6 end to end: resolve the motifs, build the reference motif
coordinates, superpose each frame on the core and score it. Writes one JSON per
run, which is what the Runs library in section 6 reads.

Safe to run on a trajectory that is still being written, which is what section
2.3's streaming analysis needs: mdtraj reads whatever frames are on disk, and a
later pass simply sees more of them.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

KLIFS_CACHE = ROOT / "data" / "cache" / "klifs"
GPCRDB_CACHE = ROOT / "data" / "cache" / "gpcrdb"

# The kinase motifs section 5.4 asks to be scored, by KLIFS region name. The
# pocket scaffold is excluded on purpose: it is the superposition reference, and
# scoring a motif against the frame it was superposed on would report zero.
KINASE_MOTIFS = ("g.l", "alphaC", "hinge", "c.l", "xDFG", "A-loop")

# Section 5.5's receptor motifs, by GPCRdb segment, plus the NPxxY stretch which
# is a span of generic numbers rather than a segment.
GPCR_MOTIFS = ("TM5", "TM6", "TM7", "H8", "NPxxY")


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def kinase_setup(start_pdb: Path, reference_pdb: str):
    """Motif residues, core and reference ensembles for a kinase run."""
    from flexappeal.refstates import Catalogue, kinase_id, kinase_references
    from flexappeal.refstates import ensembles, pocket, refs
    from flexappeal.refstates.traj import core_agreement, scaffold_core

    catalogue = Catalogue(KLIFS_CACHE)
    kid = kinase_id(catalogue)
    ids = pocket.structure_ids(catalogue, kid)
    log(f"KLIFS kinase {kid}, {len(ids)} structures")

    import gemmi

    structure = gemmi.read_structure(str(start_pdb))
    structure.setup_entities()
    chain = structure[0][0]

    own = pocket.pocket_map(catalogue, ids[reference_pdb.upper()])
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
        entries = []
        for entry in ensemble_rows[which].entries:
            code = str(entry["pdb"]).upper()
            if code not in ids:
                continue
            paths[code] = refs.fetch_mmcif(code, KLIFS_CACHE)
            entries.append({"pdb": code, "path": str(paths[code]),
                            "chain": entry.get("chain") or ""})
        rows[which] = entries
        log(f"{which}: {len(entries)} references")

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

    def core_for(code: str):
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
                                              "alignment": alignment}


def gpcr_setup(start_pdb: Path):
    """Motif residues, core and reference ensembles for a receptor run.

    Every residue comes through GPCRdb generic numbers, so the segments are the
    catalogue's own and the beta2AR numbers are never typed. NPxxY is the one
    motif that is a span of generic numbers rather than a segment, because that
    is what it is: 7.49 to 7.53.
    """
    import gemmi

    from flexappeal.refstates import (Catalogue, ensembles, receptor_references,
                                      receptor_residues, refs)
    from flexappeal.refstates.gpcr import NPXXY, build_numbering
    from flexappeal.refstates.traj import core_agreement, gpcr_core

    catalogue = Catalogue(GPCRDB_CACHE)
    numbering = build_numbering(receptor_residues(catalogue))

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
    ensemble_rows = receptor_references(catalogue)

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
        # `receptor_references` returns Ensemble objects, not lists, and the
        # GPCRdb row names the entry `pdb_code` where KLIFS uses `pdb`.
        for entry in ensemble_rows[which].entries:
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
        log(f"{which}: {len(entries)} references")

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
    alignment = ensembles.align_to_common_frame(built)
    log(f"aligned {alignment['n_aligned']} references onto "
        f"{alignment['target']}: core residual mean "
        f"{alignment['mean_core_rmsd_a']} A, max {alignment['max_core_rmsd_a']} A")

    cores = {code: motifs[ensembles.CORE_KEY]
             for state in built for code, motifs in built[state].by_code.items()}
    rigidity = core_agreement(cores)
    log(f"core agreement across {rigidity.get('n')} references: mean "
        f"{rigidity.get('mean_rmsd_a')} A, worst {rigidity.get('max_rmsd_a')} A")
    return (motif_residues, core, built, numbering,
            {"rigidity": rigidity, "alignment": alignment})


def secondary_structure(loaded, motif_residues: dict[str, list[int]]) -> dict:
    """Per-frame DSSP, and each motif's helix and strand content (section 5.6).

    Stored as the per-motif summary plus the whole-chain per-frame strings. The
    strings are what the DSSP ribbon is drawn from and they are the bulk of it,
    so they go in the run's JSON rather than into SQLite: section 7 puts
    per-residue secondary structure outside the database for exactly this
    reason, being a residues-by-frames matrix rather than a set of rows.

    `simplified=False` so the eight DSSP classes survive. The three-class form
    is derived where it is wanted; going the other way is impossible, and the
    difference between a 3-10 helix and an alpha helix is the kind of local
    change this is looking for.
    """
    from flexappeal.refstates import sse

    try:
        rows = sse.dssp_rows(loaded, simplified=False)
    except Exception as exc:                                    # noqa: BLE001
        log(f"DSSP unavailable: {str(exc)[:80]}")
        return {}
    if rows is None or len(rows) == 0:
        return {}

    protein = [r for r in loaded.topology.residues if r.is_protein]
    numbers = [int(r.resSeq) for r in protein]
    index_of = {number: position for position, number in enumerate(numbers)}

    motifs: dict[str, dict] = {}
    for name, residues in motif_residues.items():
        indices = [index_of[n] for n in residues if n in index_of]
        if not indices:
            continue
        helix, strand = [], []
        for frame in range(len(rows)):
            found = sse.motif_sse(rows[frame], name, indices,
                                  [n for n in residues if n in index_of])
            if found.helix_fraction is not None:
                helix.append(found.helix_fraction)
                strand.append(found.strand_fraction)
        if not helix:
            continue
        motifs[name] = {
            "n_residues": len(indices),
            "helix_fraction_mean": sum(helix) / len(helix),
            "helix_fraction_first": helix[0],
            "helix_fraction_last": helix[-1],
            "strand_fraction_mean": sum(strand) / len(strand),
        }

    return {
        "n_frames": int(len(rows)),
        "n_residues": len(numbers),
        "residues": numbers,
        # One character per residue per frame, which is only safe because every
        # DSSP label renders to exactly one character; mdtraj writes an
        # unassigned residue as the two-character "NA".
        "frames": ["".join(sse.UNASSIGNED_CHAR if c == sse.UNASSIGNED
                           else (c if str(c).strip() else "C")
                           for c in row) for row in rows],
        "motifs": motifs,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    parser.add_argument("--start", type=Path, required=True)
    parser.add_argument("--kind", choices=("kinase", "gpcr"), required=True)
    parser.add_argument("--reference-pdb", default="2GQG",
                        help="the entry whose pocket map resolves the "
                             "construct's own motif residues")
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "analysis")
    args = parser.parse_args()

    run_dir = args.run_dir or (ROOT / "data" / "md" / args.name)
    trajectory = run_dir / "trajectory.dcd"
    solvated = run_dir / "solvated.pdb"
    if not trajectory.is_file():
        raise SystemExit(f"no trajectory at {trajectory}")

    from flexappeal.refstates import ensembles, traj

    log(f"loading {trajectory.name} (stride {args.stride})")
    loaded = traj.load(solvated, trajectory, stride=args.stride)
    mapping = traj.deposited_numbering(args.start, loaded.topology)
    changed = traj.renumber(loaded.topology, mapping)
    log(f"{loaded.n_frames} frames, {changed} residues restored to deposited "
        f"numbering ({min(mapping.values())} to {max(mapping.values())})")

    if args.kind == "kinase":
        (motif_residues, core, built, own,
         agreement) = kinase_setup(args.start, args.reference_pdb)
    else:
        (motif_residues, core, built, own,
         agreement) = gpcr_setup(args.start)
    log(f"core: {len(core.residues)} residues, {core.source}")
    for name, residues in sorted(motif_residues.items()):
        log(f"  motif {name}: {len(residues)} residues {residues[:4]}"
            f"{'...' if len(residues) > 4 else ''}")

    if args.kind == "kinase":
        from flexappeal.refstates.pocket import r_spine, r_spine_state

        def rule_labels(chain, index):
            """The rule-based half of section 5.6, per frame.

            Only the measures one frame answers by itself. The DFG call needs
            the vendored comparison frame, which is built against a reference
            rather than per frame, so the spine is what reaches the score here.
            """
            state = r_spine_state(chain, r_spine(own, chain))
            return {"xDFG": (state["state"], "r_spine")}
    else:
        from flexappeal.refstates.gpcr import ionic_lock

        def rule_labels(chain, index):
            """The ionic lock, which is the receptor's own single-frame rule.

            Attached to TM6 because that is the helix whose motion the lock
            reports on. For beta2AR it reads broken in both reference states,
            so in practice it corroborates nothing here and the geometry
            stands alone; see FINDINGS.md.
            """
            return {"TM6": (ionic_lock(chain, own)["state"], "ionic_lock")}

    dssp_payload = secondary_structure(loaded, motif_residues)
    if dssp_payload:
        log(f"DSSP: {dssp_payload['n_residues']} residues over "
            f"{dssp_payload['n_frames']} frames; "
            + ", ".join(f"{m} {v['helix_fraction_mean']:.0%} helix"
                        for m, v in sorted(dssp_payload["motifs"].items())))

    run = traj.analyse(loaded, name=args.name, kind=args.kind, core=core,
                       motif_residues=motif_residues,
                       reference_motifs=ensembles.for_scoring(built),
                       rule_labels=rule_labels)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = run.as_row()
    if dssp_payload:
        payload["dssp"] = dssp_payload
    payload["ensembles"] = {state: e.as_row() for state, e in built.items()}
    payload["core_agreement"] = agreement
    path = out_dir / f"{args.name}.json"
    path.write_text(json.dumps(payload, indent=2))

    last = run.frames[-1] if run.frames else None
    if last is not None:
        log(f"frame {last.time_ns:.1f} ns: whole score "
            f"{last.whole_s if last.whole_s is None else round(last.whole_s, 3)} "
            f"({last.whole_label}), {last.n_conflicts} conflicts")
        for motif in last.motifs:
            if motif.s is None:
                continue
            log(f"   {motif.motif:8s} s {motif.s:+.3f}  {motif.label:12s} "
                f"d_act {motif.d_act_a:.2f} ({motif.nearest_active}) "
                f"d_inact {motif.d_inact_a:.2f} ({motif.nearest_inactive})")
    log(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
