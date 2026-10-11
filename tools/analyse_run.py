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

# The interval run_md.py saves frames at. One definition, because a mismatch
# would silently rescale every time axis in the analysis.
FRAME_PS = 100.0

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
    from flexappeal.refstates.traj import core_agreement, scaffold_core

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
                                              "sse_profiles": sse_profiles}


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
                          "values": npxxy}})


# Section 6.3: coordinates for the web at a 1 ns stride, full resolution kept
# in data/md. A 50 ns run is then 500 frames of protein and ligand only, which
# is a few megabytes rather than the gigabyte the solvated DCD is.
VIEWER_STRIDE_NS = 1.0


def periodic_shifts(traj):
    """Whole-box translations that keep a molecule in one periodic image.

    A protein sitting on a periodic boundary is written into whichever image
    the integrator last wrapped it to, and that can differ from frame to frame.
    beta2AR does exactly this: its centre of mass alternates between y = 0.5 A
    and y = 129 A in a 131.2 A box, a full box length apart, while the protein
    itself stays whole in every frame.

    Nothing measured within a frame cares, which is why the analysis is
    unaffected and the numbers on the page are right. The viewer does care. It
    fits its camera to the whole trajectory, so a receptor that occupies 54 A
    of a 187 A bounding box is drawn small and off centre, and it teleports
    across the box when the slider is dragged.

    So each frame is shifted by the whole number of box vectors that puts its
    centroid in the same image as the first frame's.

    Only the protein is treated this way, and the first attempt here got that
    wrong. Applying the protein's shift to the bilayer as well made the picture
    worse, not better: the bilayer was already in one piece, so half its frames
    were moved a box length and its y span went from 145 Å to 269. A slab
    spanning the box has a stable centroid whatever image the protein is in,
    which is why it needs the other half of the same idea instead. The lipids
    are imaged around wherever the receptor ended up, which
    `write_lipid_bundle` does with the centres the viewer bundle returns.

    Returns an (n_frames, 3) array in nanometres, to be added to every atom of
    the matching frame. Zeros when there is no unit cell to reason about, and
    zeros for a non-orthorhombic box, where a whole-box translation is not the
    right correction and guessing would be worse than leaving it alone.
    """
    import numpy as np

    if traj.unitcell_lengths is None or traj.unitcell_angles is None:
        return np.zeros((traj.n_frames, 3))
    if not np.allclose(traj.unitcell_angles, 90.0, atol=1e-3):
        return np.zeros((traj.n_frames, 3))
    box = np.asarray(traj.unitcell_lengths)
    centre = traj.xyz.mean(axis=1)
    return -np.round((centre - centre[0]) / box) * box


def write_viewer_bundle(loaded, name: str, out_dir: Path,
                        frame_ps: float = 100.0) -> dict:
    """A topology and a strided trajectory the 3D viewer can load.

    Written beside the analysis rather than served from `data/md`, because the
    solvated trajectory is most of a gigabyte and none of the solvent is drawn.
    The topology is PDB and the coordinates XTC, which is the pair Mol* already
    reads on the analysis tab, so the viewer needs no new code path.
    """
    stride = max(1, int(round(VIEWER_STRIDE_NS * 1000.0 / frame_ps)))
    sliced = loaded[::stride]
    # A copy, because `loaded` is the trajectory every metric is measured on
    # and this moves atoms. The shift is a whole number of box vectors and
    # changes nothing measurable, but a function that reaches back into the
    # analysis trajectory to tidy up a picture would be the wrong shape
    # whatever it did.
    sliced = sliced.slice(range(sliced.n_frames), copy=True)
    shifts = periodic_shifts(sliced)
    sliced.xyz = sliced.xyz + shifts[:, None, :]
    out_dir.mkdir(parents=True, exist_ok=True)
    topology = out_dir / f"{name}.pdb"
    trajectory = out_dir / f"{name}.xtc"
    sliced[0].save_pdb(str(topology))
    sliced.save_xtc(str(trajectory))
    moved = int((abs(shifts) > 1e-6).any(axis=1).sum())
    if moved:
        log(f"viewer bundle: {moved} of {sliced.n_frames} frames were in a "
            f"different periodic image and have been brought back to the first")
    return {
        "topology": topology.name, "trajectory": trajectory.name,
        "n_frames": int(sliced.n_frames), "stride": stride,
        "frame_ns": stride * frame_ps / 1000.0,
        "bytes": topology.stat().st_size + trajectory.stat().st_size,
        # Where the receptor ended up, per frame, for the bilayer to be imaged
        # around. Not a result; the caller takes it out of the payload.
        "centres": sliced.xyz.mean(axis=1),
    }


def write_lipid_bundle(solvated: Path, trajectory: Path, name: str,
                       out_dir: Path, stride: int, centres=None) -> dict:
    """The bilayer phosphorus, for section 6.2's membrane toggle.

    The viewer bundle drops the solvent and the lipids on load, and says so, on
    the grounds that a 190,000-atom membrane system is mostly water and none of
    it is drawn. That reasoning is right for the water and wrong for the
    phosphorus: a POPC bilayer around beta2AR is about 476 lipids, so one atom
    each is 476 atoms against roughly 2,600 for the protein. Section 6.2 asks
    for the membrane as transparent phosphate spheres on a toggle, and the
    cost of being able to draw it is a few per cent.

    Written as its own pair of files rather than merged into the viewer bundle.
    Merging would put non-protein atoms into the topology every downstream
    selection runs against, and that path is currently producing correct
    numbers for a live run; a separate structure the viewer loads beside the
    first cannot affect any measurement.

    Returns an empty dict for a system with no bilayer, which is what the
    template tests before offering the toggle.
    """
    import mdtraj as md

    from flexappeal.refstates import membrane

    first = md.load(str(solvated))
    phosphorus = membrane.phosphorus_indices(first.topology)
    if not phosphorus:
        return {}

    lipids = md.load(str(trajectory), top=str(solvated), stride=stride,
                     atom_indices=phosphorus)
    # Imaged around the receptor rather than shifted with it. Applying the
    # protein's own whole-box shift here was the first thing tried and it made
    # the picture worse: the bilayer was already in one piece, so half its
    # frames moved a box length and its y span went from 145 Å to 269.
    #
    # The minimum-image convention about the receptor's centre is the right
    # operation. Every phosphorus is placed in whichever periodic copy sits
    # nearest the receptor, which is a contiguous patch of bilayer around the
    # protein however the integrator happened to wrap either of them. On the
    # live run it brings the lipids' y centre onto the receptor's and collapses
    # their z span from 120 Å to 62, the second being lipids that had been
    # wrapped through the top of the box.
    if (centres is not None and len(centres) == lipids.n_frames
            and lipids.unitcell_lengths is not None):
        import numpy as np

        box = np.asarray(lipids.unitcell_lengths)[:, None, :]
        anchor = np.asarray(centres)[:, None, :]
        offset = lipids.xyz - anchor
        lipids.xyz = anchor + offset - box * np.round(offset / box)
    elif centres is not None:
        log(f"lipid bundle: {len(centres)} protein frames against "
            f"{lipids.n_frames} lipid frames, so the bilayer was not imaged "
            f"around the receptor and may be drawn away from it")
    topology = out_dir / f"{name}.lipids.pdb"
    coords = out_dir / f"{name}.lipids.xtc"
    lipids[0].save_pdb(str(topology))
    lipids.save_xtc(str(coords))
    size = topology.stat().st_size + coords.stat().st_size
    log(f"lipid bundle: {len(phosphorus)} phosphorus over "
        f"{lipids.n_frames} frames, {size / 1e6:.1f} MB")
    return {"topology": topology.name, "trajectory": coords.name,
            "n_atoms": len(phosphorus), "n_frames": int(lipids.n_frames),
            "bytes": size}


def secondary_structure(loaded, motif_residues: dict[str, list[int]],
                        profiles: dict | None = None) -> dict:
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

        # Section 5.6: per residue, the fraction of each reference set sharing
        # this frame's SSE class. Taken at the last frame, which is the one a
        # reader is asking about; the per-frame series is a motif by residue by
        # frame cube and belongs nowhere near a page.
        if profiles:
            found = sse.motif_agreement(
                sse.motif_sse(rows[-1], name, indices,
                              [n for n in residues if n in index_of]).labels,
                name, profiles)
            if found.get("n_residues"):
                motifs[name]["reference_agreement"] = found

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

    # Section 5.4 and 5.5's per-frame readings that are not motif RMSDs. The
    # reference chain for the helix comparisons is the start model, so a
    # rotation is reported against where the run began rather than against an
    # arbitrary frame.
    import gemmi as _gemmi

    _start = _gemmi.read_structure(str(args.start))
    _start.setup_entities()
    reference_chain = _start[0][0]

    if args.kind == "kinase":
        from flexappeal.refstates import kinase as _K
        from flexappeal.refstates.pocket import r_spine, r_spine_state

        def measure(chain, index, _unused):
            """Section 5.4's readings for one frame.

            The fit is onto the start model on the core, computed here: the
            transform the scoring passes maps into the reference ensemble's
            frame, which is not the frame these comparisons are against.
            """
            fit = _K.core_fit(chain, reference_chain, core.residues)
            alpha_c = own.region("alphaC")
            hinge = own.region("hinge")
            found = {
                "dfg_plane": _K.dfg_plane(chain, own),
                "alphac_salt_bridge": _K.alphac_salt_bridge(chain, own),
                "p_loop": _K.p_loop(chain, own),
                "r_spine": r_spine_state(chain, r_spine(own, chain)),
            }
            if alpha_c:
                found["alphac_displacement"] = _K.helix_displacement(
                    chain, reference_chain, alpha_c, fit)
            if hinge:
                found["hinge_backbone"] = _K.backbone_flags(
                    chain, reference_chain, hinge)
            return found

        def rule_labels(chain, index):
            """The rule-based half of section 5.6, per frame.

            Only the measures one frame answers by itself. The DFG call needs
            the vendored comparison frame, which is built against a reference
            rather than per frame, so the spine is what reaches the score here.
            """
            state = r_spine_state(chain, r_spine(own, chain))
            return {"xDFG": (state["state"], "r_spine")}
    else:
        from flexappeal.refstates import gpcr as _G
        from flexappeal.refstates.gpcr import ionic_lock

        def measure(chain, index, _unused):
            """Section 5.5's readings for one frame.

            Fitted onto the start model on the core, for the reason in the
            kinase branch above.
            """
            from flexappeal.refstates import kinase as _K

            fit = _K.core_fit(chain, reference_chain, core.residues)
            return {
                "micro_switches": _G.micro_switches(chain, own),
                "tm6_rotation": _G.tm6_rotation(chain, reference_chain, own,
                                                fit),
                "connector": _G.connector_rmsd(chain, reference_chain, own,
                                               fit),
                "helix8": _G.helix8_angle(chain, own),
                "kinks": _G.kink_angles(chain, own),
            }

        def rule_labels(chain, index):
            """No corroborating rule for the receptor, and that is a result.

            The obvious candidate is the ionic lock, and attaching it was wrong.
            Measured on the two reference models it reads `broken` in both, at
            9.95 A inactive and 15.47 A active, which is the known beta2AR
            result rather than a fault: the R3.50-E6.30 lock does not form in
            this receptor's inactive crystal forms. Mapping `broken` to active
            then made every inactive frame report a conflict on TM6, which the
            carazolol control showed at s = -1.000 against a nearest inactive
            reference of 0.00 A. A permanent conflict on a correct measurement
            is worse than no corroboration.

            So TM6 is scored on geometry alone here and the lock's distance,
            which does separate the states by 5.5 A, is reported beside it as a
            measurement rather than used as a vote. The same reasoning as D015,
            one step further: a criterion that cannot distinguish the states of
            this protein should not be asked to.
            """
            return {}

    dssp_payload = secondary_structure(loaded, motif_residues,
                                       profiles=agreement.get("sse_profiles"))
    if dssp_payload:
        log(f"DSSP: {dssp_payload['n_residues']} residues over "
            f"{dssp_payload['n_frames']} frames; "
            + ", ".join(f"{m} {v['helix_fraction_mean']:.0%} helix"
                        for m, v in sorted(dssp_payload["motifs"].items())))

    run = traj.analyse(loaded, name=args.name, kind=args.kind, core=core,
                       motif_residues=motif_residues,
                       reference_motifs=ensembles.for_scoring(built),
                       rule_labels=rule_labels, measure=measure)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    viewer = write_viewer_bundle(loaded, args.name, out_dir)
    # Taken out of the payload rather than written to it: it is a numpy array
    # for the lipid bundle's use and has no business in the analysis JSON.
    viewer_centres = viewer.pop("centres", None)
    # The bilayer, if there is one, as its own structure beside the protein.
    viewer["lipids"] = write_lipid_bundle(
        solvated, trajectory, args.name, out_dir,
        stride=max(1, int(round(VIEWER_STRIDE_NS * 1000.0 / 100.0)))
        * max(1, args.stride), centres=viewer_centres)
    # The two structures are stepped by one slider, which drives the protein,
    # so a bilayer with a different number of frames would sit at the wrong
    # one. They agree whenever `--stride` is 1, which is every run here, but a
    # disagreement should say so rather than show a bilayer a nanosecond out.
    if (viewer["lipids"]
            and viewer["lipids"]["n_frames"] != viewer["n_frames"]):
        log(f"  bilayer has {viewer['lipids']['n_frames']} frames against the "
            f"protein's {viewer['n_frames']}; the viewer steps them together, "
            f"so the toggle is withheld")
        viewer["lipids"] = {}
    log(f"viewer bundle: {viewer['n_frames']} frames at "
        f"{viewer['frame_ns']:.1f} ns, {viewer['bytes'] / 1e6:.1f} MB")
    payload = run.as_row()
    payload["viewer"] = viewer

    # Section 5.6 asks for a centred 2 ns window "for display only", with the
    # raw values stored. Both are written: the raw series is in `frames` and
    # the smoothed one beside it, so a plot can use the smooth version without
    # the stored numbers having been altered.
    from flexappeal.refstates.state import SMOOTH_WINDOW_NS, smooth

    whole = [f.whole_s for f in run.frames]
    if any(v is not None for v in whole):
        payload["whole_s_raw"] = [None if v is None else round(v, 4)
                                  for v in whole]
        # Only when the run is longer than the window. A 2 ns centred window on
        # a 1 ns run averages the whole series into one value repeated, which
        # looks like a flat trace rather than like a window too wide to apply.
        span_ns = len(whole) * FRAME_PS / 1000.0
        if span_ns > 2.0 * SMOOTH_WINDOW_NS:
            payload["whole_s_smoothed"] = [
                None if v is None else round(v, 4)
                for v in smooth(whole, frame_ps=FRAME_PS)]
            payload["smoothing_window_ns"] = SMOOTH_WINDOW_NS
        else:
            payload["smoothing_skipped"] = (
                f"run is {span_ns:.1f} ns against a "
                f"{SMOOTH_WINDOW_NS} ns window")

    if dssp_payload:
        payload["dssp"] = dssp_payload
    payload["ensembles"] = {state: e.as_row() for state, e in built.items()}
    payload["report"] = agreement
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
