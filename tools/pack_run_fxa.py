#!/usr/bin/env python3
"""Pack a reference-state MD row as a `.fxa` results file for the Runs tab.

The reference-state rows were analysed by `tools/analyse_run.py`, which measures
motifs against reference structures and writes `data/analysis/<row>.json`. That
is a different question from the one the Analysis tab answers, and the rows were
never put through the standard pass at all: there was no `.fxa` for any of them,
so they could only be read on pages of their own.

This builds the standard answer for each row, and it does it by **running the
bundle's own analysis script** rather than by reimplementing it. For each row:

  1. a configuration is assembled that describes how that row was actually run,
  2. `bundle.build` renders a real self-extracting bundle from it, exactly as the
     Prepare tab would for a user,
  3. the bundle is unpacked into a work directory, the row's finished output is
     placed where `analyse.py` looks for it, and `analyse.py` is run,
  4. the `.fxa` it writes gains the reference-state payload as extra members and
     lands in `web_runs/<row>/`.

Step 2 and 3 are the point. The Runs tab and the Example tab render whatever
`plots.build_all` finds in a results file, so the only way to be sure these four
rows show the same panels, computed the same way, as the committed lysozyme run
is for the same script to have computed them. A reimplementation here would be a
second definition of every metric, free to drift from the first.

What the rows have that a bundle run does not: the trajectory is a DCD written by
`tools/run_md.py` rather than an XTC written by `run.py`, and the topology is
`solvated.pdb` rather than `topology.pdb`. Both are handled by staging, not by
changing the script, so nothing here is special-cased inside the thing being
reused.

Usage:

    pixi run python tools/pack_run_fxa.py            # every row that needs one
    pixi run python tools/pack_run_fxa.py abl1_apo   # one row
    pixi run python tools/pack_run_fxa.py --force    # rebuild rows already done

Resumable, per rule 5: a row whose `web_runs/<row>/results.fxa` is newer than its
trajectory is skipped unless `--force`, and every attempt appends to
`data/manifests/fxa.jsonl` whether it succeeded or not.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sqlite3
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from flexappeal import bundle, options, schema  # noqa: E402

MD_ROOT = ROOT / "data" / "md"
ANALYSIS_ROOT = ROOT / "data" / "analysis"
SYSTEMS_ROOT = ROOT / "data" / "systems"
RUNS_ROOT = ROOT / "web_runs"
DB = ROOT / "data" / "refstates.db"
MANIFEST = ROOT / "data" / "manifests" / "fxa.jsonl"
WORK = ROOT / "data" / "fxa_work"


# How each row was actually run, taken from its own `run.json` where that holds
# the fact and stated here where it does not. Nothing in this table is a default
# that could be wrong unnoticed: a value that disagrees with the row's own
# record raises rather than being quietly preferred, in `_config`.
#
# `extra_metrics` is what the system earns beyond the nine the committed example
# computes. A membrane row gets the three bilayer measures and a row with a
# ligand gets the two ligand measures, so each page shows the panels its own
# chemistry supports and no others.
ROWS: dict[str, dict[str, Any]] = {
    "abl1_apo": {
        "label": "ABL1 kinase, active start (DFG-in)",
        "start": "abl1_apo_start.pdb",
        "pdb": "2GQG",
        "ligands": [],
        "membrane": False,
        "extra_metrics": [],
        "description": (
            "ABL1 kinase domain from 2GQG, the active DFG-in conformation, "
            "simulated apo in explicit water. One of four reference-state rows: "
            "each protein is run from both of its known states so a motif's "
            "drift can be read as a distance from experiment."),
    },
    "abl1_imatinib": {
        "label": "ABL1 kinase, inactive start with imatinib",
        "start": "abl1_imatinib_start.pdb",
        "pdb": "2HYY",
        "ligands": ["STI"],
        "membrane": False,
        "extra_metrics": ["ligand_rmsd", "ligand_contacts"],
        "description": (
            "ABL1 kinase domain in its inactive DFG-out conformation with "
            "imatinib bound, in explicit water. The positive control for an "
            "inactive start staying inactive, and for a drug holding a state."),
    },
    "b2ar_apo": {
        "label": "beta2 adrenergic receptor, active start",
        "start": "b2ar_apo_start.pdb",
        "pdb": "4LDE",
        "ligands": [],
        "membrane": True,
        "extra_metrics": ["membrane_apl", "membrane_scd", "membrane_thickness"],
        "description": (
            "The beta2 adrenergic receptor from 4LDE, the active state, "
            "embedded in a POPC bilayer and simulated apo. The receptor "
            "relaxes away from a ligand-bound crystal pose, which is the apo "
            "system behaving correctly rather than a method failing."),
    },
    "b2ar_carazolol": {
        "label": "beta2 adrenergic receptor, inactive start with carazolol",
        "start": "b2ar_carazolol_start.pdb",
        "pdb": "2RH1",
        "ligands": ["CAU"],
        "membrane": True,
        "extra_metrics": ["ligand_rmsd", "ligand_contacts",
                          "membrane_apl", "membrane_scd", "membrane_thickness"],
        "description": (
            "The beta2 adrenergic receptor from 2RH1, the inactive state, with "
            "carazolol bound, in a POPC bilayer. The second of the two "
            "inactive-start rows, and the one that answers P5 on a membrane."),
    },
}

# The nine the committed example computes. Named here rather than left to the
# registry's default so that a change to that default cannot silently give these
# four rows a different panel set from the run they are meant to match.
EXAMPLE_METRICS = ["rmsd", "rmsf", "rgyr", "sasa", "dssp",
                   "hbonds", "contacts", "pca", "clusters"]


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def record(row: str, **fields: Any) -> None:
    """Append one line to the stage manifest, successful or not."""
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(
            {"row": row, "at": datetime.now(timezone.utc).isoformat(), **fields}
        ) + "\n")


def _run_record(row: str) -> dict[str, Any]:
    return json.loads((MD_ROOT / row / "run.json").read_text(encoding="utf-8"))


def _config(row: str, meta: dict[str, Any]) -> dict[str, Any]:
    """The configuration that describes how this row was run.

    Read from the row's own `run.json` wherever that records the fact, so the
    manifest inside the results file says what happened rather than what this
    table guessed. The force field, the temperature, the duration and whether
    there is a membrane are all facts the queue wrote down at the time.
    """
    record_ = _run_record(row)
    forcefield = list(record_.get("forcefield") or [])
    if not forcefield:
        raise SystemExit(f"{row}: run.json records no force field")
    if bool(record_.get("membrane")) != meta["membrane"]:
        raise SystemExit(
            f"{row}: run.json says membrane={record_.get('membrane')} and this "
            f"table says {meta['membrane']}; one of them is wrong")

    cfg = options.defaults()
    cfg.update({
        "job_name": row,
        "job_author": "Marc C. Deller",
        "job_description": meta["description"],
        # The protein force field is the first file; the water model is read off
        # the second rather than assumed, because `amber14/tip3p.xml` names it.
        "protein_ff": forcefield[0],
        "water_model": ("tip3p" if len(forcefield) < 2
                        else forcefield[1].rsplit("/", 1)[-1].removesuffix(".xml")),
        "temperature": float(record_.get("temperature_k") or 300.0),
        "production_duration": float(record_.get("production_ns") or 50.0),
        "payload_tier": "standard",
        "analysis_metrics": EXAMPLE_METRICS + list(meta["extra_metrics"]),
        "analysis_selection": "backbone",
        "analysis_reference": "first",
        "traj_selection": "protein_ligand",
        "use_membrane": meta["membrane"],
        "has_ligands": bool(meta["ligands"]),
        "keep_heteroatoms": list(meta["ligands"]),
    })
    if meta["membrane"]:
        # Section 5.2's protocol, and the barostat the run actually used: a
        # membrane system with an isotropic barostat is a different simulation
        # and the validator says so.
        cfg["lipid_type"] = "POPC"
        cfg["barostat"] = "MonteCarloMembrane"

    result = schema.validate(cfg)
    if not result.ok:
        raise SystemExit(f"{row}: " + "; ".join(
            f"{issue.option_id}: {issue.message}" for issue in result.errors))
    return result.config


def _frame_interval_ps(payload: dict[str, Any]) -> float:
    """How far apart this row's frames are, measured from its own analysis.

    The DCD on disk does not know. `tools/run_md.py` writes a frame every 100 ps
    and its own comment says "frame time here is index times the interval": the
    header MDTraj reads back says one time unit per frame, so a 50 ns
    trajectory loads as 0 to 499 ps. The first results file packed from it
    reported a duration of 0.5 ns on every panel.

    Rather than restate 100 ps here, where it could drift from the runner, the
    interval is read off the time axis the reference-state analysis already
    computed for this row, and checked for being uniform. That axis is the one
    every published number for these rows was measured against.
    """
    times = [float(frame["time_ns"]) for frame in payload.get("frames", [])]
    if len(times) < 2:
        raise SystemExit("the analysis file has fewer than two frames")
    steps = {round(b - a, 6) for a, b in zip(times, times[1:])}
    if len(steps) != 1:
        raise SystemExit(f"the analysis time axis is not uniform: {sorted(steps)[:4]}")
    interval_ns = steps.pop()
    if interval_ns <= 0:
        raise SystemExit(f"the analysis time axis does not advance: {interval_ns}")
    return interval_ns * 1000.0


def _stage_trajectory(source: Path, destination: Path, interval_ps: float,
                      timestep_fs: float, topology: Path) -> int:
    """Copy the trajectory with the time axis it was produced on.

    Written as XTC, which stores a time and a step per frame, rather than as
    another DCD, whose header is what lost the axis in the first place. It is
    also what a generated bundle writes, so the analysis script is reading the
    format it was built to read.

    Streamed in chunks because these trajectories are up to 1.2 GB and the
    largest is 194,032 atoms: loading one whole would be a gigabyte of float32
    held to write a gigabyte back out.
    """
    import numpy as np
    import mdtraj as md
    from mdtraj.formats import XTCTrajectoryFile

    steps_per_frame = interval_ps / (timestep_fs / 1000.0)
    written = 0
    with XTCTrajectoryFile(str(destination), "w") as out:
        for chunk in md.iterload(str(source), top=str(topology), chunk=25):
            index = np.arange(chunk.n_frames) + written
            out.write(chunk.xyz,
                      time=index * interval_ps,
                      step=(index * steps_per_frame).astype(np.int64),
                      box=chunk.unitcell_vectors)
            written += chunk.n_frames
    return written


def _merge_ligand_residues(topology: Path, ligands: list[str]) -> None:
    """Put each ligand back into one residue.

    The built systems have every ligand split across two residues, its heavy
    atoms in one and its hydrogens in the other: imatinib is 37 and 31, and
    carazolol 22 and 22. Nothing in the simulation cares, and neither did the
    reference-state analysis, which measures the ligand by selection rather
    than by residue.

    The Analysis tab's binding measurements do care, because they build a
    molecule: with the two halves apart it took "the largest fragment", which
    is the heavy atoms with no hydrogens, and could not infer bond orders from
    them. The consequences were all the quiet kind. No two-dimensional
    depiction, so the mobility inset beside the viewer was missing on every
    bound run while the committed example had one. No strain or torsion
    measurement, so the conformation sub-score stayed pending. And a pose RMSD
    with no symmetry correction, which overstates drift on a ligand with
    equivalent atoms.

    Merged in the staged copy rather than in `data/md`, which is the
    simulation's own output and is not this script's to rewrite.

    Only where the split is the one described: residues of the same name, next
    to each other, one carrying hydrogens and the other not. A system with two
    genuine copies of a ligand has neither property and is left alone.
    """
    import re

    if not ligands or not topology.is_file():
        return
    lines = topology.read_text(encoding="utf-8").splitlines(True)

    for name in ligands:
        blocks: dict[int, dict[str, Any]] = {}
        for index, line in enumerate(lines):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[17:20].strip() != name:
                continue
            number = int(line[22:26])
            block = blocks.setdefault(number, {"rows": [], "hydrogens": 0})
            block["rows"].append(index)
            # The element column where there is one, the atom name where there
            # is not: a PDB from a builder does not always fill columns 77-78.
            element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "",
                                                    line[12:16].strip())[:1]
            if element.upper().startswith("H"):
                block["hydrogens"] += 1
        if len(blocks) != 2:
            continue
        numbers = sorted(blocks)
        if numbers[1] - numbers[0] != 1:
            continue
        carries = [bool(blocks[number]["hydrogens"]) for number in numbers]
        if carries[0] == carries[1]:
            # Both halves have hydrogens or neither does, so this is not the
            # heavy-and-hydrogen split and merging would join two molecules.
            continue
        keep = f"{numbers[0]:>4d}"
        for index in blocks[numbers[1]]["rows"]:
            line = lines[index]
            lines[index] = line[:22] + keep + line[26:]
        log(f"  {name}: {len(blocks[numbers[0]]['rows'])} + "
            f"{len(blocks[numbers[1]]['rows'])} atoms merged into one residue")
    topology.write_text("".join(lines), encoding="utf-8")


def _stage(row: str, cfg: dict[str, Any], work: Path, interval_ps: float,
           meta_ligands: list[str] | None = None) -> Path:
    """Unpack a real bundle and put the row's finished output inside it.

    The bundle arrives as the self-extracting file a user downloads and is
    unpacked here rather than executed: `bootstrap.sh` would build a pixi
    environment of its own, and the point of this step is to run the analysis
    script that bundle carries, in the environment this repository already has.
    """
    start = SYSTEMS_ROOT / ROWS[row]["start"]
    if not start.is_file():
        raise SystemExit(f"{row}: no start structure at {start}")

    built = bundle.build(cfg, start.read_bytes(),
                         structure_filename=f"{ROWS[row]['pdb']}.pdb")
    files = bundle.unpack(built.content)

    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for name, payload in files.items():
        path = work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)

    # Where `analyse.py` looks: OUTPUT/<job_name>*/ holding a topology and a
    # trajectory. The state data is a hard link, because it is small and is
    # only ever read. The topology is a **copy**, because the staging rewrites
    # it: `_merge_ligand_residues` puts the ligand back into one residue, and
    # the first version of this hard-linked the topology and then wrote through
    # the link, which edited `data/md/<row>/solvated.pdb` itself. That file is
    # the simulation's own output. The test that found it is the one asserting
    # the built systems split a ligand in two, which is the fault this staging
    # exists to work around and which that file is the evidence for.
    #
    # The trajectory is rewritten for a different reason: its time axis has to
    # be put back.
    output = (work / cfg["output_dir"]).resolve()
    run_dir = output / row
    run_dir.mkdir(parents=True, exist_ok=True)
    source = MD_ROOT / row
    if not (source / "solvated.pdb").is_file():
        raise SystemExit(f"{row}: {source / 'solvated.pdb'} is missing")
    shutil.copy2(source / "solvated.pdb", run_dir / "topology.pdb")
    state_data = source / "state_data.csv"
    if not state_data.is_file():
        raise SystemExit(f"{row}: {state_data} is missing")
    dest = run_dir / "state_data.csv"
    if dest.exists():
        dest.unlink()
    try:
        dest.hardlink_to(state_data)
    except OSError:
        shutil.copy2(state_data, dest)

    _merge_ligand_residues(run_dir / "topology.pdb", meta_ligands or [])

    log(f"{row}: restaging the trajectory at {interval_ps:.0f} ps per frame")
    frames = _stage_trajectory(
        source / "trajectory.dcd", run_dir / "trajectory.xtc",
        interval_ps, float(cfg["timestep"]), run_dir / "topology.pdb")
    log(f"{row}: {frames} frames staged, "
        f"{frames * interval_ps / 1000:.1f} ns")

    # What `run.py` would have written. The Analysis tab reads the platform, the
    # achieved speed and the wall time off this, and they are facts the queue
    # recorded, so they travel rather than being left blank.
    record_ = _run_record(row)
    wall_hours = record_.get("wall_hours")
    (run_dir / "run_manifest.json").write_text(json.dumps({
        "job_name": row,
        "platform": record_.get("platform", ""),
        "achieved_ns_per_day": record_.get("ns_per_day"),
        "wall_seconds": (float(wall_hours) * 3600
                         if isinstance(wall_hours, (int, float)) else None),
        "produced_ns": record_.get("produced_ns"),
        "equilibration_ns": record_.get("equil_ns"),
        "atoms": record_.get("atoms"),
        "machine": {"host": "Mac Studio", "note": (
            "run by tools/run_md.py under the reference-state build, not by a "
            "generated bundle; this results file was packed from its output")},
    }, indent=2), encoding="utf-8")
    return output


def _analyse(work: Path, output: Path, row: str) -> Path:
    """Run the bundle's own `analyse.py` and return the `.fxa` it wrote."""
    log(f"{row}: running the bundle's analysis pass")
    started = time.time()
    # Run the way a user runs it: `python analyse.py` from the bundle's own
    # directory. Not with `-I`, which was the first attempt: isolated mode does
    # not prepend the script's directory to `sys.path`, and a bundle for a run
    # with a ligand carries `fxalib` beside `analyse.py` and imports it for the
    # version 2 package. It failed with ModuleNotFoundError on the first ligand
    # row. The script and everything beside it were generated by this
    # repository from its own templates a moment earlier, so there is nothing
    # here that isolation would be protecting.
    result = subprocess.run(
        [sys.executable, str(work / "analyse.py")],
        cwd=str(work), capture_output=True, text=True, check=False)
    elapsed = time.time() - started
    tail = "\n".join((result.stdout or "").strip().splitlines()[-25:])
    if result.returncode != 0:
        print(tail)
        print((result.stderr or "").strip()[-2000:], file=sys.stderr)
        raise SystemExit(f"{row}: the analysis pass exited {result.returncode}")
    log(f"{row}: analysed in {elapsed / 60:.1f} min")

    written = sorted(output.glob("*.fxa"))
    # The version 2 package is a different artefact for a different reader, and
    # the Analysis tab refuses it by design, so it is not the one to pick up.
    written = [p for p in written if not p.name.endswith("_package.fxa")]
    if not written:
        print(tail)
        raise SystemExit(f"{row}: the analysis pass wrote no results file")
    return written[0]


def _references(system_id: str) -> list[dict[str, Any]]:
    """The reference ensemble for one system, for the panel that lists it.

    From the explorer database, which is where phase 1 put the catalogue state
    beside the state this build computed. That comparison is the whole content
    of the panel, so a row missing either side is carried as it is rather than
    filled in: the panel says "no opinion" where the catalogue has none.
    """
    if not DB.is_file():
        return []
    rows: list[dict[str, Any]] = []
    with sqlite3.connect(f"file:{DB}?mode=ro", uri=True) as connection:
        connection.row_factory = sqlite3.Row
        try:
            cursor = connection.execute(
                "SELECT * FROM reference_structures WHERE system_id = ? "
                "ORDER BY pdb", (system_id,))
        except sqlite3.DatabaseError:
            return []
        for row in cursor:
            item = dict(row)
            # Parsed rather than carried as a string inside JSON, and kept
            # whole: measured, it is 4.7 KB for the 20 ABL1 references and
            # 11.9 KB for the 116 receptor ones, which is nothing beside a
            # 10 MB results file. It holds the coordinates the state map plots
            # and the per-reference agreement the controls panel reports, so a
            # slimmed copy would be a third definition of which measurements a
            # panel may show.
            try:
                item["metrics"] = json.loads(item.pop("metrics_json") or "{}")
            except (KeyError, json.JSONDecodeError):
                item["metrics"] = {}
            rows.append(item)
    return rows


def _membrane_chunks(row: str) -> list[dict[str, Any]]:
    """The bilayer QC chunks, with their verdicts against section 5.8's windows."""
    if not DB.is_file():
        return []
    with sqlite3.connect(f"file:{DB}?mode=ro", uri=True) as connection:
        connection.row_factory = sqlite3.Row
        try:
            cursor = connection.execute(
                "SELECT * FROM membrane WHERE run_id = ? ORDER BY time_ns", (row,))
        except sqlite3.DatabaseError:
            return []
        return [dict(item) for item in cursor]


def _refstates_payload(row: str) -> dict[str, Any]:
    """The reference-state analysis, as it will travel inside the results file.

    The analysis file is carried whole rather than summarised: every panel this
    feeds was previously drawn from it, and a summary here would be a second
    definition of what each panel may show. Two things are changed. The `viewer`
    block named files in `data/analysis`, which do not exist beside a results
    file, so it is replaced by what the viewer members are actually called. And
    the reference ensemble and the bilayer chunks, which lived in the explorer
    database, are added, because a run in the library is one file.
    """
    path = ANALYSIS_ROOT / f"{row}.json"
    if not path.is_file():
        raise SystemExit(f"{row}: no reference-state analysis at {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    system_id = "abl1" if row.startswith("abl1") else "b2ar"

    viewer = payload.get("viewer") or {}
    overlay = viewer.get("overlay") or {}
    lipids = dict(viewer.get("lipids") or {})
    if lipids and (ANALYSIS_ROOT / f"{row}.lipids.pdb").is_file():
        # Named by member rather than by path. The viewer fetches these through
        # the session route, which serves them out of the results file.
        lipids["topology"] = "lipids.pdb"
        lipids["trajectory"] = "lipids.xtc"
    else:
        lipids = {}
    payload["viewer"] = {
        # Each ghost keeps the facts the label on its toggle is made of: which
        # reference it is and how far this run sits from it. Without them the
        # control reads "ghost the nearest active" and names nothing, and the
        # whole point of the overlay is that a reader can see what they are
        # comparing against.
        "overlay": {
            state: {"member": f"refstates_{state}.pdb",
                    "code": detail.get("code"),
                    "mean_motif_rmsd_a": detail.get("mean_motif_rmsd_a"),
                    "motifs_checked": detail.get("motifs_checked")}
            for state, detail in overlay.items()
            if (ANALYSIS_ROOT / f"{row}.ref_{state}.pdb").is_file()
        },
        "lipids": lipids,
        # Set by the caller once the topology has actually been painted, so the
        # page never claims a colouring the structure cannot support.
        "motif_painted": False,
    }
    payload["system_id"] = system_id
    payload["reference_structures"] = _references(system_id)
    payload["membrane_chunks"] = _membrane_chunks(row)
    payload["packed_at"] = datetime.now(timezone.utc).isoformat()
    return payload


def _translate_residues(topology, payload: dict[str, Any]) -> dict | None:
    """The motif residue numbers, in the numbering the packed topology uses.

    The two do not agree, and the first attempt at painting assumed they did:
    it wrote motif numbers onto 92 of 4,255 atoms instead of about 800, because
    the handful of residues whose numbers happened to coincide were the only
    ones that matched. The reference-state analysis works in the construct's own
    numbering, 233 to 494 for this ABL1 domain, resolved from KLIFS; the
    topology the standard pass writes has been through OpenMM and MDTraj and
    comes out renumbered from 1.

    Mapped by position rather than by an offset, because an offset assumes the
    span has no gaps and nothing here guarantees that. The analysis's own
    residue axis, which the secondary-structure ribbon is drawn on, is the list
    of residue numbers in order; the packed topology's protein residues are the
    same residues in the same order.

    Returns None if the two lengths disagree at all, which leaves the cartoon on
    its default theme. A mapping that is wrong by one residue would paint every
    motif one residue out, and a cartoon coloured slightly wrongly is worse than
    one not coloured: nothing on the page would look broken.
    """
    axis = ((payload.get("dssp") or {}).get("residues")) or []
    if not axis:
        log("  no residue axis in the analysis, so the cartoon is not painted")
        return None
    protein = [residue for residue in topology.residues if residue.is_protein]
    if len(protein) != len(axis):
        log(f"  {len(protein)} protein residues in the packed topology against "
            f"{len(axis)} on the analysis axis, so the cartoon is not painted")
        return None
    renumbered = {int(number): protein[index].resSeq
                  for index, number in enumerate(axis)}
    out: dict[str, list[int]] = {}
    for motif, numbers in (payload.get("motif_residues") or {}).items():
        mapped = [renumbered[int(number)] for number in numbers
                  if int(number) in renumbered]
        if mapped:
            out[motif] = mapped
    return out or None


def _repaint_topology(source: Path, payload: dict[str, Any]) -> bytes | None:
    """The packed topology with each atom's motif number in the B-factor column.

    Section 6.2 asks for the cartoon coloured by the motif tokens, and the
    viewer gets that number out of the structure rather than out of an API
    call: the vendored Mol* bundle exposes the `uncertainty` colour theme,
    which reads B-factors, and `analysis.js` indexes section 6.1's palette
    straight into it.

    Painted here rather than taken from the reference-state viewer bundle,
    which already has it. That bundle is a different selection written by a
    different pass, and pairing it with the trajectory from this one would mean
    trusting two independent atom orderings to agree. This repaints the
    topology the standard pass itself wrote, so the structure and the
    trajectory beside it cannot disagree by construction.

    `motif_bfactors` is imported from the analysis tool rather than reproduced:
    it owns the rule that an overlap belongs to the smaller motif, which is
    what keeps ABL1's four-residue xDFG from being buried inside its
    twenty-six residue activation loop.
    """
    import importlib.util

    import mdtraj as md

    spec = importlib.util.spec_from_file_location(
        "_analyse_run", ROOT / "tools" / "analyse_run.py")
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    with zipfile.ZipFile(source) as archive:
        if "topology.pdb" not in archive.namelist():
            return None
        raw = archive.read("topology.pdb")

    scratch = WORK / "_paint"
    scratch.mkdir(parents=True, exist_ok=True)
    source_pdb = scratch / "topology.pdb"
    source_pdb.write_bytes(raw)
    try:
        structure = md.load(str(source_pdb))
        motif_residues = _translate_residues(
            structure.topology, payload)
        if motif_residues is None:
            return None
        painted = module.motif_bfactors(
            structure.topology, motif_residues,
            payload.get("kind") or "kinase")
        if painted is None:
            return None
        out = scratch / "painted.pdb"
        structure.save_pdb(str(out), bfactors=painted)
        covered = int((painted > 0).sum())
        log(f"  motif numbers on {covered} of {structure.n_atoms} atoms")
        return out.read_bytes()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _frame_transform(mobile, target):
    """The rigid transform that takes one frame's atoms onto another's.

    Both arrays are the same atoms of the same system in the same order, in two
    coordinate frames that differ because they came out of two passes: the
    reference-state analysis images and superposes for its own measurements,
    and the bundle's `analyse.py` images independently before it packs. Nothing
    downstream noticed, because each file was only ever drawn beside itself.
    Put in one viewer together, the ghost sat 50 angstrom off the protein and
    the bilayer somewhere else again.

    Superposed rather than assumed to be a translation. Imaging only ever
    translates whole molecules, so a translation is what this should come out
    as, but a superposition costs nothing and cannot be wrong about it.
    `traj.superpose` is the build's own Kabsch, the one every motif distance is
    measured after, returning the transform so it can be applied to atoms that
    were not part of the fit.

    **Fit on the protein and nothing else.** The first version fitted on every
    atom the two frames share, which includes the ligand, and came out with a
    15.6 angstrom residual on a pair of frames that differ by a translation.
    The ligand is its own molecule, so imaging can and does put it in a
    different periodic copy in one frame from the other: 44 atoms of 5,288, a
    box length out, is enough to drag the fit that far. On the protein alone
    the residual is 0.005 angstrom, which is what a pure translation should
    give and is the check that this transform is the right one.
    """
    from flexappeal.refstates import traj as traj_mod

    rotation, mobile_centre, target_centre = traj_mod.superpose(mobile, target)
    residual = traj_mod.rmsd(
        traj_mod.apply_transform(mobile, rotation, mobile_centre, target_centre),
        target)
    return (rotation, mobile_centre, target_centre), residual


def _packed_frame(produced: Path):
    """The packed topology's first frame, which everything else is aligned to."""
    import mdtraj as md

    scratch = WORK / "_align"
    scratch.mkdir(parents=True, exist_ok=True)
    path = scratch / "packed.pdb"
    with zipfile.ZipFile(produced) as archive:
        path.write_bytes(archive.read("topology.pdb"))
    return md.load(str(path)), scratch


def _orientation(packed, payload: dict[str, Any],
                 phosphorus) -> dict[str, Any]:
    """Which way is up for a membrane run, measured rather than assumed.

    Two facts the viewer needs and neither is a convention it can guess. The
    bilayer normal is the axis the phosphorus atoms are least spread along: a
    leaflet is wide in the plane of the membrane and thin across it, and on
    this receptor the spread is 38 and 40 angstrom in the plane against 19
    across. The camera looks along one of the other two axes, which is what
    puts the bilayer edge on.

    Which end of that axis is extracellular comes from helix 8, which lies on
    the intracellular face: if its centroid is below the protein's along the
    normal then extracellular is the positive direction, and the extracellular
    loops are at the top. Taken from the structure rather than from the OPM
    convention the system was built with, so a system built some other way
    still comes out the right way up.

    Returns an empty dict when either fact cannot be established, and the
    viewer then leaves the camera where Mol* put it.
    """
    import numpy as np

    if phosphorus is None or len(phosphorus) < 20:
        return {}
    spread = np.asarray(phosphorus).std(axis=0)
    normal = int(np.argmin(spread))
    others = [value for index, value in enumerate(spread) if index != normal]
    if min(others) < spread[normal] * 1.5:
        # Not a bilayer lying across one axis, so there is no edge-on view to
        # offer and guessing one would point the camera somewhere arbitrary.
        return {}

    axis = payload.get("dssp", {}).get("residues") or []
    helix8 = (payload.get("motif_residues") or {}).get("H8") or []
    if not axis or not helix8:
        return {}
    protein = [residue for residue in packed.topology.residues
               if residue.is_protein]
    if len(protein) != len(axis):
        return {}
    renumbered = {int(number): protein[index].resSeq
                  for index, number in enumerate(axis)}
    wanted = {renumbered[int(n)] for n in helix8 if int(n) in renumbered}
    indices = [atom.index for atom in packed.topology.atoms
               if atom.residue.is_protein and atom.residue.resSeq in wanted]
    if not indices:
        return {}

    coordinates = packed.xyz[0] * 10.0
    towards = float(coordinates[:, normal].mean()
                    - coordinates[indices][:, normal].mean())
    if abs(towards) < 1.0:
        return {}
    sign = 1.0 if towards > 0 else -1.0
    up = [0.0, 0.0, 0.0]
    up[normal] = sign
    # The camera looks along the first axis that is not the normal, which puts
    # the normal across the screen rather than into it.
    along = [0.0, 0.0, 0.0]
    along[(normal + 1) % 3] = 1.0
    return {"up": up, "along": along, "normal_axis": normal,
             "note": ("the bilayer normal is the axis the phosphorus atoms are "
                      "least spread along, and helix 8 says which end of it is "
                      "intracellular")}


# How far from the receptor a lipid is still drawn. Eighteen angstrom is about
# two lipids out from the protein on every side: enough that the bilayer reads
# as a membrane rather than as a collar, and little enough that the camera
# frames the receptor rather than the box.
LIPID_PATCH_A = 18.0


def _viewer_members(row: str, payload: dict[str, Any], produced: Path,
                    work: Path, cfg: dict[str, Any]) -> dict[str, bytes]:
    """The ghosts and the bilayer, in the coordinate frame the viewer draws in.

    Both used to be copied straight out of `data/analysis`, where they are
    written in the reference-state analysis's own frame. That frame is not the
    one the packed trajectory is in, so in the viewer the ghost was tens of
    angstrom off the protein and the bilayer was somewhere else again. Each is
    transformed onto the packed first frame here.

    The bilayer is the whole lipid rather than its phosphorus. Phosphorus alone
    is 476 atoms and draws as a cloud of dots; the full chains are 24,752 heavy
    atoms, which is 2 MB for one frame and 48 MB for five hundred. So it is one
    frame: the bilayer is context for a receptor that sits in it, the
    measurements of the bilayer are in their own panel, and 48 MB on every page
    view to watch lipid tails wobble is not a trade worth making. Said on the
    control rather than left for a reader to notice.
    """
    import mdtraj as md
    import numpy as np

    members: dict[str, bytes] = {}
    packed, scratch = _packed_frame(produced)
    try:
        # The ghosts, from the reference-state viewer bundle's frame. Its
        # topology is the same atoms in the same order as the packed one, which
        # is checked rather than trusted: a mismatch would superpose two
        # different molecules and put the ghost anywhere at all.
        reference = ANALYSIS_ROOT / f"{row}.pdb"
        if reference.is_file():
            source = md.load(str(reference))
            same = ([(a.name, a.residue.name) for a in source.topology.atoms]
                    == [(a.name, a.residue.name) for a in packed.topology.atoms])
            if same:
                backbone = packed.topology.select("protein")
                transform, residual = _frame_transform(
                    source.xyz[0][backbone] * 10.0,
                    packed.xyz[0][backbone] * 10.0)
                log(f"  ghosts aligned onto the packed frame, {residual:.3f} "
                    f"angstrom residual")
                for state in ("active", "inactive"):
                    ghost = ANALYSIS_ROOT / f"{row}.ref_{state}.pdb"
                    if not ghost.is_file():
                        continue
                    loaded = md.load(str(ghost))
                    from flexappeal.refstates import traj as traj_mod

                    loaded.xyz[0] = (traj_mod.apply_transform(
                        loaded.xyz[0] * 10.0, *transform) / 10.0)
                    out = scratch / f"ref_{state}.pdb"
                    loaded.save_pdb(str(out))
                    members[f"refstates_{state}.pdb"] = out.read_bytes()
            else:
                log("  the reference-state bundle has different atoms from the "
                    "packed topology, so no ghost is offered")

        if not cfg.get("use_membrane"):
            return members

        # The bilayer, from the staged trajectory this run was packed from.
        staged = work / cfg["output_dir"] / row
        topology = staged / "topology.pdb"
        trajectory = staged / "trajectory.xtc"
        if not topology.is_file() or not trajectory.is_file():
            return members
        solvated = md.load(str(topology))
        lipid = solvated.topology.select("resname POP and element != H")
        if not len(lipid):
            return members

        # The protein, in both frames, which is what the transform is fitted
        # on for the reason `_frame_transform` gives.
        shared = solvated.topology.select("protein")
        backbone = packed.topology.select("protein")
        if len(shared) != len(backbone):
            log(f"  {len(shared)} protein atoms in the staged system against "
                f"{len(backbone)} packed, so no bilayer is offered")
            return members

        first = md.load_frame(str(trajectory), 0, top=str(topology))
        transform, residual = _frame_transform(
            first.xyz[0][shared] * 10.0, packed.xyz[0][backbone] * 10.0)

        # The normal is measured on the whole bilayer, before the patch is cut
        # out of it. A leaflet is wide in the plane of the membrane and thin
        # across it, which is what identifies the normal, and a patch 18
        # angstrom around the receptor is about as wide as it is thick: the
        # test correctly refused to call it a bilayer, and the camera was left
        # pointing down the normal at a membrane it could not see.
        whole_phosphorus = np.array(
            [index for index in lipid
             if solvated.topology.atom(index).name == "P"])

        # The patch of bilayer around the receptor rather than the whole box.
        #
        # All 476 lipids is 24,752 atoms, and Mol* frames the camera on
        # everything it holds: the receptor came out a speck in the middle of a
        # membrane four times its width, which is a picture of a box rather
        # than of a protein in a membrane. Restoring the camera after the load
        # does not hold, because the framing runs on Mol*'s own schedule after
        # it. A patch frames itself.
        #
        # Whole residues, chosen once at the first frame: a lipid is in or out,
        # never half drawn, and membership does not change under the reader.
        near = set(md.compute_neighbors(
            first, LIPID_PATCH_A / 10.0, shared, haystack_indices=lipid)[0])
        residues = {solvated.topology.atom(index).residue.index
                    for index in near}
        lipid = np.array([index for index in lipid
                          if solvated.topology.atom(index).residue.index
                          in residues])
        if not len(lipid):
            return members

        # Imaged about the receptor before anything is measured off it, which
        # is the step that makes a contiguous patch of bilayer out of lipids
        # scattered through whichever periodic copy the integrator left them
        # in.
        #
        # Before, not after. The orientation is read off the phosphorus spread,
        # and unimaged that spread is 35, 36 and 38 angstrom, which is a cube
        # and not a bilayer: the test refused to name a normal and the camera
        # was left pointing down it. Imaged, the same atoms are 38, 40 and 19.
        from flexappeal.refstates import traj as traj_mod

        def imaged(indices):
            coordinates = first.xyz[0][indices] * 10.0
            if first.unitcell_lengths is None:
                return coordinates
            box = np.asarray(first.unitcell_lengths)[0] * 10.0
            offset = coordinates - centre
            return centre + offset - box * np.round(offset / box)

        centre = (first.xyz[0][shared] * 10.0).mean(axis=0)
        moved = traj_mod.apply_transform(imaged(lipid), *transform)

        bilayer = first.atom_slice(lipid)
        bilayer.xyz = (moved / 10.0).reshape(1, -1, 3).astype("float32")
        out = scratch / "lipids.pdb"
        bilayer.save_pdb(str(out))
        members["lipids.pdb"] = out.read_bytes()
        log(f"  bilayer: {len(lipid)} lipid atoms aligned onto the packed "
            f"frame, {residual:.3f} angstrom residual, "
            f"{len(members['lipids.pdb']) / 1e6:.1f} MB")

        payload.setdefault("viewer", {})["lipids"] = {
            "topology": "lipids.pdb", "n_atoms": int(len(lipid)),
            "n_residues": int(len({solvated.topology.atom(i).residue.index
                                   for i in lipid})),
            "n_frames": 1, "within_a": LIPID_PATCH_A,
        }
        # Imaged and transformed the same way the patch was, so the axes are
        # the packed frame's. A rigid transform does not change a spread, but
        # the sign of "which way is extracellular" is read against the packed
        # protein, so both have to be in one frame.
        whole = traj_mod.apply_transform(imaged(whole_phosphorus), *transform)
        payload["viewer"]["orientation"] = _orientation(packed, payload, whole)
        log(f"  orientation: {payload['viewer']['orientation'] or 'not established'}"
            f" (phosphorus {len(whole_phosphorus)}, spread "
            f"{np.round(whole.std(axis=0), 1).tolist()})")
        return members
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _extra_members(row: str, payload: dict[str, Any]) -> dict[str, bytes]:
    """The payload itself. The viewer's own members are built by `_viewer_members`."""
    return {"refstates.json": json.dumps(payload).encode("utf-8")}


def _repack(source: Path, destination: Path, extra: dict[str, bytes]) -> dict[str, Any]:
    """Copy a results file, adding members and keeping the manifest honest.

    The manifest carries a size and a SHA-256 for every member and `fxa.load`
    checks them, so a member added without updating that map is a file that
    fails verification. Rewritten rather than appended to: a zip can hold two
    entries of the same name, and the reader would then see whichever came
    first.
    """
    with zipfile.ZipFile(source) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    manifest = json.loads(members.pop("manifest.json"))
    members.update(extra)

    manifest["files"] = {
        name: {"bytes": len(payload),
               "sha256": hashlib.sha256(payload).hexdigest()}
        for name, payload in sorted(members.items())
    }
    # Said in the manifest as well as implied by the members, so a reader that
    # has only the manifest can tell this file carries the extra analysis.
    manifest["refstates"] = {
        "member": "refstates.json",
        "note": ("packed from the reference-state build by "
                 "tools/pack_run_fxa.py; the standard members were written by "
                 "the bundle's own analyse.py"),
    }

    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("manifest.json", json.dumps(manifest, indent=2))
        for name, payload in sorted(members.items()):
            archive.writestr(name, payload)
    return manifest


def pack(row: str, *, force: bool = False) -> bool:
    """One row, end to end. Returns whether it packed anything."""
    meta = ROWS[row]
    stored = RUNS_ROOT / row / "results.fxa"
    trajectory = MD_ROOT / row / "trajectory.dcd"
    if (stored.is_file() and not force
            and stored.stat().st_mtime > trajectory.stat().st_mtime):
        log(f"{row}: already packed and newer than its trajectory, skipping")
        return False

    cfg = _config(row, meta)
    payload = _refstates_payload(row)
    interval_ps = _frame_interval_ps(payload)
    work = WORK / row
    try:
        output = _stage(row, cfg, work, interval_ps, meta["ligands"])
        produced = _analyse(work, output, row)
        extra = _extra_members(row, payload)
        extra.update(_viewer_members(row, payload, produced, work, cfg))
        painted = _repaint_topology(produced, payload)
        if painted:
            # Replaces the member the standard pass wrote, which is why the
            # repack rewrites the archive rather than appending to it.
            extra["topology.pdb"] = painted
        payload.setdefault("viewer", {})["motif_painted"] = bool(painted)
        # Rewritten after the flag is known, rather than before: the page reads
        # it to decide whether to ask for the motif colouring at all.
        extra["refstates.json"] = json.dumps(payload).encode("utf-8")
        manifest = _repack(produced, stored, extra)
    finally:
        # The work directory holds hard links to gigabytes of trajectory and a
        # rendered bundle; neither is worth keeping once the results file is
        # written.
        shutil.rmtree(work, ignore_errors=True)

    (RUNS_ROOT / row / "run.json").write_text(json.dumps({
        "label": meta["label"],
        "added_at": datetime.now(timezone.utc).isoformat(),
        "source": "reference-state build",
    }, indent=2), encoding="utf-8")

    size_mb = stored.stat().st_size / 1e6
    trajectory_info = manifest.get("trajectory", {})
    log(f"{row}: {size_mb:.1f} MB, {trajectory_info.get('frames_packed')} frames "
        f"of {trajectory_info.get('atoms_packed')} atoms, "
        f"{trajectory_info.get('duration_ns', 0):.1f} ns")
    record(row, ok=True, bytes=stored.stat().st_size,
           frames=trajectory_info.get("frames_packed"),
           atoms=trajectory_info.get("atoms_packed"),
           members=sorted(manifest.get("files", {})))
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # Validated below rather than with `choices`, which argparse also applies
    # to the empty list `nargs="*"` produces when no row is named: passing no
    # rows at all, which is the documented way to pack every row, was refused
    # with "invalid choice: []".
    parser.add_argument("rows", nargs="*",
                        help="which rows to pack; every row by default")
    parser.add_argument("--force", action="store_true",
                        help="rebuild a row whose results file is already current")
    args = parser.parse_args()

    rows = args.rows or sorted(ROWS)
    unknown = [row for row in rows if row not in ROWS]
    if unknown:
        parser.error(f"no such row: {', '.join(unknown)}. "
                     f"Known rows: {', '.join(sorted(ROWS))}")
    packed = 0
    for row in rows:
        try:
            packed += bool(pack(row, force=args.force))
        except SystemExit as exc:
            record(row, ok=False, error=str(exc))
            print(f"!! {exc}", file=sys.stderr)
            return 1
    log(f"{packed} of {len(rows)} row(s) packed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
