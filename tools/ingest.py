#!/usr/bin/env python3
"""Load analysed runs into the explorer's SQLite. Spec sections 6 and 7.

Reads what `tools/analyse_run.py` wrote, plus the run's own `run.json` for the
platform and speed, and writes the one database the app serves from. Idempotent,
so it can be re-run after a longer analysis pass without doubling anything.

Membrane chunks are computed here rather than in the analysis pass because they
need the full-system trajectory, lipids and water included, while the state
analysis deliberately drops the solvent on load.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DEFAULT_DB = ROOT / "data" / "refstates.db"

# Which system each run belongs to, and how to describe it. The run names are
# the queue's, so this is the one place the two vocabularies meet.
SYSTEMS = {
    "abl1_apo": {"system_id": "abl1", "family": "kinase", "environment": "water",
                 "start_state": "DFG-in", "start_model": "2GQG"},
    "abl1_imatinib": {"system_id": "abl1", "family": "kinase", "environment": "water",
                      "start_state": "DFG-out", "start_model": "2HYY"},
    "b2ar_apo": {"system_id": "b2ar", "family": "gpcr", "environment": "POPC",
                 "start_state": "active", "start_model": "4LDE"},
    "b2ar_carazolol": {"system_id": "b2ar", "family": "gpcr", "environment": "POPC",
                       "start_state": "inactive", "start_model": "5D5A"},
}


# Section 7's `systems` table has had `uniprot` from the start and nothing ever
# wrote it. These are not remembered: KLIFS `kinase_information?kinase_ID=392`
# returns `"uniprot": "P00519"` for ABL1, and GPCRdb `protein/adrb2_human/`
# returns `"accession": "P07550"`. They are pinned here rather than fetched on
# every ingest because an accession does not change and the ingest should not
# need the network, but the endpoint that gives each one is named so the value
# can be rechecked.
UNIPROT = {"abl1": "P00519", "b2ar": "P07550"}


def system_facts(run_dir: Path, system_id: str) -> dict:
    """Section 7's `uniprot`, `n_atoms` and `forcefield` for one system.

    Every key here is omitted rather than guessed when it cannot be produced,
    which is why this returns a dict to merge instead of a fixed record: a zero
    atom count reads like a measurement.

    The atom count and the force field come from the build manifest where there
    is one. `tools/build_membrane.py` writes `build.json`; the water builder
    does not, so for those rows the count is taken from the finished run's
    `run.json`, and failing that from the solvated topology, where counting the
    coordinate records is exact even above 99,999 atoms because the record name
    does not change.
    """
    facts: dict = {}
    if system_id in UNIPROT:
        facts["uniprot"] = UNIPROT[system_id]

    for manifest in ("build.json", "run.json"):
        path = run_dir / manifest
        if not path.is_file():
            continue
        try:
            built = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        if built.get("atoms"):
            facts["n_atoms"] = int(built["atoms"])
        loaded_ff = built.get("forcefield")
        if loaded_ff:
            # Both writers emit a list. Joining a string would silently give
            # "a, m, b, e, r", which is the kind of wrong that renders.
            facts["forcefield"] = (", ".join(loaded_ff)
                                   if isinstance(loaded_ff, (list, tuple))
                                   else str(loaded_ff))
        if "n_atoms" in facts and "forcefield" in facts:
            break

    if "n_atoms" not in facts:
        topology = run_dir / "solvated.pdb"
        if topology.is_file():
            with topology.open() as handle:
                facts["n_atoms"] = sum(
                    1 for line in handle
                    if line.startswith(("ATOM  ", "HETATM")))
    return facts


# The phase 1 manifests, which carry each reference's catalogue state beside the
# one this build computed. That comparison is what section 6.2's reference
# browser exists to show, and it is the positive control P1 and P2 measure.
MANIFESTS = {"abl1": ROOT / "data" / "manifests" / "refs_abl1.jsonl",
             "b2ar": ROOT / "data" / "manifests" / "refs_b2ar.jsonl"}


# The TM3-TM6 cleft distance that best separates the two states, fitted
# leave-one-out on this reference set as P3. Used here to give the receptor an
# independent state call: GPCRdb's own state is what splits the ensembles, so
# without this the "agreement" column would be comparing the catalogue with a
# copy of itself and reporting 116 of 116.
GPCR_CLEFT_THRESHOLD_A = 11.39


def our_receptor_state(row: dict) -> str | None:
    """Our own state call for a receptor reference, from its cleft distance.

    Larger means active: TM6 swings outward to open the G-protein site. None
    where the distance could not be measured, which is 16 of the 116 entries,
    because their residue at 3.50 or 6.34 is not the one the catalogue names.
    """
    distance = row.get("tm3_tm6")
    if distance is None:
        return None
    return "active" if float(distance) >= GPCR_CLEFT_THRESHOLD_A else "inactive"


def reference_rows(system_id: str) -> dict[str, dict]:
    """Per PDB code, what the catalogue said and what we computed.

    The two manifests have different shapes because the two catalogues answer
    different questions: KLIFS gives a DFG and an alphaC call per structure that
    can be compared with ours, while GPCRdb gives one state per entry and no
    independent geometry, so for the receptor the catalogue's state is recorded
    and there is nothing of ours to disagree with it.
    """
    path = MANIFESTS.get(system_id)
    if path is None:
        log(f"  {system_id}: no manifest is registered for this system")
        return {}
    if not path.is_file():
        # Said out loud rather than returned empty. The caller only logs its
        # reference count when the manifest produced rows, so a manifest that
        # had moved would drop every reference from the explorer and print
        # nothing at all about it.
        log(f"  {system_id}: no manifest at {path}; no references annotated")
        return {}
    out: dict[str, dict] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        code = str(row.get("pdb", "")).upper()
        if not code:
            continue
        # First wins: a code with several chains appears more than once and the
        # chains agree, so a later row would only overwrite with the same call.
        out.setdefault(code, row)
    return out


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def membrane_chunks(run_dir: Path, name: str,
                    expect_bilayer: bool = False) -> list[dict]:
    """Section 5.8's per-chunk bilayer health, from the full-system trajectory.

    Returns an empty list rather than raising when the run is not a membrane
    one, so the caller does not have to know which rows have a bilayer.

    `expect_bilayer` is what keeps that convenience from hiding a fault. For a
    water run an empty list is the right answer and silence is correct. For a
    POPC run it never is: it means the topology was not beside the trajectory or
    the trajectory holds no lipid, and P6 would then read `pending` on a
    finished membrane row with nothing anywhere saying why. `SYSTEMS` already
    records which environment each row is in, so the caller knows.
    """
    import mdtraj as md

    from flexappeal.refstates import membrane

    trajectory = run_dir / "trajectory.dcd"
    topology = run_dir / "solvated.pdb"
    if not trajectory.is_file():
        if expect_bilayer:
            log(f"  {name}: no trajectory yet, so no bilayer QC")
        return []
    if not topology.is_file():
        if expect_bilayer:
            log(f"  {name}: REFUSING the bilayer QC: {topology.name} is not "
                f"beside the trajectory in {run_dir}, so the full system "
                f"cannot be loaded")
        return []

    # The topology alone, to decide whether there is a bilayer at all. One
    # frame, so this costs nothing next to the trajectory.
    reference = md.load(str(topology))
    if not membrane.phosphorus_indices(reference.topology):
        if expect_bilayer:
            log(f"  {name}: REFUSING the bilayer QC: this row is a membrane "
                f"row but its topology has no lipid phosphorus")
        return []

    frame_ns = 0.1          # section 5.2 saves every 100 ps
    per_chunk = max(1, int(membrane.CHUNK_NS / frame_ns))

    # Streamed a chunk at a time rather than loaded whole. This pass needs the
    # solvent, so it cannot drop atoms the way `traj.load` does, and a finished
    # 50 ns membrane row is 500 frames of 194,032 atoms: about 1.2 GB held at
    # once for data that is used one 10 ns window at a time and then discarded.
    # `iterload` caps it at one chunk, around 230 MB, and the chunks were
    # already independent of each other.
    chunks = []
    seen = 0
    for piece in md.iterload(str(trajectory), top=str(topology),
                             chunk=per_chunk):
        if piece.n_frames == 0:
            continue
        chunk = membrane.analyse_chunk(
            piece, start_ns=seen * frame_ns,
            end_ns=(seen + piece.n_frames) * frame_ns)
        chunks.append(chunk)
        seen += piece.n_frames
    gate = membrane.gate_six(chunks)
    log(f"  {len(chunks)} membrane chunks; G6 {'OPEN' if gate['open'] else 'closed'}: "
        f"{gate.get('reason')}")
    return [c.as_row() for c in chunks]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", action="append", required=True,
                        help="a run name, repeatable")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--analysis", type=Path,
                        default=ROOT / "data" / "analysis")
    parser.add_argument("--md", type=Path, default=ROOT / "data" / "md")
    parser.add_argument("--skip-membrane", action="store_true",
                        help="skip the bilayer pass, which reloads the full "
                             "solvated trajectory and is the slow part")
    args = parser.parse_args()

    from flexappeal.refstates import store

    for name in args.name:
        meta = SYSTEMS.get(name)
        if meta is None:
            log(f"{name}: not a known run name; skipping")
            continue
        payload_path = Path(args.analysis) / f"{name}.json"
        if not payload_path.is_file():
            log(f"{name}: no analysis at {payload_path}; skipping")
            continue
        analysis = json.loads(payload_path.read_text())

        run_dir = Path(args.md) / name
        fields = {"start_state": meta["start_state"],
                  "start_model": meta["start_model"], "tier": 2, "replica": 1}
        run_json = run_dir / "run.json"
        if run_json.is_file():
            try:
                done = json.loads(run_json.read_text())
                fields["platform"] = done.get("platform")
                fields["ns_per_day"] = done.get("ns_per_day")
                fields["status"] = done.get("status", "complete")
            except json.JSONDecodeError:
                pass
        else:
            # A run still in flight is recorded as such rather than left out:
            # section 2.3's streaming analysis means the explorer shows partial
            # runs, and a row with no status would look complete.
            fields["status"] = "running"

        manifest = reference_rows(meta["system_id"])
        facts = system_facts(run_dir, meta["system_id"])
        counts = store.ingest_analysis(
            Path(args.db), analysis, system_id=meta["system_id"], run_id=name,
            family=meta["family"], environment=meta["environment"],
            run_fields=fields, system_fields=facts)
        if facts.get("n_atoms"):
            log(f"  system {meta['system_id']}: {facts['n_atoms']:,} atoms"
                + (f", {facts['forcefield']}" if facts.get("forcefield") else "")
                + (f", UniProt {facts['uniprot']}" if facts.get("uniprot") else ""))
        log(f"{name}: {counts['motifs']} motifs, {counts['references']} "
            f"references, {counts['frames']} frame rows ({fields['status']})")

        # Section 5.6's second map axis, computed during the analysis and
        # merged into the manifest rows so both axes live on the same record.
        map_axis = ((analysis.get("report") or {}).get("map_axis")
                    or (analysis.get("core_agreement") or {}).get("map_axis")
                    or {})
        axis_values = map_axis.get("values") or {}
        if axis_values:
            log(f"  map axis {map_axis.get('motif')} to nearest "
                f"{map_axis.get('against')}: {len(axis_values)} references")

        if manifest:
            agreed = annotated = 0
            with store.connect(Path(args.db)) as connection:
                for code, row in manifest.items():
                    theirs = row.get("their_dfg") or row.get("state")
                    ours = (row.get("our_dfg") if "our_dfg" in row
                            else our_receptor_state(row))
                    metrics = {k: v for k, v in row.items()
                               if k not in ("pdb", "chain", "note", "status")}
                    if code in axis_values:
                        metrics["npxxy_to_inactive"] = axis_values[code]
                    store.put_reference(
                        connection, f"{meta['system_id']}:{code}",
                        meta["system_id"], pdb=code,
                        chain=row.get("chain") or "",
                        source="KLIFS" if meta["family"] == "kinase" else "GPCRdb",
                        source_state=theirs, our_state=ours, metrics=metrics,
                        # Omitted rather than nulled when the manifest predates
                        # the field, so an old manifest ingests without writing
                        # a stand-in over a column it knows nothing about.
                        **({"resolution": row["resolution"]}
                           if row.get("resolution") is not None else {}))
                    annotated += 1
                    # Agreement only where both sides have an opinion. A row
                    # this build could not label is not a disagreement and is
                    # not an agreement either.
                    if ours and theirs and str(ours).lower() == str(theirs).lower():
                        agreed += 1
            comparable = sum(
                1 for row in manifest.values()
                if (row.get("our_dfg") if "our_dfg" in row
                    else our_receptor_state(row))
                and (row.get("their_dfg") or row.get("state")))
            log(f"  {annotated} references annotated from the manifest, "
                f"{agreed} of {comparable} agreeing with the catalogue")

        if not args.skip_membrane:
            # The row's environment is already known here, so the bilayer pass
            # can be told whether an empty result is expected or a fault.
            expect_bilayer = meta["environment"] != "water"
            chunks = membrane_chunks(run_dir, name,
                                     expect_bilayer=expect_bilayer)
            if chunks:
                with store.connect(Path(args.db)) as connection:
                    written = store.put_membrane(connection, name, chunks)
                log(f"  {written} membrane rows")
            elif expect_bilayer:
                log(f"  {name}: no membrane rows written for a "
                    f"{meta['environment']} row; P6 will read as pending")

    log(f"wrote {args.db}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
