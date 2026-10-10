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


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def membrane_chunks(run_dir: Path, name: str) -> list[dict]:
    """Section 5.8's per-chunk bilayer health, from the full-system trajectory.

    Returns an empty list rather than raising when the run is not a membrane
    one, so the caller does not have to know which rows have a bilayer.
    """
    import mdtraj as md

    from flexappeal.refstates import membrane

    trajectory = run_dir / "trajectory.dcd"
    topology = run_dir / "solvated.pdb"
    if not trajectory.is_file() or not topology.is_file():
        return []

    loaded = md.load(str(trajectory), top=str(topology))
    if not membrane.phosphorus_indices(loaded.topology):
        return []

    frame_ns = 0.1          # section 5.2 saves every 100 ps
    per_chunk = max(1, int(membrane.CHUNK_NS / frame_ns))
    chunks = []
    for start in range(0, loaded.n_frames, per_chunk):
        piece = loaded[start:start + per_chunk]
        if piece.n_frames == 0:
            continue
        chunk = membrane.analyse_chunk(
            piece, start_ns=start * frame_ns,
            end_ns=(start + piece.n_frames) * frame_ns)
        chunks.append(chunk)
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

        counts = store.ingest_analysis(
            Path(args.db), analysis, system_id=meta["system_id"], run_id=name,
            family=meta["family"], environment=meta["environment"],
            run_fields=fields)
        log(f"{name}: {counts['motifs']} motifs, {counts['references']} "
            f"references, {counts['frames']} frame rows ({fields['status']})")

        if not args.skip_membrane:
            chunks = membrane_chunks(run_dir, name)
            if chunks:
                with store.connect(Path(args.db)) as connection:
                    written = store.put_membrane(connection, name, chunks)
                log(f"  {written} membrane rows")

    log(f"wrote {args.db}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
