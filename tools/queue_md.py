#!/usr/bin/env python3
"""Run the reference-state MD queue, one system at a time.

There is one GPU, and running two systems on it at once helps neither: the
beta2AR smoke test and the ABL1 production run shared it for fifteen minutes and
the smoke test got 8% of a core and never reached its first frame, while ABL1
slowed from 213 to 205 ns/day. So the queue is strictly sequential.

Order is section 5.3's priority order, as far as D009's two pairs allow it:
beta2AR active-start first, because its relaxation is the clearest expected
result, then the rest. A row already carrying `run.json` is complete and is
skipped, and `tools/run_md.py` resumes from its own checkpoint, so this can be
killed and restarted at any point. That matters more than it sounds: four rows at
50 ns is most of a day, and the machine sleeps.

Tier 1 first. Section 5.3 requires 1 ns of every row before any long run starts,
and the membrane systems have never been stepped: `--smoke` runs that and stops,
which is also how the membrane ns/day gets measured before the budget is
committed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Section 5.3's priority order. Each row is what `run_md.py` needs to run it.
QUEUE = [
    {"name": "b2ar_apo", "start": "b2ar_apo_start.pdb", "membrane": True,
     "note": "active, apo: the clearest expected relaxation"},
    {"name": "abl1_apo", "start": "abl1_apo_start.pdb", "membrane": False,
     "note": "DFG-in, apo"},
    {"name": "b2ar_carazolol", "start": "b2ar_carazolol_start.pdb", "membrane": True,
     "note": "inactive, carazolol bound"},
    {"name": "abl1_imatinib", "start": "abl1_imatinib_start.pdb", "membrane": False,
     "note": "DFG-out, imatinib bound"},
]

SMOKE_NS = 1.0
SMOKE_EQUIL_NS = 0.2
PRODUCTION_NS = 50.0


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def complete(out_dir: Path, name: str, target_ns: float) -> bool:
    """Whether this row has already reached `target_ns`.

    Read from `run.json`, which `run_md.py` writes only after production
    finishes, rather than from the trajectory length: a trajectory that is being
    appended to right now has a length and is not complete.
    """
    path = out_dir / name / "run.json"
    if not path.is_file():
        return False
    try:
        done = json.loads(path.read_text())
    except json.JSONDecodeError:
        return False
    return float(done.get("production_ns", 0.0)) >= target_ns


def run_row(row: dict, *, out_dir: Path, production_ns: float,
            equil_ns: float | None, smoke: bool) -> dict:
    """One row, as a subprocess so a crash cannot take the queue with it."""
    name = row["name"] + ("_smoke" if smoke else "")
    command = [
        "pixi", "run", "python", str(ROOT / "tools" / "run_md.py"),
        "--name", name,
        "--start", str(ROOT / "data" / "systems" / row["start"]),
        "--production-ns", str(production_ns),
        "--out", str(out_dir),
    ]
    if row["membrane"]:
        command.append("--membrane")
        # Built once by tools/build_membrane.py, under the row's own name. The
        # smoke run writes elsewhere and must still use that same system.
        command += ["--system-dir", str(out_dir / row["name"])]
    if equil_ns is not None:
        command += ["--equil-ns", str(equil_ns)]

    # Streamed line by line to this queue's own log rather than captured and
    # printed at the end. A captured 24-hour run shows nothing at all until it
    # finishes, which is the one time progress is worth having: there is no way
    # to tell a run that is slow from one that has hung.
    row_log = out_dir / f"{name}.log"
    log(f"starting {name}: {row['note']}  (log: {row_log})")
    started = time.time()
    with row_log.open("w") as handle:
        process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True,
                                   bufsize=1)
        for line in process.stdout:
            handle.write(line)
            handle.flush()
            # Echoed selectively: the queue log is the one a human reads, and
            # a full StateDataReporter stream would bury the stage changes.
            if any(key in line for key in
                   ("minimis", "heating", "NPT", "production", "platform",
                    "released", "resumed", "done in", "Error", "error",
                    "Traceback", "NaN", "reusing", "membrane protocol")):
                print(f"    {line.rstrip()}", flush=True)
        returncode = process.wait()
    elapsed = time.time() - started

    if returncode != 0:
        log(f"{name} FAILED after {elapsed / 60:.1f} min; see {row_log}")
        for line in row_log.read_text().strip().splitlines()[-8:]:
            print(f"    {line}", flush=True)
    else:
        log(f"{name} finished in {elapsed / 60:.1f} min")

    speed = None
    run_json = out_dir / name / "run.json"
    if run_json.is_file():
        try:
            speed = json.loads(run_json.read_text()).get("ns_per_day")
        except json.JSONDecodeError:
            pass
    return {"name": name, "ok": returncode == 0,
            "minutes": round(elapsed / 60, 1), "ns_per_day": speed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true",
                        help="Tier 1: 1 ns of every row, then stop")
    parser.add_argument("--only", default=None,
                        help="run one named row rather than the whole queue")
    parser.add_argument("--production-ns", type=float, default=PRODUCTION_NS)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "md")
    parser.add_argument("--wait-for", default=None,
                        help="a run name to wait for before starting, so this "
                             "can be queued behind a run already in flight")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.wait_for:
        waiting = out_dir / args.wait_for / "run.json"
        log(f"waiting for {args.wait_for} to finish before taking the GPU")
        while not waiting.is_file():
            time.sleep(60)
        log(f"{args.wait_for} is done; starting")

    rows = [r for r in QUEUE if args.only is None or r["name"] == args.only]
    if not rows:
        log(f"no row named {args.only}")
        return 1

    target = SMOKE_NS if args.smoke else args.production_ns
    results = []
    for row in rows:
        name = row["name"] + ("_smoke" if args.smoke else "")
        if complete(out_dir, name, target):
            log(f"{name} already at {target} ns; skipping")
            continue
        results.append(run_row(
            row, out_dir=out_dir, production_ns=target,
            equil_ns=SMOKE_EQUIL_NS if args.smoke else None, smoke=args.smoke))

    log("queue done")
    for result in results:
        speed = f"{result['ns_per_day']} ns/day" if result["ns_per_day"] else "no speed recorded"
        log(f"  {result['name']}: {'ok' if result['ok'] else 'FAILED'}, "
            f"{result['minutes']} min, {speed}")

    # The budget check section 0.2 ties to G1, done on measured speed rather
    # than on a guess scaled from the water system's atom count.
    measured = [r for r in results if r["ns_per_day"]]
    if measured and args.smoke:
        total = 0.0
        for result in measured:
            base = result["name"].removesuffix("_smoke")
            total += args.production_ns / result["ns_per_day"] * 24.0
        log(f"projected {args.production_ns} ns per row on measured speed: "
            f"{total:.1f} GPU-hours for {len(measured)} rows")
    return 0 if all(r["ok"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
