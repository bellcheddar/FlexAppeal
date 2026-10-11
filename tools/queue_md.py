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
#
# `ligand` is the CCD residue name of a bound ligand, and a row that has one
# cannot be built without it: amber14 has no template for imatinib, so the
# ligand needs a SMIRNOFF one generated from its chemical definition. This
# field was missing and `abl1_imatinib` would have failed at build time, after
# the three rows ahead of it had already spent their GPU hours. The membrane
# rows carry their ligand too, for the refusal below rather than for the build:
# their system is built once by `tools/build_membrane.py`, which does its own
# ligand handling, and `run_md.py` reloads it.
QUEUE = [
    {"name": "b2ar_apo", "start": "b2ar_apo_start.pdb", "membrane": True,
     "ligand": None,
     "note": "active, apo: the clearest expected relaxation"},
    {"name": "abl1_apo", "start": "abl1_apo_start.pdb", "membrane": False,
     "ligand": None,
     "note": "DFG-in, apo"},
    {"name": "b2ar_carazolol", "start": "b2ar_carazolol_start.pdb", "membrane": True,
     "ligand": "CAU",
     "note": "inactive, carazolol bound"},
    {"name": "abl1_imatinib", "start": "abl1_imatinib_start.pdb", "membrane": False,
     "ligand": "STI",
     "note": "DFG-out, imatinib bound"},
]

SMOKE_NS = 1.0
SMOKE_EQUIL_NS = 0.2
PRODUCTION_NS = 50.0
# Matches `run_md.SHORT_ROW_TOLERANCE_NS`. Kept as a literal rather than
# imported because `run_md` builds an OpenMM system at import time on some
# paths, and the queue should not pay that to read one float.
SHORT_ROW_TOLERANCE_NS = 0.05


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def complete(out_dir: Path, name: str, target_ns: float) -> bool:
    """Whether this row has already reached `target_ns`.

    Read from `run.json`, which `run_md.py` writes only after production
    finishes, rather than from the trajectory length: a trajectory that is being
    appended to right now has a length and is not complete.

    Judged on `produced_ns`, where the row's context actually ended, in
    preference to `production_ns`, the target it was asked for. Comparing the
    target against itself made this unfalsifiable: a row that stopped 10 ns
    short still wrote 50 into `production_ns` and was certified by it. A row
    found short returns False and is resumed to its target on the next pass.
    Rows written before `produced_ns` existed carry only the target and are
    taken at their word. Section 4ag.
    """
    path = out_dir / name / "run.json"
    if not path.is_file():
        return False
    try:
        done = json.loads(path.read_text())
    except json.JSONDecodeError:
        return False
    reached = done.get("produced_ns", done.get("production_ns"))
    if reached is None:
        return False
    return float(reached) >= target_ns - SHORT_ROW_TOLERANCE_NS


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
    # A water row with a ligand needs its chemistry, because `run_md.py` builds
    # the system itself and amber14 has no template for a drug. The CCD code and
    # not a path to the SDF: the chemistry comes from the dictionary and the
    # coordinates from the deposition, and `run_md.py` fetches the definition
    # itself if it is not already on disk.
    if row.get("ligand") and not row["membrane"]:
        command += ["--ligand", row["ligand"]]
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
    note("started", name, production_ns=production_ns,
         membrane=bool(row["membrane"]), ligand=row.get("ligand"))
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
        tail = row_log.read_text().strip().splitlines()[-8:]
        for line in tail:
            print(f"    {line}", flush=True)
        record_failure(out_dir / name, name, returncode, elapsed, tail)
        note("failed", name, returncode=returncode,
             minutes=round(elapsed / 60, 1))
    else:
        log(f"{name} finished in {elapsed / 60:.1f} min")
        note("finished", name, minutes=round(elapsed / 60, 1))

    speed = None
    run_json = out_dir / name / "run.json"
    if run_json.is_file():
        try:
            speed = json.loads(run_json.read_text()).get("ns_per_day")
        except json.JSONDecodeError:
            pass
    return {"name": name, "ok": returncode == 0,
            "minutes": round(elapsed / 60, 1), "ns_per_day": speed}


# Section 0 asks for a manifest per stage in `data/manifests/<stage>.jsonl`,
# and section 2.3 names this one. There was none: the only record of what the
# queue had done was `production.log`, which is prose for a reader, and the
# per-row `run.json`, which exists only once a row has finished. A row that is
# running, or that died, left nothing a later pass could read.
MANIFEST = ROOT / "data" / "manifests" / "md.jsonl"


def note(event: str, name: str, **fields) -> None:
    """One line in the MD manifest. Append-only, and never fatal.

    Append-only because it is a log of what happened and not a view of what is
    true now: a row that started and then failed is two lines, and collapsing
    them would throw away the thing worth knowing. `tools/ingest.py` reads the
    run directories for current state; this is the history beside it.

    Wrapped, because a queue that cannot write its own notes should still run
    the simulations. Losing the record is bad; losing 24 hours of GPU time
    because a directory was not writable would be worse.
    """
    import json

    try:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        with MANIFEST.open("a") as handle:
            handle.write(json.dumps({
                "at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "event": event, "name": name, **fields}) + "\n")
    except OSError as exc:                                     # noqa: BLE001
        log(f"  could not write {MANIFEST.name}: {exc}")


def record_failure(run_dir: Path, name: str, returncode: int,
                   elapsed: float, tail: list[str]) -> None:
    """Write the `run.json` the runner did not get to write.

    Section 0's fallback ladder ends with degrading gracefully and recording
    `status = "failed:<reason>"`. Nothing did. `run_md.py` writes `run.json`
    only when it finishes, and `tools/ingest.py` reads an absent `run.json` as
    a row still in flight, so a row that died at three in the morning showed as
    `running` on the explorer for as long as anyone cared to look, next to a
    partial trajectory that made it look plausible.

    Written here rather than in the runner because the queue is what survives
    the failure. A runner killed by the OOM reaper or by a segfault in the
    platform layer gets no chance to write anything, and those are exactly the
    failures worth recording.

    The reason comes from the last line of the row log that looks like an
    error, falling back to the exit code, and is kept short because it is shown
    in a table cell.
    """
    reason = f"exit {returncode}"
    for line in reversed(tail):
        stripped = line.strip()
        if any(mark in stripped for mark in
               ("Error", "error", "Exception", "NaN", "nan", "Killed",
                "Traceback", "assert")):
            # The log line is timestamped; the message is what is wanted.
            reason = stripped.split("] ", 1)[-1][:80]
            break

    existing = run_dir / "run.json"
    if existing.is_file():
        # The runner did finish and then something else failed. Its own record
        # is the better one, so it is left alone rather than overwritten with
        # a worse summary of the same run.
        return
    run_dir.mkdir(parents=True, exist_ok=True)
    existing.write_text(json.dumps({
        "name": name, "status": f"failed:{reason}",
        "wall_hours": round(elapsed / 3600, 2), "returncode": returncode,
    }, indent=2) + "\n")
    log(f"  recorded {name} as failed:{reason}")


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
