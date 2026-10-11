#!/usr/bin/env python3
"""Section 2.1's hardware probe, and the settings section 2.2 derives from it.

Writes three files under `config/`:

  hardware.toml    what this machine is, measured rather than assumed
  tuning.toml      section 2.2's table, derived from the probe
  benchmarks.toml  measured ns/day per system and platform

Section 2.2 ends "Every compute budget in Section 5.6 is derived from these
measurements, never assumed", and until now none of the three files existed.
The settings were right, because `tools/run_md.py` works them out at run time
and logs what it chose, but nothing wrote down what the machine is or what it
managed, so a budget could not be checked against anything.

Safe to re-run. Every value is either measured here or read from a run that has
already finished, and a key that cannot be produced is omitted rather than
written as a zero, which is section 0's rule and matters most in a file whose
whole job is to be trusted later.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Section 2.2's fixed rules. The ones that are not derived from a measurement
# are here so the file says where every number came from.
IO_WORKERS = 8
CONCURRENT_MD = 1
MEMORY_HEADROOM_GB = 24
CPU_WORKERS_RESERVED = 2
CPU_WORKERS_MINIMUM = 2


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def shown(path: Path) -> str:
    """A path for the log, relative to the project where that makes sense.

    `--out` is allowed to point anywhere, and `relative_to` raises rather than
    giving up when it does, so the short form is an attempt and not an
    assumption. `tools/make_figures.py` has the same helper for the same
    reason, and this file was written with the same bug in it.
    """
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def sysctl(name: str) -> str | None:
    """One sysctl value, or nothing if this machine does not have it."""
    try:
        found = subprocess.run(["sysctl", "-n", name], capture_output=True,
                               text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    value = found.stdout.strip()
    return value or None


def gpu_description() -> dict:
    """The GPU, from `system_profiler`, which is slow so it is called once."""
    try:
        found = subprocess.run(
            ["system_profiler", "-json", "SPDisplaysDataType"],
            capture_output=True, text=True, timeout=60)
        payload = json.loads(found.stdout)
    except (OSError, subprocess.SubprocessError, json.JSONDecodeError):
        return {}
    cards = payload.get("SPDisplaysDataType") or []
    if not cards:
        return {}
    card = cards[0]
    out = {}
    for key, name in (("sppci_model", "model"),
                      ("sppci_cores", "cores"),
                      ("spdisplays_vendor", "vendor")):
        if card.get(key):
            out[name] = str(card[key])
    return out


def disk_free_gb(path: Path) -> float | None:
    """Free space where the simulations are written, not where / is."""
    try:
        found = subprocess.run(["df", "-k", str(path)], capture_output=True,
                               text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    lines = found.stdout.strip().splitlines()
    if len(lines) < 2:
        return None
    parts = lines[-1].split()
    # df -k: Filesystem 1024-blocks Used Available ...
    try:
        return round(int(parts[3]) / 1024 / 1024, 1)
    except (IndexError, ValueError):
        return None


def openmm_platforms() -> list[dict]:
    """Every platform OpenMM has here, with the speed each one claims.

    The speed is OpenMM's own estimate and not a measurement, which is why
    `tools/run_md.py` tries each candidate for real instead of trusting the
    order. Recorded anyway, because what is installed is a fact about the
    machine and the next person should not have to find out by running it.
    """
    try:
        import openmm as mm
    except ImportError:
        return []
    found = []
    for index in range(mm.Platform.getNumPlatforms()):
        platform = mm.Platform.getPlatform(index)
        found.append({"name": platform.getName(),
                      "claimed_speed": round(float(platform.getSpeed()), 3)})
    return found


def probe() -> dict:
    """Section 2.1's facts about this machine."""
    memory = sysctl("hw.memsize")
    hardware: dict = {"probed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    for key, name in (("machdep.cpu.brand_string", "cpu"),
                      ("hw.perflevel0.physicalcpu", "performance_cores"),
                      ("hw.perflevel1.physicalcpu", "efficiency_cores")):
        value = sysctl(key)
        if value is None:
            continue
        hardware[name] = int(value) if value.isdigit() else value
    if memory and memory.isdigit():
        hardware["memory_gb"] = round(int(memory) / 1024 ** 3, 1)
    free = disk_free_gb(ROOT)
    if free is not None:
        hardware["disk_free_gb"] = free
    gpu = gpu_description()
    if gpu:
        hardware["gpu"] = gpu
    platforms = openmm_platforms()
    if platforms:
        hardware["openmm_platforms"] = platforms
    return hardware


def chosen_platform() -> dict:
    """What the runs actually used, from the rows that have finished.

    Read rather than re-benchmarked. `run_md.py` already tries each candidate
    against the real system and records what it kept, and benchmarking a 20 ps
    run again here would compete with the queue for the one GPU.
    """
    out: dict = {}
    for run_json in sorted((ROOT / "data" / "md").glob("*/run.json")):
        try:
            done = json.loads(run_json.read_text())
        except json.JSONDecodeError:
            continue
        if done.get("platform"):
            out.setdefault(done["platform"], 0)
            out[done["platform"]] += 1
    return out


def tuning(hardware: dict, platforms: dict) -> dict:
    """Section 2.2's table, with each value's origin recorded beside it."""
    settings: dict = {}

    if platforms:
        # The one the finished rows used, not the one OpenMM ranks first.
        settings["md_platform"] = max(platforms, key=platforms.get)
    # Section 5.2 says never single for production. D010 records that Apple's
    # OpenCL is single precision only, so `run_md.py` drops to it and measures
    # the drift instead. Recording what is actually used, not what is wanted.
    settings["md_precision"] = ("single" if settings.get("md_platform") == "OpenCL"
                                else "mixed")

    cores = hardware.get("performance_cores")
    if isinstance(cores, int):
        settings["cpu_workers"] = max(CPU_WORKERS_MINIMUM,
                                      cores - CPU_WORKERS_RESERVED)
    settings["io_workers"] = IO_WORKERS
    settings["concurrent_md"] = CONCURRENT_MD
    settings["memory_headroom_gb"] = MEMORY_HEADROOM_GB
    return settings


def benchmarks() -> dict:
    """Measured ns/day per system and platform, from the runs themselves.

    `run.json` is written when a row finishes, so this fills in as the queue
    does. A row still in flight contributes its current rate from the state
    table instead, marked as partial, because a 50 ns row takes most of a day
    and a benchmarks file that says nothing until then is no use to anyone
    planning the next one.
    """
    out: dict = {}
    for directory in sorted((ROOT / "data" / "md").glob("*/")):
        name = directory.name
        run_json = directory / "run.json"
        if run_json.is_file():
            try:
                done = json.loads(run_json.read_text())
            except json.JSONDecodeError:
                done = {}
            row = {k: done[k] for k in
                   ("platform", "ns_per_day", "atoms", "production_ns",
                    "wall_hours")
                   if done.get(k) is not None}
            if row:
                row["complete"] = done.get("status") == "complete"
                out[name] = row
            continue

        state = directory / "state_data.csv"
        if not state.is_file():
            continue
        try:
            last = state.read_text().strip().splitlines()[-1].split(",")
            rate, done_ns = float(last[-1]), float(last[1]) / 1000.0
        except (IndexError, ValueError):
            continue
        if rate > 0:
            out[name] = {"ns_per_day": round(rate, 1),
                         "produced_ns": round(done_ns, 2),
                         "complete": False}
    return out


def as_toml(table: dict, *, header: str) -> str:
    """A small TOML writer, because `tomli_w` is not in this environment.

    Only what these three files need: scalars, arrays of tables, and one level
    of nesting. Anything else raises rather than writing something that parses
    differently from what was meant.
    """
    def scalar(value) -> str:
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (int, float)):
            return repr(value)
        return json.dumps(str(value))

    lines = [f"# {line}" for line in header.strip().splitlines()]
    lines.append("")
    nested, arrays = {}, {}
    for key, value in table.items():
        if isinstance(value, dict):
            nested[key] = value
        elif isinstance(value, list) and value and isinstance(value[0], dict):
            arrays[key] = value
        elif isinstance(value, list):
            lines.append(f"{key} = [{', '.join(scalar(v) for v in value)}]")
        else:
            lines.append(f"{key} = {scalar(value)}")
    for key, value in nested.items():
        lines.append("")
        lines.append(f"[{key}]")
        for inner, found in value.items():
            if isinstance(found, (dict, list)):
                raise ValueError(f"{key}.{inner} nests deeper than this writer goes")
            lines.append(f"{inner} = {scalar(found)}")
    for key, rows in arrays.items():
        for row in rows:
            lines.append("")
            lines.append(f"[[{key}]]")
            for inner, found in row.items():
                lines.append(f"{inner} = {scalar(found)}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "config")
    parser.add_argument("--print", action="store_true",
                        help="show what would be written and write nothing")
    args = parser.parse_args()

    hardware = probe()
    platforms = chosen_platform()
    settings = tuning(hardware, platforms)
    measured = benchmarks()

    files = {
        "hardware.toml": (hardware, "Section 2.1's probe. Regenerate with "
                                    "tools/hwprobe.py; do not edit by hand."),
        "tuning.toml": (settings, "Section 2.2's settings, derived from "
                                  "hardware.toml and from what the finished "
                                  "runs actually used."),
        "benchmarks.toml": ({"runs": [dict(name=k, **v)
                                      for k, v in measured.items()]},
                            "Measured ns/day per system and platform. Section "
                            "2.2: every compute budget is derived from these, "
                            "never assumed."),
    }

    args.out.mkdir(parents=True, exist_ok=True)
    for name, (table, header) in files.items():
        text = as_toml(table, header=header)
        if args.print:
            print(f"--- {name} ---")
            print(text)
            continue
        (args.out / name).write_text(text)
        log(f"wrote {shown(args.out / name)}")

    if not args.print:
        log(f"{hardware.get('cpu', 'unknown cpu')}, "
            f"{hardware.get('performance_cores', '?')} performance cores, "
            f"{hardware.get('memory_gb', '?')} GB, "
            f"{len(hardware.get('openmm_platforms', []))} OpenMM platforms")
        log(f"{len(measured)} runs benchmarked, "
            f"{sum(1 for r in measured.values() if r.get('complete'))} finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
