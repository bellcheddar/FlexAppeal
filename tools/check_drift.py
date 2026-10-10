#!/usr/bin/env python3
"""Energy drift and temperature stability, per run. Spec section 5.2.

Section 5.2 sets a 4 fs step with hydrogen mass repartitioning and names a
fallback: "2 fs without HMR if energy drift exceeds 0.01 kT/ns/dof". This
measures that quantity so the fallback is a decision rather than a worry.

**What the number means here, stated because it is easy to over-read.** These
runs use a Langevin middle integrator and, in production, a barostat. Both
exchange energy with a bath by design, so the total energy is not a conserved
quantity and a non-zero slope is not by itself a fault. What the threshold
usefully catches in a thermostatted run is a *systematic* trend: an integrator
losing ground at 4 fs shows up as a steady climb that the thermostat cannot
absorb, alongside a temperature that sits off its setpoint.

So three things are reported together, and any one of them alone would mislead:

    drift        the least-squares slope of total energy, in kT/ns/dof
    temperature  mean and slope against the setpoint
    noise        the residual spread, so a slope can be read against it

A slope smaller than its own uncertainty is reported as "within noise" rather
than as a number, because quoting 0.0003 kT/ns/dof from data that cannot
resolve it is a false precision.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Section 5.2's limit.
DRIFT_LIMIT_KT_PER_NS_PER_DOF = 0.01
# Boltzmann in the units StateDataReporter writes, kJ/mol/K.
BOLTZMANN = 0.00831446261815324
# How many standard errors a slope must clear before it is called a measured
# drift. One sigma is the obvious choice and is wrong for a check like this:
# with pure noise a one-sigma test calls about a third of draws resolved, so a
# QC gate built on it would flag a third of perfectly good runs. Three sigma
# brings that to roughly 0.3 per cent, which is the right trade for something
# whose job is to be believed when it fires.
RESOLVED_SIGMA = 3.0

# How far the mean temperature may sit from the setpoint before it is worth
# saying so. Three kelvin on 300 is the ordinary fluctuation of a system this
# size, and a thermostat holding 10 K off its target is not holding it.
TEMPERATURE_TOLERANCE_K = 3.0


def degrees_of_freedom(system_xml: Path) -> int | None:
    """3N minus constraints minus the three the COM remover takes.

    Read from the serialised system rather than from the atom count, because
    HBonds and rigid water remove a large fraction of them: for the membrane
    system it is 167,040 constraints against 194,032 particles, so using 3N
    would understate the drift per degree of freedom by a factor of three.
    """
    if not system_xml.is_file():
        return None
    import openmm as mm

    system = mm.XmlSerializer.deserialize(system_xml.read_text())
    massless = sum(1 for i in range(system.getNumParticles())
                   if system.getParticleMass(i)._value == 0.0)
    moving = system.getNumParticles() - massless
    removes_com = any(isinstance(f, mm.CMMotionRemover)
                      for f in system.getForces())
    return 3 * moving - system.getNumConstraints() - (3 if removes_com else 0)


def _series(path: Path) -> dict[str, list[float]]:
    """The reporter's columns, by the name it wrote rather than by position."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}

    def column(*wanted: str) -> list[float]:
        for key in rows[0]:
            clean = key.strip().strip('"#')
            if any(clean.startswith(name) for name in wanted):
                out = []
                for row in rows:
                    try:
                        out.append(float(row[key]))
                    except (TypeError, ValueError):
                        pass
                return out
        return []

    return {
        "time_ps": column("Time (ps)"),
        "total": column("Total Energy"),
        "temperature": column("Temperature"),
    }


def _fit(x: list[float], y: list[float]) -> tuple[float, float]:
    """Least-squares slope, and the standard error of that slope."""
    import numpy as np

    xs, ys = np.asarray(x, float), np.asarray(y, float)
    if xs.size < 3:
        return 0.0, float("inf")
    slope, intercept = np.polyfit(xs, ys, 1)
    residual = ys - (slope * xs + intercept)
    spread = float(np.sqrt((residual ** 2).sum() / (xs.size - 2)))
    denominator = float(((xs - xs.mean()) ** 2).sum())
    error = spread / np.sqrt(denominator) if denominator > 0 else float("inf")
    return float(slope), float(error)


def check(run_dir: Path, setpoint_k: float | None = None) -> dict:
    """Drift, temperature and noise for one run."""
    import numpy as np

    state = run_dir / "state_data.csv"
    if not state.is_file():
        return {"run": run_dir.name, "status": "no state_data.csv"}
    series = _series(state)
    if not series.get("total") or len(series["total"]) < 3:
        return {"run": run_dir.name, "status": "too few frames to fit"}

    dof = degrees_of_freedom(run_dir / "system.xml")
    if setpoint_k is None:
        meta = run_dir / "run.json"
        if meta.is_file():
            try:
                setpoint_k = json.loads(meta.read_text()).get("temperature_k")
            except json.JSONDecodeError:
                setpoint_k = None

    time_ns = [t / 1000.0 for t in series["time_ps"]]
    slope, error = _fit(time_ns, series["total"])      # kJ/mol per ns
    temperature = np.asarray(series["temperature"], float)
    mean_k = float(temperature.mean())
    kt = BOLTZMANN * (setpoint_k or mean_k)

    out: dict = {
        "run": run_dir.name, "status": "ok",
        "frames": len(series["total"]),
        "span_ns": round(time_ns[-1] - time_ns[0], 3),
        "dof": dof,
        "mean_temperature_k": round(mean_k, 2),
        "temperature_sd_k": round(float(temperature.std()), 2),
        "setpoint_k": setpoint_k,
        "drift_kj_per_mol_per_ns": round(slope, 1),
        "drift_error_kj_per_mol_per_ns": round(error, 1),
    }
    temp_slope, temp_error = _fit(time_ns, list(temperature))
    out["temperature_slope_k_per_ns"] = round(temp_slope, 3)

    if dof and kt:
        per_dof = slope / (kt * dof)
        per_dof_error = error / (kt * dof)
        out["drift_kt_per_ns_per_dof"] = round(per_dof, 6)
        out["drift_error_kt_per_ns_per_dof"] = round(per_dof_error, 6)
        out["limit"] = DRIFT_LIMIT_KT_PER_NS_PER_DOF
        # A slope inside RESOLVED_SIGMA of its own uncertainty is not a
        # measured drift.
        out["resolved"] = abs(per_dof) > RESOLVED_SIGMA * per_dof_error
        out["resolved_sigma"] = RESOLVED_SIGMA
        out["within_limit"] = abs(per_dof) <= DRIFT_LIMIT_KT_PER_NS_PER_DOF

    if setpoint_k:
        out["temperature_on_setpoint"] = abs(mean_k - setpoint_k) <= TEMPERATURE_TOLERANCE_K
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", default=None,
                        help="a run directory name, repeatable; default is all")
    parser.add_argument("--md", type=Path, default=ROOT / "data" / "md")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    directories = ([Path(args.md) / name for name in args.run] if args.run
                   else sorted(p for p in Path(args.md).iterdir() if p.is_dir()))
    results = [check(d) for d in directories]

    print(f"\nSection 5.2's limit is {DRIFT_LIMIT_KT_PER_NS_PER_DOF} kT/ns/dof.")
    print("A Langevin thermostat exchanges energy with a bath by design, so a "
          "non-zero\nslope is not a fault on its own; what matters is a "
          "systematic trend alongside a\ntemperature that has left its "
          "setpoint.\n")
    print(f"{'run':20s} {'frames':>6s} {'span':>7s} {'mean T':>8s} "
          f"{'T slope':>9s} {'drift':>12s}  verdict")
    for row in results:
        if row.get("status") != "ok":
            print(f"{row['run']:20s} {row.get('status')}")
            continue
        drift = (f"{row['drift_kt_per_ns_per_dof']:+.5f}"
                 if row.get("resolved") else "within noise")
        verdict = ("ok" if row.get("within_limit", True) else "OVER LIMIT")
        if row.get("temperature_on_setpoint") is False:
            verdict += ", T off setpoint"
        print(f"{row['run']:20s} {row['frames']:6d} {row['span_ns']:6.1f}n "
              f"{row['mean_temperature_k']:7.2f}K "
              f"{row['temperature_slope_k_per_ns']:+8.3f} {drift:>12s}  {verdict}")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwrote {args.out}")
    over = [r for r in results if r.get("within_limit") is False]
    return 1 if over else 0


if __name__ == "__main__":
    sys.exit(main())
