#!/usr/bin/env python3
"""Solvate, equilibrate and run one reference-state system.

Protocol from spec section 5.2, water systems only for now: PME at a 1.0 nm
cutoff, HBonds constrained, hydrogen mass repartitioned to 1.5 amu for a 4 fs
step, Langevin middle integrator at 1 ps^-1.

Resumable by construction (spec 0.1 rule 4). A checkpoint lands every nanosecond
and a restart picks up from it, because this is expected to be killed: the build
runs unattended, the machine sleeps, and one replica per row means a run that
cannot resume costs its whole row rather than a third of it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Section 5.2. Lengths are in nanoseconds; D009 set production to 50.
HEAT_NS, EQUIL_NS, PRODUCTION_NS = 0.5, 1.0, 50.0
TIMESTEP_FS, FRAME_PS, CHECKPOINT_NS = 4.0, 100.0, 1.0
TEMPERATURE_K, FRICTION_PER_PS, PRESSURE_BAR = 300.0, 1.0, 1.0
PADDING_NM, IONIC_STRENGTH_M = 1.2, 0.15
RESTRAINT_KCAL = 10.0


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def build(start_pdb: Path, out_dir: Path, ligand_sdf: Path | None):
    """Solvate and parameterise, or reload what a previous run already built."""
    from openmm import app, unit
    import openmm as mm

    system_xml = out_dir / "system.xml"
    solvated_pdb = out_dir / "solvated.pdb"
    if system_xml.is_file() and solvated_pdb.is_file():
        log("reusing the solvated system from a previous run")
        pdb = app.PDBFile(str(solvated_pdb))
        system = mm.XmlSerializer.deserialize(system_xml.read_text())
        return pdb.topology, pdb.positions, system

    from pdbfixer import PDBFixer

    log("PDBFixer: missing atoms, termini, protonation at pH 7.4")
    fixer = PDBFixer(filename=str(start_pdb))
    fixer.findMissingResidues()
    # Only internal gaps are built. Missing *terminal* residues are left off:
    # modelling a tail that was never resolved adds a flexible decoration that
    # then dominates the RMSF and is not evidence of anything.
    chains = list(fixer.topology.chains())
    keys = list(fixer.missingResidues)
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    forcefield_files = ["amber14-all.xml", "amber14/tip3p.xml"]
    generator = None
    if ligand_sdf is not None:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator

        log(f"parameterising {ligand_sdf.name} with openff-2.2.1 (AM1-BCC, a minute or two)")
        generator = SMIRNOFFTemplateGenerator(
            molecules=[Molecule.from_file(str(ligand_sdf))], forcefield="openff-2.2.1")

    forcefield = app.ForceField(*forcefield_files)
    if generator is not None:
        forcefield.registerTemplateGenerator(generator.generator)

    modeller = app.Modeller(fixer.topology, fixer.positions)
    log(f"solvating: {PADDING_NM} nm padding, {IONIC_STRENGTH_M} M NaCl")
    modeller.addSolvent(forcefield, model="tip3p",
                        padding=PADDING_NM * unit.nanometer,
                        ionicStrength=IONIC_STRENGTH_M * unit.molar,
                        positiveIon="Na+", negativeIon="Cl-")
    log(f"solvated: {modeller.topology.getNumAtoms():,} atoms")

    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
        rigidWater=True, hydrogenMass=1.5 * unit.amu)

    out_dir.mkdir(parents=True, exist_ok=True)
    system_xml.write_text(mm.XmlSerializer.serialize(system))
    with solvated_pdb.open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    return modeller.topology, modeller.positions, system


def platform(system=None):
    """The fastest platform that will actually build a context here.

    Section 2.2 expects OpenCL at mixed precision on Apple Silicon and that is
    not available: Apple's OpenCL implementation is single precision only, and
    asking for mixed fails context creation with "No compatible OpenCL platform
    is available", which names the platform rather than the property that caused
    it. FlexAppeal hit the same thing before and the measured cost of giving up
    and falling back to the CPU was 29.5 ns/day against 206.

    So each candidate is tried for real rather than assumed, and the precision is
    dropped to what the device supports. Section 5.2's "never single for
    production" cannot be honoured on this hardware with a GPU, so the deviation
    is logged (D010) and the energy drift is measured rather than hoped for.
    """
    import openmm as mm

    attempts = [("CUDA", {"Precision": "mixed"}), ("OpenCL", {"Precision": "mixed"}),
                ("OpenCL", {"Precision": "single"}), ("CPU", {})]
    integrator = mm.VerletIntegrator(0.001)
    for name, props in attempts:
        try:
            found = mm.Platform.getPlatformByName(name)
            if system is not None:
                mm.Context(system, mm.VerletIntegrator(0.001), found, props)
            return found, props
        except Exception as exc:                               # noqa: BLE001
            log(f"  {name} {props or '{}'}: unavailable ({str(exc)[:60]})")
    raise SystemExit("no OpenMM platform available")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True, help="run id, e.g. abl1_apo")
    parser.add_argument("--start", type=Path, required=True)
    parser.add_argument("--ligand", type=Path, default=None)
    parser.add_argument("--production-ns", type=float, default=PRODUCTION_NS)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "md")
    args = parser.parse_args()

    from openmm import app, unit
    import openmm as mm

    out_dir = Path(args.out) / args.name
    out_dir.mkdir(parents=True, exist_ok=True)
    topology, positions, system = build(args.start, out_dir, args.ligand)

    integrator = mm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin, FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds)
    plat, props = platform(system)
    simulation = app.Simulation(topology, system, integrator, plat, props)
    log(f"platform: {plat.getName()} {props}")

    checkpoint = out_dir / "state.chk"
    steps_per_ns = int(1000 / (TIMESTEP_FS / 1000))
    done_ns = 0.0
    if checkpoint.is_file():
        simulation.loadCheckpoint(str(checkpoint))
        done_ns = simulation.context.getState().getTime().value_in_unit(unit.nanoseconds)
        log(f"resumed from checkpoint at {done_ns:.1f} ns")
    else:
        simulation.context.setPositions(positions)
        log("minimising")
        simulation.minimizeEnergy()
        log(f"heating to {TEMPERATURE_K} K over {HEAT_NS} ns (NVT)")
        simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
        simulation.step(int(HEAT_NS * steps_per_ns))
        log(f"equilibrating {EQUIL_NS} ns (NPT)")
        system.addForce(mm.MonteCarloBarostat(PRESSURE_BAR * unit.bar,
                                              TEMPERATURE_K * unit.kelvin))
        simulation.context.reinitialize(preserveState=True)
        simulation.step(int(EQUIL_NS * steps_per_ns))
        simulation.saveCheckpoint(str(checkpoint))
        simulation.context.setTime(0.0 * unit.picoseconds)
        simulation.currentStep = 0
        log("equilibration done, production starts")

    target_ns = args.production_ns
    if done_ns >= target_ns:
        log(f"already at {done_ns:.1f} ns of {target_ns}; nothing to do")
        return 0

    simulation.reporters.append(app.DCDReporter(
        str(out_dir / "trajectory.dcd"), int(FRAME_PS / (TIMESTEP_FS / 1000)),
        append=(out_dir / "trajectory.dcd").is_file()))
    simulation.reporters.append(app.StateDataReporter(
        str(out_dir / "state_data.csv"), int(FRAME_PS / (TIMESTEP_FS / 1000)),
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True,
        speed=True, append=(out_dir / "state_data.csv").is_file()))
    simulation.reporters.append(app.CheckpointReporter(
        str(checkpoint), int(CHECKPOINT_NS * steps_per_ns)))

    remaining = int((target_ns - done_ns) * steps_per_ns)
    log(f"producing {target_ns - done_ns:.1f} ns ({remaining:,} steps)")
    started = time.time()
    simulation.step(remaining)
    elapsed = time.time() - started
    ns_per_day = (target_ns - done_ns) / (elapsed / 86400)
    log(f"done in {elapsed / 3600:.2f} h ({ns_per_day:.1f} ns/day)")

    (out_dir / "run.json").write_text(json.dumps({
        "name": args.name, "production_ns": target_ns,
        "platform": plat.getName(), "ns_per_day": round(ns_per_day, 1),
        "wall_hours": round(elapsed / 3600, 2),
        "atoms": topology.getNumAtoms(), "status": "complete",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
