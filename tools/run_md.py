#!/usr/bin/env python3
"""Solvate, equilibrate and run one reference-state system.

Protocol from spec section 5.2. Common to both system kinds: PME at a 1.0 nm
cutoff, HBonds constrained, hydrogen mass repartitioned to 1.5 amu for a 4 fs
step, Langevin middle integrator at 1 ps^-1.

Water and membrane systems then diverge, and section 5.2 gives them separate
columns rather than one protocol with a flag, so they are kept apart here too:

    water     300 K, MonteCarloBarostat, 1 ns NPT equilibration
    membrane  303.15 K because POPC is fluid there and not at 300 K,
              MonteCarloMembraneBarostat at zero surface tension with XY
              isotropic and Z free, and 10 ns of NPT during which the protein
              restraints come off in six steps

A membrane system is not solvated here: `tools/build_membrane.py` has already
built it, and `build()` reloads that rather than rebuilding, which is also what
makes the bilayer a CPU job that can run while the GPU is busy elsewhere.

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

# Section 5.2's membrane column.
MEMBRANE_TEMPERATURE_K = 303.15
MEMBRANE_EQUIL_NS = 10.0
# Six steps, 10 down to zero, as section 5.2 lists them. The last is zero rather
# than a small number so production starts from a genuinely unrestrained system.
RESTRAINT_LADDER_KCAL = (10.0, 5.0, 2.0, 1.0, 0.5, 0.0)
SURFACE_TENSION_BAR_NM = 0.0


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def build(start_pdb: Path, out_dir: Path, ligand_sdf: Path | None,
          system_dir: Path | None = None):
    """Solvate and parameterise, or reload what a previous run already built.

    `system_dir` is where a prebuilt system is looked for, and defaults to the
    output directory. They are separate because a membrane system is built once
    by `tools/build_membrane.py` and then run more than once: the Tier 1 smoke
    run writes its trajectory somewhere else but has to use the same system.

    Keeping them joined caused exactly the failure that motivated the split. The
    smoke run was given its own output directory, found no `system.xml` there,
    and silently fell through to `addSolvent`, producing a 191,648-atom water box
    with no lipids in it that looked like the right size and was the wrong
    system.
    """
    from openmm import app, unit
    import openmm as mm

    source = Path(system_dir) if system_dir is not None else out_dir
    system_xml = source / "system.xml"
    solvated_pdb = source / "solvated.pdb"
    if system_xml.is_file() and solvated_pdb.is_file():
        log(f"reusing the built system in {source}")
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
    (out_dir / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    with (out_dir / "solvated.pdb").open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    return modeller.topology, modeller.positions, system


def positional_restraint(system, topology, positions, selection):
    """A flat-bottom-free harmonic restraint on chosen atoms, strength settable.

    Added once with a global `k` that is then lowered through section 5.2's
    ladder, rather than removed and re-added per step: changing a global
    parameter needs only `setParameter` on the context, while adding a force
    needs `reinitialize`, which at six steps would mean six full context
    rebuilds of a 190,000-atom system.

    Returns the force's index and the parameter name, so the caller can lower it.
    """
    import openmm as mm
    from openmm import unit

    force = mm.CustomExternalForce(
        "k_restraint*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k_restraint",
                             0.0 * unit.kilocalorie_per_mole / unit.angstrom**2)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    for index in selection:
        force.addParticle(int(index), positions[index].value_in_unit(unit.nanometer))
    return system.addForce(force)


def lipid_z_restraint(system, topology, positions, strength_kcal):
    """Hold lipid phosphorus at its starting z during minimisation (section 5.2).

    In z only. Restraining the headgroups in the plane would stop the bilayer
    finding its own area per lipid, which is the one bilayer property the run has
    to get right; restraining them in z stops a lipid being pulled through the
    midplane while the protein settles.
    """
    import openmm as mm
    from openmm import unit

    force = mm.CustomExternalForce("k_lipid*(z - z0)^2")
    force.addGlobalParameter(
        "k_lipid", strength_kcal * unit.kilocalorie_per_mole / unit.angstrom**2)
    force.addPerParticleParameter("z0")
    count = 0
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol != "P":
            continue
        if not (atom.residue.name.startswith("PO") or atom.residue.name in
                ("POPC", "PC", "OL", "PA")):
            continue
        force.addParticle(atom.index,
                          [positions[atom.index][2].value_in_unit(unit.nanometer)])
        count += 1
    if count == 0:
        return None, 0
    return system.addForce(force), count


def protein_heavy_atoms(topology):
    """Protein heavy atoms, which are what the restraints hold."""
    from openmm.app import element as elem

    standard = set()
    for atom in topology.atoms():
        residue = atom.residue
        if residue.name in ("HOH", "WAT", "NA", "CL", "POPC", "PC", "OL", "PA"):
            continue
        if residue.name.startswith("PO"):
            continue
        if atom.element is None or atom.element == elem.hydrogen:
            continue
        standard.add(atom.index)
    return sorted(standard)


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
    parser.add_argument("--system-dir", type=Path, default=None,
                        help="reuse the system.xml and solvated.pdb in this "
                             "directory instead of building one. A membrane "
                             "system is built once and run more than once.")
    parser.add_argument("--equil-ns", type=float, default=None,
                        help="override the equilibration length. Tier 1's smoke "
                             "run uses this: it exists to prove the system steps "
                             "without blowing up, and 10 ns of equilibration "
                             "would make proving that cost more than the run it "
                             "is checking.")
    parser.add_argument("--membrane", action="store_true",
                        help="use section 5.2's membrane column: 303.15 K, "
                             "membrane barostat, 10 ns staged equilibration")
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "md")
    args = parser.parse_args()

    from openmm import app, unit
    import openmm as mm

    out_dir = Path(args.out) / args.name
    out_dir.mkdir(parents=True, exist_ok=True)
    topology, positions, system = build(args.start, out_dir, args.ligand,
                                        args.system_dir)
    if args.membrane and not any(
            r.name.startswith("PO") or r.name in ("POPC", "PC", "OL", "PA")
            for r in topology.residues()):
        # The membrane protocol on a system with no lipids in it would run, and
        # every number it produced would be wrong in a way nothing downstream
        # could detect. Refused here rather than discovered in the analysis.
        raise SystemExit(
            "--membrane was given but this system contains no lipids; "
            "build it with tools/build_membrane.py first")

    temperature = MEMBRANE_TEMPERATURE_K if args.membrane else TEMPERATURE_K
    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin, FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds)

    # Restraints and the barostat go on before the Simulation is built, so the
    # context is created once. `k_restraint` starts at zero and is raised for
    # equilibration, which keeps production free of any restraint without
    # needing a second context.
    restrained = []
    if args.membrane:
        restrained = protein_heavy_atoms(topology)
        positional_restraint(system, topology, positions, restrained)
        lipid_force, n_lipid_p = lipid_z_restraint(
            system, topology, positions, RESTRAINT_LADDER_KCAL[0])
        system.addForce(mm.MonteCarloMembraneBarostat(
            PRESSURE_BAR * unit.bar,
            SURFACE_TENSION_BAR_NM * unit.bar * unit.nanometer,
            temperature * unit.kelvin,
            mm.MonteCarloMembraneBarostat.XYIsotropic,
            mm.MonteCarloMembraneBarostat.ZFree))
        log(f"membrane protocol: {temperature} K, {len(restrained):,} restrained "
            f"protein heavy atoms, {n_lipid_p} lipid P held in z")
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
    elif args.membrane:
        simulation.context.setPositions(positions)
        simulation.context.setParameter(
            "k_restraint",
            RESTRAINT_LADDER_KCAL[0] * unit.kilocalorie_per_mole / unit.angstrom**2)
        log(f"minimising with protein restrained at {RESTRAINT_LADDER_KCAL[0]} "
            "kcal/mol/A^2 and lipid P held in z")
        simulation.minimizeEnergy()
        log(f"heating to {temperature} K over {HEAT_NS} ns (NVT)")
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        simulation.step(int(HEAT_NS * steps_per_ns))

        # Section 5.2: 10 ns NPT, restraints down the ladder, and the lipid z
        # restraint off after the second step so the bilayer can relax its area
        # while the protein is still held.
        total_equil_ns = args.equil_ns if args.equil_ns is not None else MEMBRANE_EQUIL_NS
        per_step_ns = total_equil_ns / len(RESTRAINT_LADDER_KCAL)
        for index, strength in enumerate(RESTRAINT_LADDER_KCAL):
            simulation.context.setParameter(
                "k_restraint",
                strength * unit.kilocalorie_per_mole / unit.angstrom**2)
            if index == 2 and n_lipid_p:
                simulation.context.setParameter(
                    "k_lipid", 0.0 * unit.kilocalorie_per_mole / unit.angstrom**2)
                log("  lipid z restraint released")
            log(f"  NPT {per_step_ns:.2f} ns at {strength} kcal/mol/A^2")
            simulation.step(int(per_step_ns * steps_per_ns))
            simulation.saveCheckpoint(str(checkpoint))

        simulation.context.setTime(0.0 * unit.picoseconds)
        simulation.currentStep = 0
        simulation.saveCheckpoint(str(checkpoint))
        log("equilibration done, production starts unrestrained")
    else:
        simulation.context.setPositions(positions)
        log("minimising")
        simulation.minimizeEnergy()
        log(f"heating to {temperature} K over {HEAT_NS} ns (NVT)")
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        simulation.step(int(HEAT_NS * steps_per_ns))
        log(f"equilibrating {EQUIL_NS} ns (NPT)")
        system.addForce(mm.MonteCarloBarostat(PRESSURE_BAR * unit.bar,
                                              temperature * unit.kelvin))
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
        "temperature_k": temperature, "membrane": bool(args.membrane),
        "equil_ns": args.equil_ns if args.equil_ns is not None
                    else (MEMBRANE_EQUIL_NS if args.membrane else EQUIL_NS),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
