"""MM/GBSA binding free energies, computed by the reference implementation.

This fills the one measurement FlexAppeal could not make. Until now the energy
sub-score was pending for two reasons at once -- no cohort to take a percentile
in, and no number to take it of -- and `--energies` was the only way one arrived.
This removes the second reason.

**It is AmberTools' MMPBSA.py, not a reimplementation, and that is the whole
point.** The energies this produces join a cohort whose existing members came out
of `MMPBSA_igb_2` on a water-stripped GAFF2 topology, and a percentile across a
cohort measured two different ways is not a percentile of anything. So the
protocol is pinned to match: single trajectory, `igb=2`, mbondi2 radii, salt at
the simulation's own ionic strength. Where a choice existed it was made by
looking at what produced the numbers these will be ranked against.

The same argument, made twice already today, is why this module shells out to a
program instead of adding a `GBSAOBCForce` and computing it in OpenMM. That would
have been less code and a day's less care, and it would have put FlexAppeal on its
own ruler again.

What it does, per complex:

1. Rebuilds the parameterised system **unconstrained**. The simulated system in
   `system.xml` constrains every bond to hydrogen, so those bonds carry no
   force-field term and a topology written from it has bonds with no type at all,
   which is where the first attempt at this stopped. An energy calculation wants
   the full potential anyway.
2. Writes an Amber topology through ParmEd and sets the GB radii. A topology
   written from an OpenMM system has `RADII` and `SCREEN` all zero, because the
   system it came from had no implicit-solvent force to carry them. MMPBSA.py
   then divides by them and reports every `EGB` as `nan`, which it does without
   failing: the first run here produced a complete, well-formatted results file
   whose every total was `nan`. `changeRadii mbondi2` fixes it, and mbondi2 is the
   set `igb=2` expects.
3. Splits receptor and ligand topologies from the complex with `ante-MMPBSA.py`,
   so all three are guaranteed consistent with each other.
4. Runs MMPBSA.py over a subsample of frames and parses the difference block.

How to tell it worked. In a single-trajectory calculation the bonded terms cancel
**exactly**, because the receptor and the ligand are read out of the complex's own
coordinates. `BOND`, `ANGLE`, `DIHED` and both 1-4 terms must come back as
-0.0000, and `_check_cancellation` refuses the result if they do not: a bonded
term that fails to cancel means the three topologies disagree about the molecule,
which is the one error in this pipeline that produces a plausible number rather
than an obvious failure.

What this is not. There is no entropy term, which is the usual choice and the
usual caveat: single-trajectory MM/GBSA with no `-nmode` is systematically
over-negative, typically by several kcal/mol. Benzene in T4 lysozyme L99A comes
out around -10.7 kcal/mol here against roughly -5.2 from experiment, which is the
expected direction and size of that error. These numbers rank ligands against each
other; they are not predictions of an affinity.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# --- the protocol, pinned to the cohort these numbers join -------------------
# igb=2 is OBC(I). It is what `MMPBSA_igb_2` in the existing package used, and
# mbondi2 is the radius set it is parameterised against: pairing igb=2 with the
# default Bondi radii is a known way to get numbers that look fine and are not.
IGB = 2
RADII_SET = "mbondi2"
# Salt screening, mol/L. Defaults to the ionic strength the system was actually
# solvated at rather than to a convention, and the caller passes it.
DEFAULT_SALT = 0.15

# Frames to score. MMPBSA.py takes roughly a third of a second per frame for a
# few thousand atoms, so this is about the analysis staying a few minutes rather
# than about the statistics: the standard error over 200 frames of a well-behaved
# complex is already a few tenths of a kcal/mol.
MAX_ENERGY_FRAMES = 200
# Below this there is no point reporting a standard deviation.
MIN_ENERGY_FRAMES = 10

# The terms MMPBSA.py reports in its difference block, in its own spelling,
# because `package._energy_terms` carries them through to DeepFlex unaltered.
DIFFERENCE_TERMS = (
    "BOND", "ANGLE", "DIHED", "VDWAALS", "EEL", "1-4 VDW", "1-4 EEL",
    "EGB", "ESURF", "DELTA G gas", "DELTA G solv", "DELTA TOTAL",
)
# The ones that must cancel in a single-trajectory calculation, and the slack
# allowed. Anything above this is a topology disagreement, not rounding.
CANCELLING_TERMS = ("BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL")
CANCELLATION_TOLERANCE = 0.01

_ROW = re.compile(
    r"^(?P<term>BOND|ANGLE|DIHED|VDWAALS|EEL|1-4 VDW|1-4 EEL|EGB|ESURF|"
    r"DELTA G gas|DELTA G solv|DELTA TOTAL)\s+"
    r"(?P<mean>-?\d+\.\d+|nan)\s+(?P<sd>-?\d+\.\d+|nan)\s+"
    r"(?P<sem>-?\d+\.\d+|nan)\s*$")


class EnergyError(RuntimeError):
    """The calculation could not be run or its result could not be trusted."""


@dataclass
class Energies:
    """One complex's binding energy, shaped for `package.assemble(energies=...)`."""

    compound: str
    receptor: str
    dg: float
    sd: float
    sem: float
    n_frames: int
    terms: dict[str, list[float]] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def as_supplied(self) -> dict[str, Any]:
        """The `--energies` shape, so nothing in `package` needs to know these
        were computed here rather than handed in from a study."""
        entry = {
            "mmpbsa_gb_dg": round(self.dg, 4),
            "mmpbsa_gb_sd": round(self.sd, 4),
            "terms": {"gb": self.terms},
        }
        # Under both keys, with the same value. `<compound>/<receptor>` is the
        # specific one and wins where a compound was run against several
        # receptors, but the receptor token is derived from the structure
        # filename or the accession and a caller that derives it differently
        # would silently lose the energy. Two identical entries cost nothing and
        # make the lookup unable to miss.
        return {f"{self.compound}/{self.receptor}": entry, self.compound: entry}


def available() -> tuple[bool, str | None]:
    """Whether the toolchain is here, and what is missing if not.

    Checked rather than assumed: AmberTools reaches a ligand bundle because
    AM1-BCC needs it, so it is present for the runs that can use this and absent
    for a run charged by Gasteiger. A missing tool costs the energy sub-score and
    must not cost the analysis.
    """
    for tool in ("MMPBSA.py", "ante-MMPBSA.py"):
        if shutil.which(tool) is None:
            return False, (
                f"{tool} is not on PATH, so MM/GBSA cannot be computed. It comes "
                f"with AmberTools, which a bundle installs when the ligand needs "
                f"AM1-BCC charges")
    try:
        import parmed  # noqa: F401, PLC0415
    except ImportError:
        return False, "ParmEd is not installed, so no Amber topology can be written"
    return True, None


# ===========================================================================
#  Topologies
# ===========================================================================


def _unconstrained_system(topology, forcefield):
    """The complex's potential with nothing held rigid.

    `constraints=None` and `rigidWater=False` on purpose. The simulated system
    constrains bonds to hydrogen, which removes their force-field terms, and a
    topology written from that has typeless bonds that ParmEd refuses. It is also
    the right physics here: a single-point energy wants every term present, and
    the bonded ones cancel in the difference either way.
    """
    from openmm import app  # noqa: PLC0415

    return forcefield.createSystem(
        topology, nonbondedMethod=app.NoCutoff,
        constraints=None, rigidWater=False, removeCMMotion=False)


def _write_topologies(work: Path, topology, system, positions,
                      ligand_resnames: list[str]) -> list[str]:
    """Write complex, receptor and ligand Amber topologies into `work`.

    Returns any warnings. Raises EnergyError if the split fails, because every
    number downstream depends on the three agreeing.
    """
    import parmed  # noqa: PLC0415
    import parmed.tools as parmed_tools  # noqa: PLC0415

    warnings: list[str] = []
    structure = parmed.openmm.load_topology(topology, system=system, xyz=positions)
    complex_parm = work / "complex.prmtop"
    structure.save(str(complex_parm), overwrite=True)

    # The radii, without which every EGB is silently nan. See the module note.
    parm = parmed.amber.AmberParm(str(complex_parm))
    if all(value == 0 for value in parm.parm_data["RADII"]):
        parmed_tools.changeRadii(parm, RADII_SET).execute()
        parm.save(str(complex_parm), overwrite=True)
    else:                                                     # pragma: no cover
        warnings.append(
            f"the topology already carried GB radii "
            f"({parm.parm_data.get('RADIUS_SET')}), so {RADII_SET} was not applied")

    mask = ",".join(f":{name}" for name in ligand_resnames)
    result = subprocess.run(
        ["ante-MMPBSA.py", "-p", str(complex_parm),
         "-r", str(work / "receptor.prmtop"), "-l", str(work / "ligand.prmtop"),
         "-n", mask],
        cwd=work, capture_output=True, text=True, check=False)
    if result.returncode != 0 or not (work / "ligand.prmtop").is_file():
        raise EnergyError(
            f"the receptor and ligand topologies could not be split from the "
            f"complex on mask {mask!r}: "
            f"{(result.stderr or result.stdout).strip()[-400:]}")
    return warnings


# ===========================================================================
#  Parsing
# ===========================================================================


def _difference_block(text: str) -> dict[str, list[float]]:
    """The `Complex - Receptor - Ligand` block, as {term: [mean, sd, sem]}.

    Only that block. The three absolute blocks above it are each thousands of
    kcal/mol and only their difference is a binding quantity, which is also what
    `package._energies_for` refuses values the size of.
    """
    marker = "Differences (Complex - Receptor - Ligand)"
    if marker not in text:
        raise EnergyError(
            "MMPBSA.py produced no difference block, so the calculation did not "
            "complete. Its output is kept beside the results file.")
    block = text.split(marker, 1)[1]

    terms: dict[str, list[float]] = {}
    for line in block.splitlines():
        match = _ROW.match(line.strip())
        if match is None:
            continue
        values = [match.group("mean"), match.group("sd"), match.group("sem")]
        if any(value == "nan" for value in values):
            raise EnergyError(
                f"MMPBSA.py reported {match.group('term')} as nan. The usual "
                f"cause is a topology with no GB radii, which makes every "
                f"solvation term undefined while the run still reports success.")
        terms[match.group("term")] = [float(value) for value in values]
        if match.group("term") == "DELTA TOTAL":
            break

    if "DELTA TOTAL" not in terms:
        raise EnergyError("MMPBSA.py reported no total for the difference block")
    return terms


def _check_cancellation(terms: dict[str, list[float]]) -> None:
    """Refuse a result whose bonded terms did not cancel.

    In a single-trajectory calculation they must, exactly: the receptor and the
    ligand are the complex's own coordinates, so every internal term appears on
    both sides of the subtraction. A residual means the three topologies describe
    different molecules -- a mis-split mask, or a prmtop paired with the wrong
    trajectory -- and that error produces a plausible-looking number rather than
    a failure, which is why it is checked rather than trusted.
    """
    offenders = {
        term: terms[term][0] for term in CANCELLING_TERMS
        if term in terms and abs(terms[term][0]) > CANCELLATION_TOLERANCE
    }
    if offenders:
        detail = ", ".join(f"{k} {v:+.3f}" for k, v in offenders.items())
        raise EnergyError(
            f"the bonded terms did not cancel ({detail} kcal/mol), which in a "
            f"single-trajectory calculation means the complex, receptor and "
            f"ligand topologies do not describe the same molecule. Refusing the "
            f"energy rather than reporting it.")


# ===========================================================================
#  The calculation
# ===========================================================================


def mmgbsa(topology_pdb: Path, trajectory: Path, *, ligand_resnames: list[str],
           receptor: str, compound: str, forcefield, salt: float = DEFAULT_SALT,
           max_frames: int = MAX_ENERGY_FRAMES, start_frame: int = 0,
           keep_in: Path | None = None, progress=None) -> tuple[Energies | None,
                                                                list[str]]:
    """Compute one complex's MM/GBSA binding energy.

    `forcefield` is an already-built `openmm.app.ForceField` with the ligand's
    template generator registered, because building one means parameterising the
    ligand again and the caller has already paid for that.

    Returns (Energies, warnings), or (None, reasons) where it could not be run.
    Never raises for an absent toolchain or a refused result: the energy is one
    sub-score of six and losing it must not cost the other five.
    """
    import mdtraj as md  # noqa: PLC0415
    from openmm import app  # noqa: PLC0415

    ok, reason = available()
    if not ok:
        return None, [reason]

    warnings: list[str] = []
    work = Path(tempfile.mkdtemp(prefix="flexappeal-mmgbsa-"))
    try:
        if progress is not None:
            progress("building the Amber topologies")
        pdb = app.PDBFile(str(topology_pdb))
        system = _unconstrained_system(pdb.topology, forcefield)
        warnings.extend(_write_topologies(
            work, pdb.topology, system, pdb.positions, ligand_resnames))

        # Subsample, and say by how much. The whole production run is read, which
        # is the right window here: the restrained equilibration is a separate
        # phase before it, so there is no discarded head in this trajectory.
        total = len(md.open(str(trajectory)))
        stride = max(1, -(-(total - start_frame) // max_frames))
        traj = md.load(str(trajectory), top=str(topology_pdb), stride=stride)
        if start_frame:
            traj = traj[max(0, start_frame // stride):]
        if traj.n_atoms != pdb.topology.getNumAtoms():
            raise EnergyError(
                f"the trajectory has {traj.n_atoms} atoms and the topology "
                f"{pdb.topology.getNumAtoms()}; they are not the same selection")
        if traj.n_frames < MIN_ENERGY_FRAMES:
            return None, [
                f"only {traj.n_frames} frame(s) to score, below the "
                f"{MIN_ENERGY_FRAMES} needed for a usable mean and spread"]

        traj.save_netcdf(str(work / "traj.nc"))
        if stride > 1:
            warnings.append(
                f"the MM/GBSA mean is over {traj.n_frames} frames, every "
                f"{stride} of {total}: the calculation is about a third of a "
                f"second a frame and the standard error is already a few tenths "
                f"of a kcal/mol at this many")

        recipe = (
            f"Single-trajectory MM/GBSA written by FlexAppeal. igb={IGB} with "
            f"{RADII_SET} radii, matching the protocol the cohort these numbers\n"
            f"join was measured with. No entropy term, so the result is "
            f"systematically over-negative and ranks rather than predicts.\n"
            f"&general\n"
            f"  startframe=1, endframe=99999, interval=1, verbose=2, keep_files=0,\n"
            f"/\n&gb\n  igb={IGB}, saltcon={salt},\n/\n")
        (work / "mmpbsa.in").write_text(recipe)

        if progress is not None:
            progress(f"MM/GBSA over {traj.n_frames} frames")
        result = subprocess.run(
            ["MMPBSA.py", "-O", "-i", "mmpbsa.in", "-o", "FINAL_RESULTS_MMPBSA.dat",
             "-cp", "complex.prmtop", "-rp", "receptor.prmtop",
             "-lp", "ligand.prmtop", "-y", "traj.nc"],
            cwd=work, capture_output=True, text=True, check=False)
        report = work / "FINAL_RESULTS_MMPBSA.dat"
        if result.returncode != 0 or not report.is_file():
            raise EnergyError(
                f"MMPBSA.py exited {result.returncode}: "
                f"{(result.stderr or result.stdout).strip()[-400:]}")

        text = report.read_text()
        if keep_in is not None:
            keep_in.mkdir(parents=True, exist_ok=True)
            (keep_in / "FINAL_RESULTS_MMPBSA.dat").write_text(text)
            (keep_in / "mmpbsa.in").write_text(recipe)

        terms = _difference_block(text)
        _check_cancellation(terms)
        total_term = terms["DELTA TOTAL"]

        return Energies(
            compound=compound, receptor=receptor,
            dg=total_term[0], sd=total_term[1], sem=total_term[2],
            n_frames=traj.n_frames,
            terms={k: v for k, v in terms.items() if k in DIFFERENCE_TERMS},
            notes=warnings,
        ), warnings

    except (EnergyError, OSError, ValueError) as exc:
        return None, [f"MM/GBSA was not computed: {exc}"]
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ===========================================================================
#  Building the force field a run used
# ===========================================================================

# Ligand force field prefix -> the openmmforcefields generator that provides it.
#
# A second copy of `bundle._LIGAND_GENERATOR`, deliberately, because this module
# travels inside a bundle and `bundle` does not: importing it there would need
# the templates, the schema and the option registry for a three-entry mapping.
# `tests/test_energy.py` asserts the two agree, which is the same arrangement
# `console.py`'s theme has for the same reason.
LIGAND_GENERATOR = {
    "openff": "SMIRNOFFTemplateGenerator",
    "gaff": "GAFFTemplateGenerator",
    "espaloma": "EspalomaTemplateGenerator",
}


def forcefield_for_run(run_dir: Path) -> tuple[Any, list[str], list[str]]:
    """Rebuild the force field a run used, for the energy calculation.

    Returns (ForceField, ligand resnames, warnings). The protein force field is
    loaded on its own: the complex topology has no water or ions in it and the
    calculation is NoCutoff with implicit solvent supplied by MMPBSA.py, so a
    water model would be parameters for atoms that are not there.

    The ligand's chemistry comes from the bundle's own `ligands/<NAME>.sdf`, the
    same file the run was parameterised from, so the charges here are the charges
    that were simulated.
    """
    import json  # noqa: PLC0415

    from openmm import app  # noqa: PLC0415

    config_path = run_dir / "config.json"
    if not config_path.is_file():
        raise EnergyError(
            f"{run_dir} has no config.json, so the force field the run used "
            f"cannot be reconstructed")
    cfg = json.loads(config_path.read_text()).get("config", {})

    warnings: list[str] = []
    forcefield = app.ForceField(cfg.get("protein_ff") or "amber14-all.xml")

    sdfs = sorted((run_dir / "ligands").glob("*.sdf")) if (
        run_dir / "ligands").is_dir() else []
    if not sdfs:
        raise EnergyError(
            f"{run_dir}/ligands holds no .sdf, so the ligand cannot be "
            f"parameterised and no binding energy can be computed")

    ligand_ff = str(cfg.get("ligand_ff") or "openff-2.2.1")
    generator_name = next(
        (value for prefix, value in LIGAND_GENERATOR.items()
         if ligand_ff.startswith(prefix)), None)
    if generator_name is None:
        raise EnergyError(
            f"{ligand_ff} is not a ligand force field this can rebuild "
            f"({', '.join(sorted(LIGAND_GENERATOR))})")

    import openmmforcefields.generators as generators  # noqa: PLC0415
    from openff.toolkit import Molecule  # noqa: PLC0415

    molecules, resnames = [], []
    for path in sdfs:
        molecules.append(Molecule.from_file(str(path)))
        resnames.append(path.stem)

    generator = getattr(generators, generator_name)(
        molecules=molecules, forcefield=ligand_ff)
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield, resnames, warnings


def for_run_directory(run_dir: Path, output_dir: Path, *, receptor: str,
                      salt: float | None = None,
                      max_frames: int = MAX_ENERGY_FRAMES,
                      keep_in: Path | None = None,
                      progress=None) -> tuple[dict[str, Any], list[str]]:
    """Compute MM/GBSA for every ligand in one run, in `--energies` shape.

    The shape is the point: `package.assemble` already accepts energies supplied
    from a study, so computing them here needs nothing new on that side and the
    two paths cannot disagree about how a number reaches the file.
    """
    import json  # noqa: PLC0415

    ok, reason = available()
    if not ok:
        return {}, [reason]

    topology = output_dir / "topology.pdb"
    trajectory = next(iter(sorted(output_dir.glob("trajectory.*"))), None)
    if not topology.is_file() or trajectory is None:
        return {}, [
            f"{output_dir} has no topology.pdb and trajectory to score, so no "
            f"binding energy was computed"]

    try:
        forcefield, resnames, warnings = forcefield_for_run(run_dir)
    except EnergyError as exc:
        return {}, [f"MM/GBSA was not computed: {exc}"]

    if salt is None:
        cfg = json.loads((run_dir / "config.json").read_text()).get("config", {})
        salt = float(cfg.get("ionic_strength") or DEFAULT_SALT)

    # One calculation covering every kept ligand together. A job simulates one
    # compound, and `package` takes the first kept heteroatom as its identity, so
    # splitting them here would invent complexes the rest of the file does not
    # have.
    result, more = mmgbsa(
        topology, trajectory, ligand_resnames=resnames, receptor=receptor,
        compound=resnames[0], forcefield=forcefield, salt=salt,
        max_frames=max_frames, keep_in=keep_in, progress=progress)
    warnings.extend(more)
    if result is None:
        return {}, warnings
    return result.as_supplied(), warnings
