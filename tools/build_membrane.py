#!/usr/bin/env python3
"""Embed an oriented receptor in a POPC bilayer.

Spec section 5.1 step 2. OpenMM's `Modeller.addMembrane` is the primary builder
and `packmol-memgen` the documented fallback; POPC only in v1, because
`addMembrane` offers no cholesterol and a POPC:cholesterol bilayer means Lipid21
and a different route, which section 5.1 puts out of scope.

The receptor must already be in the membrane frame. `tools/build_b2ar.py` does
that by superposing onto OPM, and this refuses to build if the check does not
hold: a receptor embedded at the wrong tilt gives a bilayer that looks fine in
the atom count and is wrong in every measurement taken from it.

Separate from the runner because this is CPU work that can proceed while another
system has the GPU, which is the concurrency section 2.3 asks for.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Section 5.1. POPC is fluid at 303.15 K, which is why the membrane runs sit
# there rather than at the 300 K the water systems use.
LIPID = "POPC"
PADDING_NM = 1.5
FALLBACK_PADDING_NM = 1.2
IONIC_STRENGTH_M = 0.15
EXPECTED_ATOMS = (80_000, 110_000)
MAX_ATOMS = 150_000


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def check_orientation(pdb_path: Path) -> dict:
    """Refuse to build on a receptor that is not in the membrane frame."""
    import gemmi

    from flexappeal.refstates import Catalogue, receptor_residues
    from flexappeal.refstates.gpcr import build_numbering

    numbering = build_numbering(
        receptor_residues(Catalogue(ROOT / "data" / "cache" / "gpcrdb")))
    structure = gemmi.read_structure(str(pdb_path))
    structure.setup_entities()
    z = {}
    for chain in structure[0]:
        for residue in chain:
            atom = residue.find_atom("CA", "*")
            if atom is not None:
                z[residue.seqid.num] = atom.pos.z

    import numpy as np

    belt = [abs(z[n]) for n, letter in numbering.residue_letter.items()
            if letter in ("W", "Y") and n in z
            and numbering.segment.get(n, "").startswith("TM")]
    basics = [z[n] for n, letter in numbering.residue_letter.items()
              if letter in ("R", "K") and n in z
              and numbering.segment.get(n, "") in ("ICL1", "ICL2", "ICL3", "H8")]
    belt_median = float(np.median(belt)) if belt else 0.0
    basic_median = float(np.median(basics)) if basics else 0.0
    passes = 10.0 <= belt_median <= 17.0 and basic_median < -13.0
    return {"belt_median": belt_median, "basic_median": basic_median,
            "passes": passes}


def ligand_molecule(start_pdb: Path, resname: str, sdf_path: Path):
    """Recover the ligand as a fully specified, hydrogen-bearing molecule.

    `addMembrane` calls `createSystem` internally, so the SMIRNOFF generator has
    to match the residue that is actually in the topology. Handing it the CCD
    ideal SDF does not: the deposited ligand is heavy atoms only, and the
    generator's template then carries twelve hydrogens the topology has never
    heard of, which surfaces as "No template found for residue 322 (CAU)".

    The fix is the one `flexappeal/runtime/run.py.j2` already uses for bundles:
    take the chemistry from the CCD definition, the coordinates from the
    deposition, add the hydrogens here, and merge that molecule into the system
    so topology and template are built from the same object. PDBFixer cannot do
    it, having no template for a ligand in the first place.
    """
    from openff.toolkit import Molecule
    from openmm import app
    from rdkit import Chem
    from rdkit.Chem import AllChem

    source = app.PDBFile(str(start_pdb))
    present = [r for r in source.topology.residues() if r.name == resname]
    if len(present) != 1:
        raise SystemExit(f"expected exactly one {resname}, found {len(present)}")

    keep = {a.index for a in present[0].atoms()}
    fragment = app.Modeller(source.topology, source.positions)
    fragment.delete([a for a in fragment.topology.atoms() if a.index not in keep])
    block = io.StringIO()
    app.PDBFile.writeFile(fragment.topology, fragment.positions, block)

    rd = Chem.MolFromPDBBlock(block.getvalue(), removeHs=True)
    # RemoveAllHs, not RemoveHs: the latter keeps stereo-defining hydrogens and
    # the atom counts then differ from the PDB fragment for no visible reason.
    template = Chem.RemoveAllHs(Chem.MolFromMolFile(str(sdf_path)))
    if rd is None or template is None:
        raise SystemExit(f"could not read {resname} from the structure or the CCD definition")
    rd = AllChem.AssignBondOrdersFromTemplate(template, rd)
    rd = Chem.AddHs(rd, addCoords=True)

    # The added hydrogens carry no PDB residue metadata, so without this they
    # land in their own UNK residue and every later `resname CAU` selection
    # silently returns the heavy atoms only.
    reference = next((a.GetPDBResidueInfo() for a in rd.GetAtoms()
                      if a.GetPDBResidueInfo() is not None), None)
    counter = 0
    for atom in rd.GetAtoms():
        if atom.GetPDBResidueInfo() is not None:
            continue
        counter += 1
        stamp = Chem.AtomPDBResidueInfo()
        stamp.SetResidueName(reference.GetResidueName())
        stamp.SetResidueNumber(reference.GetResidueNumber())
        stamp.SetChainId(reference.GetChainId())
        stamp.SetIsHeteroAtom(True)
        stamp.SetName(f"H{counter}".ljust(4)[:4])
        atom.SetMonomerInfo(stamp)

    molecule = Molecule.from_rdkit(rd, allow_undefined_stereo=True)
    molecule.name = resname
    return molecule


SOLVENT_RESIDUES = {"HOH", "WAT", "NA", "CL"}
# Two heavy atoms of different molecules closer than this are overlapping rather
# than in contact. A hydrogen bond is 2.7 A between oxygens and the shortest real
# non-bonded heavy-atom contacts are around 2.4, so 2.2 catches only overlaps.
CLASH_CUTOFF_A = 2.2


def remove_periodic_clashes(modeller, *, cutoff_a: float = CLASH_CUTOFF_A) -> dict:
    """Delete solvent molecules that overlap another molecule's periodic image.

    `addMembrane` leaves these behind, and they are the reason the membrane runs
    would not start. Whole lipid molecules are kept when the patch is tiled, so
    the bilayer legitimately overhangs the periodic box by about 9 A on each
    side; `addSolvent` is then called to extend the box in z and fills the box
    interior, placing water where those overhanging atoms already sit once the
    boundary is applied.

    On the apo system this is 235 overlapping pairs, 231 of them water against
    water, with individual water oxygens carrying 200,000 kJ/mol/nm. Deleting
    the later molecule of each pair keeps the pre-equilibrated patch intact and
    costs 246 molecules of 42,022, about 1 per cent of the solvent, and takes
    the as-built energy from -2,056,143 to -2,213,385 kJ/mol.

    A correction to what this comment first claimed: these overlaps were *not*
    why the membrane runs would not start. That was the lipid z restraint, which
    was written non-periodically and saw the wrapped lower leaflet as a 117 A
    deviation. The +136,000,000 kJ/mol figure attributed here to single
    precision was that restraint, and the comparison that produced it was
    invalid: a run with no restraints on the CPU against a run with restraints
    on OpenCL. The repair is still worth doing, because the overlaps are real
    and the energy improves, but it is a tidy-up rather than the fix.

    A clash between two molecules that are not solvent is reported and not
    deleted: removing a lipid or part of the receptor is not a repair, and this
    build has never seen one that was not a disulphide.
    """
    import numpy as np
    from openmm import unit
    from scipy.spatial import cKDTree

    topology, positions = modeller.topology, modeller.positions
    box = topology.getUnitCellDimensions()
    if box is None:
        return {"removed": 0, "note": "no periodic box; nothing to check"}
    box = np.array(box.value_in_unit(unit.angstrom))

    coords = np.array(positions.value_in_unit(unit.angstrom))
    atoms = list(topology.atoms())
    heavy = np.array([i for i, a in enumerate(atoms)
                      if a.element is not None and a.element.symbol != "H"])
    wrapped = np.mod(coords[heavy], box)

    bonded = set()
    for bond in topology.bonds():
        first, second = bond[0].residue.index, bond[1].residue.index
        bonded.add((first, second))
        bonded.add((second, first))

    # boxsize makes the tree periodic, which is the whole point: these pairs are
    # only close once the boundary is applied and a plain tree does not see them.
    pairs = cKDTree(wrapped, boxsize=box).query_pairs(cutoff_a, output_type="ndarray")

    residues = list(topology.residues())
    victims: set[int] = set()
    structural: list[str] = []
    found = 0
    for i, j in pairs:
        one, two = atoms[heavy[i]], atoms[heavy[j]]
        first, second = one.residue.index, two.residue.index
        if first == second or (first, second) in bonded:
            continue
        found += 1
        solvent = [r for r in (first, second)
                   if residues[r].name in SOLVENT_RESIDUES]
        if not solvent:
            structural.append(
                f"{one.residue.name}{one.residue.id}:{one.name} against "
                f"{two.residue.name}{two.residue.id}:{two.name}")
            continue
        # Water before ions: deleting an ion would leave the system
        # non-neutral, and `addMembrane` has already chosen the counts that
        # neutralise it. Otherwise the later residue, which is the one
        # addSolvent added, so the pre-equilibrated patch is left as it was.
        waters = [r for r in solvent if residues[r].name in ("HOH", "WAT")]
        victims.add(max(waters) if waters else max(solvent))

    if victims:
        modeller.delete([residues[index] for index in sorted(victims)])
    return {"clashes": found, "removed": len(victims),
            "structural": structural[:5], "n_structural": len(structural)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--which", choices=("apo", "carazolol"), required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "md")
    args = parser.parse_args()

    from openmm import app, unit
    import openmm as mm

    start = ROOT / "data" / "systems" / f"b2ar_{args.which}_start.pdb"
    meta_path = ROOT / "data" / "systems" / f"b2ar_{args.which}_start.json"
    meta = json.loads(meta_path.read_text())

    check = check_orientation(start)
    log(f"orientation: belt |z| {check['belt_median']:.1f} A, "
        f"intracellular basics {check['basic_median']:.1f} A -> "
        f"{'ok' if check['passes'] else 'FAILED'}")
    if not check["passes"]:
        raise SystemExit("receptor is not in the membrane frame; refusing to embed")

    out_dir = Path(args.out) / f"b2ar_{args.which}"
    out_dir.mkdir(parents=True, exist_ok=True)

    from pdbfixer import PDBFixer

    log("PDBFixer: missing atoms and hydrogens at pH 7.4")
    fixer = PDBFixer(filename=str(start))
    fixer.findMissingResidues()
    # The construct is already trimmed to TM1-3 .. H8+3, so any missing residue
    # at an end is the trim rather than a gap, and building it back would undo
    # the construct the spec defined.
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # Taken out before the hydrogens go on, and merged back below. PDBFixer has
    # no "keep these heteroatoms" switch and cannot protonate a residue it has
    # no template for, so the ligand has to travel this way round.
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    # Section 5.1 names "Amber ff14SB + Lipid17 + TIP3P" as three things, but
    # OpenMM's amber14-all.xml already contains lipid17: loading it again fails
    # with "Residue template DLPC with the same override level 0 already
    # exists". The lipids are present either way, which a POPC count on the
    # built system confirms rather than assumes.
    forcefield_files = ["amber14-all.xml", "amber14/tip3p.xml"]
    generator = None
    ligand = None
    if meta.get("ligand"):
        from flexappeal import sources
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator

        code = meta["ligand"]
        sdf_path = ROOT / "data" / "systems" / f"{code}.sdf"
        if not sdf_path.is_file():
            log(f"fetching the chemical definition for {code} from the RCSB")
            sdf_path.write_bytes(sources.fetch_ccd_ideal(code).data)
        ligand = ligand_molecule(start, code, sdf_path)
        log(f"{code}: {ligand.n_atoms} atoms with hydrogens, formal charge "
            f"{ligand.total_charge.m:+.0f}")
        log(f"parameterising {code} with openff-2.2.1 (AM1-BCC, a minute or two)")
        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand], forcefield="openff-2.2.1")

    forcefield = app.ForceField(*forcefield_files)
    if generator is not None:
        forcefield.registerTemplateGenerator(generator.generator)

    def assemble():
        """The receptor plus its ligand, before any bilayer."""
        built = app.Modeller(fixer.topology, fixer.positions)
        if ligand is not None:
            built.add(ligand.to_topology().to_openmm(),
                      ligand.conformers[0].to_openmm())
        # gemmi writes a placeholder CRYST1 of 1 1 1 for a P 1 structure with no
        # cell, and `addMembrane` believes it: having a unit cell at all, it
        # takes the box height as max(patch z, cell z + 2 * padding) = the patch,
        # about 7 nm, rather than from the receptor, which is 9.4 nm tall with
        # ICL3 and H8 modelled. The receptor then wraps through the periodic
        # boundary onto itself and the relaxation produces NaN on its first
        # step. Dropping the cell makes `addMembrane` size z from the protein,
        # which is what it does when no cell is present.
        built.topology.setUnitCellDimensions(None)
        return built

    # `addMembrane` builds its own Context for the relaxation that inflates the
    # protein back into the patch, and with no platform given it takes the
    # fastest, which here is OpenCL. Apple's OpenCL is single precision only,
    # and that relaxation does not survive it: both systems died at
    # `integrator.step(20)` with "Particle coordinate is NaN" while the
    # assembly going in minimised cleanly to 8,286 kJ/mol. On the CPU platform,
    # at double precision, the same build completes. It also keeps the GPU free
    # for the production run, which is the concurrency section 2.3 asks for.
    relax_platform = mm.Platform.getPlatformByName("CPU")

    modeller = assemble()
    if ligand is not None:
        # openff writes chain and residue ids as ints where OpenMM uses strings,
        # and PDBFile.writeFile(keepIds=True) then fails on len(res.id).
        for chain in modeller.topology.chains():
            chain.id = str(chain.id)
        for residue in modeller.topology.residues():
            residue.id = str(residue.id)
        log(f"merged {code} back into the receptor: "
            f"{modeller.topology.getNumAtoms():,} atoms before the bilayer")
    padding = PADDING_NM
    log(f"addMembrane: {LIPID}, {padding} nm padding, {IONIC_STRENGTH_M} M NaCl")
    started = time.time()
    modeller.addMembrane(forcefield, lipidType=LIPID,
                         minimumPadding=padding * unit.nanometer,
                         ionicStrength=IONIC_STRENGTH_M * unit.molar,
                         positiveIon="Na+", negativeIon="Cl-",
                         platform=relax_platform)
    n_atoms = modeller.topology.getNumAtoms()
    log(f"built in {time.time() - started:.0f} s: {n_atoms:,} atoms")

    if n_atoms > MAX_ATOMS:
        log(f"over {MAX_ATOMS:,}; rebuilding at {FALLBACK_PADDING_NM} nm (section 5.1)")
        modeller = assemble()
        padding = FALLBACK_PADDING_NM
        modeller.addMembrane(forcefield, lipidType=LIPID,
                             minimumPadding=padding * unit.nanometer,
                             ionicStrength=IONIC_STRENGTH_M * unit.molar,
                             positiveIon="Na+", negativeIon="Cl-",
                             platform=relax_platform)
        n_atoms = modeller.topology.getNumAtoms()
        log(f"rebuilt: {n_atoms:,} atoms")

    lipids = sum(1 for r in modeller.topology.residues() if r.name in ("POPC", "POP"))
    waters = sum(1 for r in modeller.topology.residues() if r.name in ("HOH", "WAT"))
    inside = EXPECTED_ATOMS[0] <= n_atoms <= EXPECTED_ATOMS[1]
    log(f"{lipids:,} lipids, {waters:,} waters; expected "
        f"{EXPECTED_ATOMS[0]:,} to {EXPECTED_ATOMS[1]:,} atoms -> "
        f"{'within' if inside else 'outside, logged'}")

    repair = remove_periodic_clashes(modeller)
    log(f"periodic clash repair: {repair.get('clashes', 0)} overlapping pairs, "
        f"{repair.get('removed', 0)} solvent molecules deleted")
    if repair.get("n_structural"):
        log(f"  {repair['n_structural']} clashes between non-solvent molecules, "
            f"left in place: {'; '.join(repair['structural'])}")
    n_atoms = modeller.topology.getNumAtoms()

    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
        rigidWater=True, hydrogenMass=1.5 * unit.amu)

    (out_dir / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    with (out_dir / "solvated.pdb").open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    box = modeller.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    (out_dir / "build.json").write_text(json.dumps({
        "which": args.which, "state": meta["state"], "ligand": meta.get("ligand"),
        "atoms": n_atoms, "lipids": lipids, "waters": waters,
        "padding_nm": padding, "box_nm": [round(v, 2) for v in box],
        "within_expected": inside, "orientation": check,
        "clash_repair": repair,
        "forcefield": forcefield_files,
    }, indent=2))
    log(f"wrote {out_dir}/system.xml and solvated.pdb; box {[round(v,1) for v in box]} nm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
