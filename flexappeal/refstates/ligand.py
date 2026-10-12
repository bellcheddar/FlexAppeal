"""Ligand pose drift and contact persistence, for section 5.7's bound rows.

Section 5.7 asks the ligand-bound runs for "the same analysis plus ligand RMSD
and contact persistence". Two of the four production rows are bound,
`b2ar_carazolol` with carazolol and `abl1_imatinib` with imatinib, so this is
what those rows need beyond what every row already gets.

Both measurements are made **after superposing the frame on the same core the
state scores use**. That is the whole point of them: a ligand that holds its
pose is one that has not moved relative to the pocket, and a ligand that drifts
while the protein drifts with it has not come loose. Measuring in the
trajectory's own frame would report the protein's tumbling as ligand motion.

The pose RMSD is minimised over the ligand's own symmetry, reusing
`flexappeal.measure`'s graph automorphisms rather than growing a second copy.
Without it, a phenyl ring that flips by 180 degrees reads as one to two
angstroms of drift while being the same molecule in the same place, and over
50 ns that is the difference between reporting a stable pose and a lost one.
"""

from __future__ import annotations

from typing import Any

# Section 5.4's own contact definition for the R-spine is a heavy-atom
# distance, and the analysis tab counts a receptor residue as in contact with
# the ligand at 4.0 A of closest approach. The same number is used here so one
# build does not carry two definitions of "touching".
CONTACT_CUTOFF_A = 4.0

# A residue is reported as a contact when it is within the cutoff for at least
# this fraction of frames. Below it the residue is brushing past rather than
# holding on, and a table of every residue the ligand ever touched is a table
# nobody reads.
PERSISTENT_FRACTION = 0.5


def ligand_residues(topology) -> list[int]:
    """Residue indices that are neither protein nor anything `load` drops.

    `traj.kept_atoms` keeps the protein and whatever is not solvent, ion or
    bilayer, so what is left over in a kept trajectory is the ligand. Found by
    exclusion for the same reason `kept_atoms` is: a ligand this code has not
    been told about is exactly the case that has to work.
    """
    from .traj import NOT_THE_SYSTEM

    protein = set()
    for atom in topology.atoms:
        if atom.residue.is_protein:
            protein.add(atom.residue.index)
    out = []
    for residue in topology.residues:
        name = residue.name.upper()
        if residue.index in protein or name in NOT_THE_SYSTEM:
            continue
        if name.startswith("PO") and name != "PRO":
            continue
        out.append(residue.index)
    return out


def _heavy(topology, residue_index: int) -> list[int]:
    return [atom.index for atom in topology.residue(residue_index).atoms
            if atom.element is not None and atom.element.symbol != "H"]


def _imaged_ligand(traj, heavy: list[int]) -> tuple:
    """The ligand's coordinates in the periodic copy nearest the receptor.

    Returns `(xyz_in_nm, note)`. Without this the first real trajectory
    reported carazolol 110 Å from its own receptor and holding no contacts at
    all: the integrator had wrapped the ligand through the edge of the box, so
    it sat at y = 0.18 nm against a receptor at 12.83 in a box 13.11 tall. The
    two are touching. The arithmetic is not.

    It is the same minimum-image convention `write_lipid_bundle` applies to the
    bilayer, and the third time this build has met periodic wrapping: the
    viewer bundle had 17 of 50 frames in another image, and the lipids were
    shifted rather than imaged and grew a 269 Å span. A ligand is the case
    where it matters most, because the failure is silent and reads as a
    result: a pose RMSD of 5.4 Å and zero contacts is exactly what a ligand
    that had come out of its pocket would look like.

    One shift per frame for the whole molecule, taken from its centroid, not
    per atom. Imaging atoms independently would place a ligand sitting across
    the boundary into two pieces, which is a worse answer than the one being
    fixed and would not look wrong in the numbers.
    """
    import numpy as np

    xyz = traj.xyz[:, heavy, :]
    if traj.unitcell_lengths is None:
        return xyz, ("the trajectory carries no box, so the ligand could not "
                     "be imaged about the receptor")
    protein = [atom.index for atom in traj.topology.atoms
               if atom.residue.is_protein and atom.element is not None
               and atom.element.symbol != "H"]
    if len(protein) < 3:
        return xyz, "too little protein to image the ligand about"

    box = np.asarray(traj.unitcell_lengths)
    if not np.all(box > 0):
        return xyz, "the trajectory's box is degenerate, so no imaging was done"
    anchor = traj.xyz[:, protein, :].mean(axis=1)
    centre = xyz.mean(axis=1)
    shift = box * np.round((centre - anchor) / box)
    return xyz - shift[:, None, :], ""


def pose_and_contacts(traj, core_residues: list[int], *,
                      cutoff_a: float = CONTACT_CUTOFF_A) -> dict[str, Any]:
    """Pose RMSD per frame and per-residue contact persistence.

    Returns an empty dict for a run with no ligand, which is both apo rows, so
    the caller can store it unconditionally and the page can test for it.
    """
    import numpy as np

    from .traj import ca_coordinates, superpose, apply_transform

    found = ligand_residues(traj.topology)
    if not found:
        return {}
    # Grouped by residue name, not taken as the largest single residue.
    #
    # `abl1_imatinib`'s solvated system writes STI as *two* residues, 37 heavy
    # atoms and 31 hydrogens, so the largest-residue rule happened to return
    # imatinib's heavy atoms and happened to be right. It also handed the graph
    # inference a molecule with no hydrogens, which is why the symmetry mapping
    # came back empty. One molecule split across residues is one ligand.
    by_name: dict[str, list[int]] = {}
    for index in found:
        by_name.setdefault(traj.topology.residue(index).name.upper(),
                           []).append(index)
    name = max(by_name, key=lambda key: sum(
        traj.topology.residue(i).n_atoms for i in by_name[key]))
    members = by_name[name]
    heavy = [index for residue in members
             for index in _heavy(traj.topology, residue)]
    if len(heavy) < 3:
        return {"status": "ligand has fewer than three heavy atoms",
                "resname": name}

    core_coords, core_found = ca_coordinates(traj, core_residues)
    if len(core_found) < 3:
        return {"status": "core too small to superpose on", "resname": name}

    # Every frame onto the first, on the core. The reference is this run's own
    # first frame rather than a deposited pose: the question section 5.7 asks
    # is whether the ligand stayed where this simulation started it.
    reference = core_coords[0]
    # Imaged about the receptor before anything is measured. Both readings
    # depend on it: a wrapped ligand is a hundred angstrom from its own pocket,
    # which reads as a lost pose and no contacts rather than as a fault.
    imaged, imaging_note = _imaged_ligand(traj, heavy)
    ligand_xyz = imaged * 10.0
    moved = np.empty_like(ligand_xyz)
    for index in range(traj.n_frames):
        transform = superpose(core_coords[index], reference)
        moved[index] = apply_transform(ligand_xyz[index], *transform)

    symmetry, note = _symmetry_for(traj, members, heavy)
    drift = _drift(moved, symmetry)

    contacts = _contacts(traj, set(members), heavy, cutoff_a, imaged)
    out: dict[str, Any] = {
        "resname": name,
        "n_residues": len(members),
        "n_heavy_atoms": len(heavy),
        "pose_rmsd_a": [round(float(v), 3) for v in drift],
        "final_pose_rmsd_a": round(float(drift[-1]), 3),
        "max_pose_rmsd_a": round(float(drift.max()), 3),
        "symmetry_corrected": symmetry is not None,
        "contact_cutoff_a": cutoff_a,
        "persistent_contacts": contacts["persistent"],
        "n_contacts_per_frame": contacts["per_frame"],
        "mean_contacts": round(float(np.mean(contacts["per_frame"])), 2),
    }
    for extra in (note, imaging_note):
        if extra:
            out["note"] = (out["note"] + "; " + extra) if out.get("note") else extra
    return out


def _symmetry_for(traj, members: list[int], heavy: list[int]):
    """The ligand's heavy-atom automorphisms, or None with a reason.

    Reused from `flexappeal.measure` rather than reimplemented. A note is
    returned rather than a silent fallback, because an uncorrected pose RMSD on
    a symmetric ligand overstates drift and the reader has to be able to tell
    which they are looking at.
    """
    import numpy as np

    try:
        from ..measure import heavy_automorphisms, ligand_mol
    except Exception:                                           # noqa: BLE001
        return None, "the symmetry helpers could not be imported"

    atoms = np.array(sorted(atom.index for residue in members
                            for atom in traj.topology.residue(residue).atoms),
                     dtype=np.int64)
    try:
        mol = ligand_mol(traj, atoms)
        symmetry = heavy_automorphisms(mol, len(atoms))
    except Exception:                                           # noqa: BLE001
        return None, ("the ligand's graph could not be inferred, so the pose "
                      "RMSD matches atoms by index and overstates drift for a "
                      "symmetric ligand")
    if symmetry is None:
        return None, ("no symmetry mapping was available, so the pose RMSD "
                      "matches atoms by index")

    # The mappings index the molecule's heavy atoms; `heavy` indexes the
    # topology's. Both are ascending, so they must be the same list. Checked
    # rather than assumed: a mapping applied to a differently ordered array
    # silently measures the wrong thing, which is the fault `measure` records
    # at the same point.
    in_mol = [slot for slot, index in enumerate(atoms) if int(index) in set(heavy)]
    if len(in_mol) != len(heavy) or not np.array_equal(
            atoms[np.array(in_mol, dtype=np.int64)],
            np.asarray(heavy, dtype=np.int64)):
        return None, ("the ligand's inferred graph did not line up with its "
                      "heavy atoms, so the pose RMSD matches atoms by index")
    return symmetry, ""


def _drift(coords, symmetry):
    """Per-frame RMSD to frame 0, minimised over the ligand's symmetry."""
    import numpy as np

    try:
        from ..measure import _symmetric_drift

        return _symmetric_drift(coords, symmetry)
    except Exception:                                           # noqa: BLE001
        reference = coords[0]
        return np.sqrt(((coords - reference) ** 2).sum(axis=2).mean(axis=1))


def _contacts(traj, members: set, heavy: list[int],
              cutoff_a: float, ligand_xyz=None) -> dict[str, Any]:
    """Which receptor residues the ligand holds on to, and for how much of the run.

    Closest heavy-atom approach per residue, so a residue counts once however
    many of its atoms are near, which is what makes the number a count of
    contacts rather than of atom pairs.
    """
    import numpy as np

    cutoff_nm = cutoff_a / 10.0
    candidates = []
    for residue in traj.topology.residues:
        if residue.index in members or not residue.is_protein:
            continue
        atoms = [a.index for a in residue.atoms
                 if a.element is not None and a.element.symbol != "H"]
        if atoms:
            candidates.append((residue, np.array(atoms, dtype=np.int64)))

    # The imaged coordinates, not the raw ones. `ligand_xyz` is already in the
    # periodic copy nearest the receptor; falling back to the raw slice keeps
    # the function usable on its own, which the tests rely on.
    coords = (traj.xyz[:, np.array(heavy, dtype=np.int64), :]
              if ligand_xyz is None else ligand_xyz)
    per_frame = np.zeros(traj.n_frames, dtype=int)
    persistent = []
    for residue, atoms in candidates:
        # (frames, ligand heavy, residue heavy) is the honest shape and is too
        # large to hold for a whole protein, so this is done residue by
        # residue and reduced immediately.
        diff = (coords[:, :, None, :] - traj.xyz[:, None, atoms, :])
        closest = np.sqrt((diff ** 2).sum(axis=3)).min(axis=(1, 2))
        touching = closest <= cutoff_nm
        per_frame += touching.astype(int)
        fraction = float(touching.mean())
        if fraction >= PERSISTENT_FRACTION:
            persistent.append({
                "residue": int(residue.resSeq),
                "name": residue.name,
                "fraction_of_frames": round(fraction, 3),
                "closest_a": round(float(closest.min() * 10.0), 2),
            })
    persistent.sort(key=lambda row: (-row["fraction_of_frames"],
                                     row["residue"]))
    return {"persistent": persistent,
            "per_frame": [int(v) for v in per_frame]}
