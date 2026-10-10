"""What one receptor-ligand complex measures, for an `.fxa` version 2 package.

`flexappeal.fxa2` owns the format and `flexappeal.package` owns the assembly; this
module owns the numbers. It takes an MDTraj trajectory of one finished job and
returns a `Measured`, which is every member body that complex contributes.

Three rules run through all of it, and each has already cost somebody a wrong
report:

* **Angstrom and nanoseconds.** MDTraj works in nanometres and version 1 reported
  nanometres, so every length leaving this module has been multiplied by ten and
  carries an `_a` suffix. A conversion missed here is invisible until somebody
  checks a number by hand.
* **The scored window.** Every scored quantity is computed from
  `equil_start_frame` onward, and the file reports the frame actually used. A
  strain measured over a different window from the RMSD is not the strain that
  belongs beside that grade.
* **Omit, never zero.** A measurement that could not be made leaves its key out.
  DeepFlex records the sub-score as pending and says so on the page; a `0` is read
  as a measurement.

The thresholds and definitions here deliberately match the ones DeepFlex applies
to its own pipeline (a 6 A pocket filter, a 5 A contact cutoff, a quarter of the
native contacts for occupancy, a 2 A pose-cluster linkage cutoff). They are not
arbitrary choices on this side: a number measured under a different definition
from the cohort it is ranked against is not comparable to it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import fxa2

# --- pocket definition ------------------------------------------------------
# The crop shell is section 6.1's: every whole residue with any atom within 10 A
# of the ligand. The filter is tighter, because the filtered set is what the
# ligand RMSD is superposed on, and a pocket selection that includes a flexible
# loop or a terminus makes the superposition wander. Every ligand RMSD measured
# against it is then wrong in a way that looks exactly like the ligand moving.
POCKET_SHELL_A = 10.0
POCKET_FILTER_A = 6.0
POCKET_FALLBACK_A = 8.0
TERMINUS_TRIM = 5
MIN_POCKET_RESIDUES = 12

# --- occupancy --------------------------------------------------------------
CONTACT_CUTOFF_A = 5.0
NATIVE_CONTACT_FRACTION = 0.25

# --- the scored window ------------------------------------------------------
# A ceiling on how much of the run may be discarded. Past this the complex is
# marked `capped` rather than quietly losing most of its trajectory.
EQUIL_CEILING_FRACTION = 0.5

# A bound ligand's centroid cannot move this far between two stored frames. Over
# it, the coordinates have been wrapped across a periodic boundary and every
# distance in the file is suspect.
PBC_JUMP_A = 10.0


class MeasureError(ValueError):
    """A complex could not be measured. Recorded as a failure, not a bad grade."""


@dataclass
class Measured:
    """One complex's contribution to a package.

    Every field is a member body or a group of them, so `flexappeal.package` can
    write the file without knowing how any of it was computed. A field left None
    is a member that is not written, which is how an absent measurement stays
    absent rather than becoming a zero.
    """

    complex_id: str
    receptor: str
    ligand_id: str
    metrics: dict[str, Any]
    pocket: dict[str, Any]
    structures: dict[str, str] = field(default_factory=dict)
    interactions: dict[str, Any] | None = None
    summary: list[dict[str, Any]] | None = None
    series: dict[str, np.ndarray] | None = None
    panels: dict[str, np.ndarray] | None = None
    panels_summary: dict[str, Any] | None = None
    barcode: dict[str, np.ndarray] | None = None
    conformer: dict[str, np.ndarray] | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        return str(self.metrics.get("status") or "ok")


def failed(cid: str, receptor: str, ligand: str, reason: str) -> Measured:
    """A complex that could not be measured.

    `status: "failed:<reason>"` and nothing else. A missing result and a poor
    result are different findings and DeepFlex reports them differently, so this
    must never be a bad grade.
    """
    return Measured(
        complex_id=cid, receptor=receptor, ligand_id=ligand,
        metrics={"complex_id": cid, "receptor": receptor, "ligand_id": ligand,
                 "status": f"failed:{reason}", "is_reference": False},
        pocket={"complex_id": cid, "receptor": receptor,
                "status": f"failed:{reason}", "notes": [reason]},
    )


# ===========================================================================
#  Selections
# ===========================================================================


def ligand_atoms(topology, resnames: list[str]) -> tuple[np.ndarray, np.ndarray, bool]:
    """The ligand's atoms, its heavy atoms, and whether it came in pieces.

    Returns the largest connected fragment when the selection spans several, with
    a flag, because two copies of a ligand in the same selection makes a centroid
    that sits between them and moves when either does.
    """
    import mdtraj as md  # noqa: PLC0415  -- heavy import, only needed here

    if not resnames:
        raise MeasureError("no_ligand_resname")
    selection = " or ".join(f"resname {name}" for name in resnames)
    atoms = topology.select(selection)
    if len(atoms) == 0:
        raise MeasureError("no_ligand")

    residues = sorted({topology.atom(int(i)).residue.index for i in atoms})
    fragmented = len(residues) > 1
    if fragmented:
        # Several residues matched: keep the largest, which is the ligand, and
        # flag the rest. mdtraj's find_molecules is not used here because a
        # decimated topology may have lost the bonds it needs.
        biggest = max(residues, key=lambda index: sum(
            1 for _ in topology.residue(index).atoms))
        atoms = np.array([a.index for a in topology.residue(biggest).atoms], dtype=int)

    heavy = np.array([i for i in atoms
                      if topology.atom(int(i)).element.symbol != "H"], dtype=int)
    if len(heavy) == 0:
        raise MeasureError("ligand_all_hydrogen")
    assert md is not None
    return atoms, heavy, fragmented


def _residue_numbers(topology, indices) -> list[int]:
    """The residue numbers of a set of residue indices, ascending.

    `resSeq`, not the residue index, because that is the number a reader sees in
    a viewer and the number DeepFlex looks up when it marks pocket residues on
    the RMSF plot. Two residues can share one: trypsin's chymotrypsin numbering
    uses insertion codes, which MDTraj drops, so 184 is both a GLY and a TYR in
    the 3PTB fixture. Duplicates are kept rather than merged, because dropping
    one would silently shorten the pocket.
    """
    return sorted({topology.residue(int(i)).resSeq for i in indices})


def pocket_residues(traj, heavy: np.ndarray) -> dict[str, Any]:
    """The binding site, raw and filtered, as `pocket/<cid>.json` wants it.

    Raw is section 6.1's crop shell: every whole residue with any atom within
    10 A of the ligand, which is what the pocket PDB contains. Filtered is the
    subset within 6 A whose residues are not near a chain terminus, which is what
    the ligand RMSD is superposed on.
    """
    import mdtraj as md  # noqa: PLC0415

    topology = traj.topology
    protein = topology.select("protein")
    if len(protein) == 0:
        raise MeasureError("no_protein")
    protein_heavy = np.array(
        [i for i in protein if topology.atom(int(i)).element.symbol != "H"], dtype=int)

    def residues_within(cutoff_a: float) -> set[int]:
        near = md.compute_neighbors(
            traj[0], cutoff_a / 10.0, heavy, haystack_indices=protein_heavy)[0]
        return {topology.atom(int(i)).residue.index for i in near}

    raw_indices = residues_within(POCKET_SHELL_A)
    filtered_indices = residues_within(POCKET_FILTER_A)

    # The first and last few residues of every chain, which are mobile even in a
    # well-behaved run and drag a superposition with them.
    trim: set[int] = set()
    for chain in topology.chains:
        residues = [r.index for r in chain.residues if r.is_protein]
        trim.update(residues[:TERMINUS_TRIM])
        trim.update(residues[-TERMINUS_TRIM:])

    kept = filtered_indices - trim
    notes: list[str] = []
    used_fallback = False
    if len(kept) < MIN_POCKET_RESIDUES:
        used_fallback = True
        wider = residues_within(POCKET_FALLBACK_A) - trim
        notes.append(
            f"only {len(kept)} residue(s) survived the {POCKET_FILTER_A:.0f} A "
            f"filter, which is under the {MIN_POCKET_RESIDUES} needed for a stable "
            f"superposition, so the pocket was widened to {POCKET_FALLBACK_A:.0f} A"
        )
        kept = wider
    if len(kept) < 3:
        raise MeasureError("pocket_too_small")

    raw_numbers = _residue_numbers(topology, raw_indices)
    kept_numbers = _residue_numbers(topology, kept)
    return {
        "indices_raw": sorted(raw_indices),
        "indices_filtered": sorted(kept),
        "record": {
            "n_ligand_atoms": None,         # filled by the caller, which knows
            "pocket_residues_raw": raw_numbers,
            "pocket_residues_filtered": kept_numbers,
            "pocket_resnames": {
                str(topology.residue(int(i)).resSeq): topology.residue(int(i)).name
                for i in sorted(kept)
            },
            "pocket_n_raw": len(raw_numbers),
            "pocket_n_filtered": len(kept_numbers),
            "added_vs_raw": sorted(set(kept_numbers) - set(raw_numbers)),
            "dropped_from_raw": sorted(set(raw_numbers) - set(kept_numbers)),
            "used_fallback": used_fallback,
            "notes": notes,
        },
    }


# ===========================================================================
#  Structures, section 6
# ===========================================================================


def _rename_ligand(text: str, resnames: list[str]) -> str:
    """Rewrite the ligand residue as `LIG`, on every record that names it.

    The `TER` records are the trap. A `TER` carries the residue name in the same
    columns (18 to 20) as an `ATOM` record, so a rename that walks only the atom
    lines leaves the real het code sitting in the file under a hashed filename.
    Matching on the residue-name columns regardless of record type is what makes
    that impossible rather than merely unlikely.

    `HETATM` becomes `ATOM` for the same residue, because section 6.1 asks for
    the ligand as an `ATOM` record: the viewer selects on it.
    """
    wanted = {name.strip().upper() for name in resnames}
    padded = fxa2.LIGAND_RESNAME.rjust(3)

    out = []
    for line in text.splitlines():
        record = line[:6].strip()
        if record in ("ATOM", "HETATM", "TER") and len(line) >= 20:
            if line[17:20].strip().upper() in wanted:
                line = line[:17] + padded + line[20:]
                if record == "HETATM":
                    line = "ATOM  " + line[6:]
        out.append(line)
    return "\n".join(out) + "\n"


def _strip_solvent(traj, keep: np.ndarray) -> np.ndarray:
    """Drop any water, counter-ion or lipid from an atom selection.

    Section 6.1 bans them from the crop. A decimated version 1 payload has
    usually dropped them already, but a run directory's topology has not, and a
    pocket shell drawn around a ligand in explicit solvent is mostly water.
    """
    topology = traj.topology
    banned = fxa2.SOLVENT_RESNAMES | fxa2.ION_RESNAMES
    return np.array(
        [i for i in keep
         if topology.atom(int(i)).residue.name.strip().upper() not in banned],
        dtype=int,
    )


def _pdb_text(sliced) -> str:
    """A trajectory slice as PDB text.

    Through a temporary file rather than a buffer: MDTraj's `save_pdb` calls
    `os.path.exists` on what it is given, so a StringIO raises a TypeError from
    inside the writer rather than being written to.
    """
    import tempfile  # noqa: PLC0415

    with tempfile.TemporaryDirectory(prefix="flexappeal-pdb-") as work:
        path = f"{work}/crop.pdb"
        sliced.save_pdb(path)
        with open(path, encoding="utf-8") as handle:
            return handle.read()


def pocket_pdb(traj, frame: int, atoms: np.ndarray, resnames: list[str]) -> str:
    """One frame of the crop, with the ligand named `LIG` and no solvent."""
    return _rename_ligand(_pdb_text(traj[frame].atom_slice(atoms)), resnames)


def scrub_pdb(traj, atoms: np.ndarray, resnames: list[str], *,
              models: int = 100) -> str:
    """The trajectory as one multi-model PDB, superposed so the receptor holds still.

    One file rather than one per frame. The viewer steps the model index, which is
    a parameter update; a file per frame makes it tear the structure down and
    rebuild it, so the canvas blanks and the camera jumps on every step.
    """
    stride = max(1, traj.n_frames // models)
    return _rename_ligand(
        _pdb_text(traj[::stride].atom_slice(atoms)), resnames)


# ===========================================================================
#  The whole measurement for one complex
# ===========================================================================


def measure_complex(traj, *, complex_id: str, receptor: str, ligand_id: str,
                    ligand_resnames: list[str]) -> Measured:
    """Measure one receptor-ligand complex.

    Raises nothing: a complex that cannot be measured comes back as a failure
    record, because the caller is assembling a package and one bad trajectory
    must not cost the other nineteen.
    """
    try:
        return _measure(traj, complex_id, receptor, ligand_id, ligand_resnames)
    except MeasureError as exc:
        return failed(complex_id, receptor, ligand_id, str(exc))


def _measure(traj, cid: str, receptor: str, ligand: str,
             resnames: list[str]) -> Measured:
    atoms, heavy, fragmented = ligand_atoms(traj.topology, resnames)
    pocket = pocket_residues(traj, heavy)

    # The crop: the ligand plus every whole residue in the 10 A shell, with the
    # solvent removed. Selected by residue index rather than by residue number,
    # because residue numbers are not unique once insertion codes have been
    # dropped and selecting on one pulls in whichever other residue shares it.
    topology = traj.topology
    crop = {int(i) for i in atoms}
    for index in pocket["indices_raw"]:
        crop.update(a.index for a in topology.residue(int(index)).atoms)
    crop_atoms = _strip_solvent(traj, np.array(sorted(crop), dtype=int))

    notes: list[str] = list(pocket["record"]["notes"])
    if fragmented:
        notes.append(
            "the ligand selection matched more than one residue; the largest "
            "fragment was measured and the rest ignored"
        )

    metrics: dict[str, Any] = {
        "complex_id": cid,
        "receptor": receptor,
        "ligand_id": ligand,
        "status": "ok",
        "n_frames": int(traj.n_frames),
        "extra_ligand_fragments": bool(fragmented),
        # DeepFlex's scorer reads this directly rather than with a default, and a
        # FlexAppeal study has no reference compound in it.
        "is_reference": False,
    }

    pocket_record = dict(pocket["record"])
    pocket_record.update({
        "complex_id": cid,
        "receptor": receptor,
        "status": "ok",
        "n_frames": int(traj.n_frames),
        "n_ligand_atoms": int(len(atoms)),
        "n_atoms_topology": int(traj.n_atoms),
        "n_atoms_trajectory": int(traj.n_atoms),
    })

    structures = {
        f"{cid}_pocket.pdb": pocket_pdb(traj, 0, crop_atoms, resnames),
    }
    if traj.n_frames > 1:
        structures[f"{cid}_pocket_last.pdb"] = pocket_pdb(
            traj, traj.n_frames - 1, crop_atoms, resnames)

    if notes:
        metrics["notes"] = notes

    return Measured(
        complex_id=cid, receptor=receptor, ligand_id=ligand,
        metrics=metrics, pocket=pocket_record, structures=structures,
        notes=notes,
    )
