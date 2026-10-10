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

# --- convergence and pose clusters ------------------------------------------
# Linkage distance for the pose clusters, in angstrom, and the histogram
# resolution for the block-split divergence. DeepFlex's own numbers again, so a
# cluster count and a divergence here mean what they mean there.
POSE_CLUSTER_CUTOFF_A = 2.0
CONV_BINS = 20
CONV_BLOCKS = 4
# Scored frames per histogram bin. A finite-sample histogram divergence is biased
# upward by about the bin count over the sample size, so a fixed bin count makes
# the bias depend on the run length and a short, perfectly converged run scores
# as a diverging one. This divisor holds the bias constant and puts a long run on
# CONV_BINS bins, which is the resolution the spec's thresholds were calibrated
# against.
CONV_FRAMES_PER_BIN = 90
# Fewer than this many frames per block and the four histograms are too sparse to
# compare: the divergence then measures the sampling rather than the trajectory.
MIN_FRAMES_PER_BLOCK = 6
# The pairwise distance matrix is this many frames on a side at most. It is
# quadratic in the frame count, and a cluster assignment does not improve past a
# few hundred samples of the same window.
MAX_CLUSTER_FRAMES = 400

# --- ligand strain ----------------------------------------------------------
# Frames of the scored window sampled for the strain and the torsions. Strain is
# a slowly varying quantity and each frame costs two MMFF minimisations, so a
# finer sampling buys noise and a bigger member. The spec asks for around this
# many for the panel.
MAX_STRAIN_FRAMES = 120
MAX_STRAIN_ITERATIONS = 500
# Force constant on the torsion restraints during the held minimisation, in
# kcal/mol/rad^2. Stiff enough that a dihedral does not leave its well while
# bonds and angles relax, soft enough not to fight the minimiser into
# non-convergence.
TORSION_RESTRAINT_K = 1000.0
# Rotamer binning: three wells at -180..-60, -60..60, 60..180. A dihedral has to
# clear a boundary by the hysteresis before the well it is assigned to changes,
# or a torsion rattling on a boundary registers a flip on every frame.
WELL_EDGES = (-60.0, 60.0)
HYSTERESIS_DEG = 20.0
# Total charges tried when inferring bond orders, in order. Benzamidine is an
# amidinium at pH 7 and a neutral amidine on paper, and nothing in a trajectory
# says which was simulated, so the first assignment that sanitises is taken.
CANDIDATE_CHARGES = (0, 1, -1, 2, -2)

# --- receptor panels --------------------------------------------------------
# Both matrices are reduced to at most this on a side before being written. A
# 2,000-frame 2D RMSD is four million floats and the panel is a few hundred
# pixels wide, so the reduction belongs here rather than at the far end.
MAX_PANEL_SIDE = fxa2.MAX_MATRIX_SIDE


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
    conformer_summary: dict[str, Any] | None = None
    depiction_svg: str | None = None
    depiction_painted: bool = False
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

    Section 6.1 bans all three from the crop. `flexappeal.package` already
    reduces a source to protein and ligand before measuring it, so this is the
    second line rather than the first: a version 1 `.fxa` written at the `full`
    payload tier keeps every atom, and a 10 A shell drawn around a ligand in a
    bilayer is mostly lipid.
    """
    topology = traj.topology
    banned = fxa2.SOLVENT_RESNAMES | fxa2.ION_RESNAMES | fxa2.LIPID_RESNAMES
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
#  The scored window, section 4.1
# ===========================================================================


def statistical_inefficiency(series: np.ndarray) -> float:
    """The number of correlated samples per independent one.

    The normalised autocorrelation summed until it first turns negative, which is
    the standard truncation: past that point the estimate is noise about zero, so
    continuing to sum it adds variance rather than information.
    """
    n = len(series)
    if n < 3:
        return 1.0
    deviation = series - series.mean()
    variance = float((deviation * deviation).mean())
    if variance <= 0.0:
        return 1.0

    g = 1.0
    for lag in range(1, n - 1):
        correlation = float((deviation[: n - lag] * deviation[lag:]).mean()) / variance
        if correlation <= 0.0 and lag > 1:
            break
        g += 2.0 * correlation * (1.0 - lag / n)
    return max(1.0, g)


def detect_equilibration(series: np.ndarray, *,
                         candidates: int = 50) -> tuple[int, str]:
    """The first frame to score from, and how it was chosen.

    Chodera's rule: take the truncation that leaves the most *effectively
    uncorrelated* samples, `(N - t) / g(t)`. Discarding more frames lowers the
    sample count and also lowers the correlation, and this is where those two
    trade off. It is pymbar's `detect_equilibration` without the dependency
    rather than a cheaper substitute for it: the same quantity is maximised.

    Evaluated at about fifty candidate truncations rather than at every frame.
    The curve is smooth, and the search is otherwise cubic in the frame count,
    which on a 2,000-frame run costs minutes to produce one frame index.

    Capped at half the run. Over the cap the method is reported as `ceiling`,
    which DeepFlex renders as `capped` on the page: honest rather than hidden.
    """
    n = len(series)
    if n < 8:
        # Too short to tell equilibrated from not. Score all of it and say so,
        # rather than discard a fifth of eight frames on no evidence.
        return 0, "all_frames"

    # Searched unconstrained and clamped afterwards, rather than searched inside
    # the cap. Clamping a search that was already bounded cannot tell a run that
    # equilibrated just under the cap from one that never equilibrated at all,
    # and those two deserve different words on the page.
    step = max(1, (n - 8) // candidates)
    best_start, best_effective = 0, -1.0
    for start in range(0, n - 4, step):
        tail = series[start:]
        if len(tail) < 4:
            break
        effective = len(tail) / statistical_inefficiency(tail)
        if effective > best_effective:
            best_effective, best_start = effective, start

    ceiling = max(1, int(n * EQUIL_CEILING_FRACTION))
    if best_start > ceiling:
        return ceiling, "detect_equilibration->ceiling"
    return best_start, "detect_equilibration"


# ===========================================================================
#  Geometry
# ===========================================================================


def _copy(traj, frames: slice | None = None):
    """An independent copy, because `superpose` rewrites coordinates in place.

    Three different superpositions are taken below and each has to fit against
    the original frame 0. Fitting one on top of another gives numbers that look
    plausible and are not what their names say, and MDTraj's slicing does not
    reliably hand back a copy.
    """
    import mdtraj as md  # noqa: PLC0415

    window = slice(None) if frames is None else frames
    return md.Trajectory(
        xyz=traj.xyz[window].copy(), topology=traj.topology, time=traj.time[window],
        unitcell_lengths=(None if traj.unitcell_lengths is None
                          else traj.unitcell_lengths[window].copy()),
        unitcell_angles=(None if traj.unitcell_angles is None
                         else traj.unitcell_angles[window].copy()),
    )


def _fit(traj, fit_on: np.ndarray):
    """A copy superposed on frame 0 over the given atoms.

    The reference is a *separate* object, not the trajectory being fitted.
    MDTraj's `superpose` rewrites coordinates in place, so passing a trajectory
    as its own reference moves frame 0 first and then fits every later frame
    against the moved version. A rigid 2 A translation came out as 1.88 A that
    way: wrong, plausible, and a unit test away from shipping.

    One frame of reference, not a whole copy of it. `superpose` reads only
    `reference.xyz[frame]`, so a full second copy doubled the peak memory of
    every fit for nothing; on a long run that is the difference between a
    gigabyte and two. Checked rather than assumed: the fitted RMSDs are
    bit-identical either way.
    """
    reference = _copy(traj, slice(0, 1))
    fitted = _copy(traj)
    fitted.superpose(reference, frame=0, atom_indices=fit_on)
    return fitted


def _drift_rmsd(fitted, measure_on: np.ndarray) -> np.ndarray:
    """Per-frame RMSD of an atom set to its own frame-0 position, in angstrom.

    Computed by hand rather than with `md.rmsd`, which superposes on exactly the
    atoms it then measures. The whole point of the pose RMSD is that the fit and
    the measurement are different sets: fit the pocket, measure the ligand, and
    what comes out is the ligand moving inside a receptor held still.
    """
    coords = fitted.xyz[:, measure_on, :] * 10.0
    return np.sqrt(((coords - coords[0]) ** 2).sum(axis=2).mean(axis=1))


def _atoms_of(topology, residue_indices, *, name: str | None = None,
              heavy: bool = False) -> np.ndarray:
    """Atoms of the given protein residues, optionally one atom name or heavy only."""
    out: list[int] = []
    for index in residue_indices:
        residue = topology.residue(int(index))
        if not residue.is_protein:
            continue
        for atom in residue.atoms:
            if name is not None and atom.name != name:
                continue
            if heavy and atom.element.symbol == "H":
                continue
            out.append(atom.index)
    return np.array(sorted(out), dtype=int)


def time_axis(traj, frame_spacing_ps: float | None) -> tuple[np.ndarray, str | None]:
    """Nanoseconds per frame, from the trajectory or from the configured spacing.

    An `.xtc` carries its own times and they are the right answer. Some writers
    leave them at zero, and a time axis of all zeros collapses every chart to a
    single point without anything raising, so the configured frame spacing is the
    fallback and the frame index is the last resort. Returns the axis and a note
    when it is not the trajectory's own.
    """
    time_ps = np.asarray(traj.time, dtype=np.float64)
    if len(time_ps) > 1 and np.all(np.diff(time_ps) > 0):
        return time_ps / 1000.0, None
    if frame_spacing_ps:
        return (np.arange(traj.n_frames, dtype=np.float64) * float(frame_spacing_ps)
                / 1000.0), (
            f"the trajectory carries no usable time stamps, so the time axis is "
            f"built from the configured frame spacing of {frame_spacing_ps:g} ps")
    return np.arange(traj.n_frames, dtype=np.float64), (
        "the trajectory carries no usable time stamps and no frame spacing was "
        "configured, so every time here is a frame index and not nanoseconds")


def _slope_per_ns(values: np.ndarray, time_ns: np.ndarray) -> float | None:
    """A least-squares slope, or None when there is nothing to fit."""
    if len(values) < 3:
        return None
    good = np.isfinite(values) & np.isfinite(time_ns)
    if good.sum() < 3 or np.ptp(time_ns[good]) <= 0:
        return None
    return float(np.polyfit(time_ns[good], values[good], 1)[0])


def residue_contacts(traj, ligand_residue: int,
                     residue_indices: list[int]) -> np.ndarray:
    """Per-frame count of receptor residues in contact with the ligand.

    Residues rather than atom pairs, which is what section 4 asks for
    (`n_native_contacts`: "receptor residues within the contact cutoff of the
    ligand"). Closest-heavy-atom distance per residue, under the minimum image
    convention, so a residue counts once however many of its atoms are near.
    """
    import mdtraj as md  # noqa: PLC0415

    pairs = np.array([[ligand_residue, int(index)] for index in residue_indices],
                     dtype=int)
    if not len(pairs):
        raise MeasureError("no_pocket_residues")
    distances, _ = md.compute_contacts(traj, contacts=pairs, scheme="closest-heavy",
                                       periodic=traj.unitcell_vectors is not None,
                                       ignore_nonprotein=False)
    return (distances <= CONTACT_CUTOFF_A / 10.0).sum(axis=1).astype(np.int64)


# ===========================================================================
#  Interactions, section 5.3, 5.5 and the S_int inputs
# ===========================================================================

# Typed by geometry and element, with no chemistry toolkit involved. FlexAppeal
# has no bond orders, no formal charges and no ring perception for a ligand at
# this point: the SDF the bundle fetched stays on the user's machine and never
# enters the package, and re-deriving chemistry from coordinates alone is how a
# fingerprint acquires types it cannot justify.
#
# So the three types here are the ones that heavy-atom geometry supports, and
# nothing claims to be a pi-stack. That is a real loss against a ProLIF
# fingerprint and it is recorded in the interaction record rather than glossed:
# the occupancies are honest, the type vocabulary is coarser.
#
# `HBond` is the standard heavy-atom hydrogen-bond proxy, a polar contact. The
# angle term is deliberately absent: a decimated payload may carry no hydrogens
# at all, a ligand's hydrogens have no template in an MDTraj topology, and a
# criterion that silently changes definition depending on which source the
# trajectory came from is worse than one that is coarse everywhere.
HBOND_CUTOFF_A = 3.5
IONIC_CUTOFF_A = 4.0
HYDROPHOBIC_CUTOFF_A = 4.5

POLAR_ELEMENTS = frozenset({"N", "O", "S"})
HYDROPHOBIC_RESIDUES = frozenset({
    "ALA", "VAL", "LEU", "ILE", "PHE", "MET", "PRO", "TRP", "CYS", "TYR",
})
# Side-chain atoms that carry a formal charge at physiological pH.
ANIONIC_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATIONIC_ATOMS = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
                  ("HIS", "ND1"), ("HIS", "NE2"), ("HIP", "ND1"), ("HIP", "NE2")}

# An interaction holding this fraction of the scored window is treated as key.
# DeepFlex's own threshold, so the two agree on which contacts carry the grade.
KEY_INTERACTION_MIN_OCCUPANCY = 0.3

# The fingerprint is sampled rather than computed on every frame of a long run.
# The barcode is capped at 400 columns anyway, and an occupancy over 600 frames
# and one over 2,000 are the same number to more digits than anybody reads.
MAX_FINGERPRINT_FRAMES = 600
# Distances are thresholded in blocks of pairs so that a large ligand against a
# large pocket never materialises a frames-by-pairs float array.
PAIR_BLOCK = 4096


def _candidate_pairs(topology, ligand_heavy: np.ndarray,
                     residue_indices) -> tuple[np.ndarray, list[tuple[int, str]], dict]:
    """Every (ligand atom, receptor atom) pair that could be a typed interaction.

    Returns the pairs, the (residue index, type) each one would contribute to,
    and the cutoff per type. Built once from the topology, because which type a
    pair could be is a property of the atoms rather than of the frame.
    """
    pairs: list[tuple[int, int]] = []
    labels: list[tuple[int, str]] = []
    cutoffs: dict[str, float] = {
        "HBond": HBOND_CUTOFF_A,
        "Ionic": IONIC_CUTOFF_A,
        "Hydrophobic": HYDROPHOBIC_CUTOFF_A,
    }

    ligand = [(int(i), topology.atom(int(i)).element.symbol) for i in ligand_heavy]
    for index in residue_indices:
        residue = topology.residue(int(index))
        if not residue.is_protein:
            continue
        name = residue.name.strip().upper()
        hydrophobic_residue = name in HYDROPHOBIC_RESIDUES
        for atom in residue.atoms:
            element = atom.element.symbol
            if element == "H":
                continue
            anionic = (name, atom.name) in ANIONIC_ATOMS
            cationic = (name, atom.name) in CATIONIC_ATOMS
            for ligand_atom, ligand_element in ligand:
                if element in POLAR_ELEMENTS and ligand_element in POLAR_ELEMENTS:
                    pairs.append((ligand_atom, atom.index))
                    labels.append((index, "HBond"))
                if (anionic and ligand_element == "N") or \
                        (cationic and ligand_element == "O"):
                    pairs.append((ligand_atom, atom.index))
                    labels.append((index, "Ionic"))
                if hydrophobic_residue and element == "C" and ligand_element == "C":
                    pairs.append((ligand_atom, atom.index))
                    labels.append((index, "Hydrophobic"))

    return np.array(pairs, dtype=int).reshape(-1, 2), labels, cutoffs


def fingerprint(traj, ligand_heavy: np.ndarray, residue_indices
                ) -> tuple[np.ndarray, list[tuple[int, str]], int]:
    """Per-frame interaction bits, one column per (residue, type).

    Returns the bit matrix (frames by interactions), the (residue index, type)
    of each column, and the frame stride used.
    """
    import mdtraj as md  # noqa: PLC0415

    pairs, labels, cutoffs = _candidate_pairs(traj.topology, ligand_heavy,
                                              residue_indices)
    step = max(1, int(np.ceil(traj.n_frames / MAX_FINGERPRINT_FRAMES)))
    sampled = traj[::step]
    if not len(pairs):
        return np.zeros((sampled.n_frames, 0), dtype=bool), [], step

    columns: dict[tuple[int, str], np.ndarray] = {}
    periodic = sampled.unitcell_vectors is not None
    for begin in range(0, len(pairs), PAIR_BLOCK):
        block = pairs[begin:begin + PAIR_BLOCK]
        distances = md.compute_distances(sampled, block, periodic=periodic) * 10.0
        for offset, (residue, itype) in enumerate(labels[begin:begin + len(block)]):
            close = distances[:, offset] <= cutoffs[itype]
            key = (residue, itype)
            if key in columns:
                np.logical_or(columns[key], close, out=columns[key])
            else:
                columns[key] = close

    # Only interactions that were actually seen. The candidate set above is every
    # pair the topology makes geometrically possible, which for a pocket is
    # mostly pairs that never come within their cutoff, and keeping those as
    # all-zero columns quietly ruins two of the numbers below:
    # `hbond_occupancy_mean` becomes the mean over dozens of zeros rather than
    # over the hydrogen bonds that exist, and the occupancy table fills with rows
    # claiming a contact at 0%. A ProLIF fingerprint has no such columns, so
    # dropping them is also what keeps these numbers comparable to one.
    ordered = sorted((key for key, column in columns.items() if column.any()),
                     key=lambda key: (key[0], key[1]))
    bits = np.stack([columns[key] for key in ordered], axis=1) if ordered \
        else np.zeros((sampled.n_frames, 0), dtype=bool)
    return bits, ordered, step


def _mean_lifetime(column: np.ndarray, frame_ns: float) -> float | None:
    """The mean length of an uninterrupted run of presence, in nanoseconds.

    This is what separates one long contact from a hundred flickers at the same
    occupancy, which is why the spec asks for it rather than only the fraction.
    """
    if not column.any() or frame_ns <= 0:
        return None
    padded = np.concatenate([[0], column.astype(np.int8), [0]])
    edges = np.diff(padded)
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return float(np.mean(ends - starts) * frame_ns)


def interaction_records(traj, cid: str, ligand_heavy: np.ndarray, residue_indices,
                        start: int, time_ns: np.ndarray
                        ) -> tuple[dict | None, list[dict] | None, dict | None]:
    """The interaction record, the occupancy table and the barcode.

    `S_int` is the mean of the key-interaction occupancy and the Tanimoto to
    frame 0, so both have to come from one fingerprint: computing them from two
    different contact definitions would average two different claims.
    """
    bits, columns, step = fingerprint(traj, ligand_heavy, residue_indices)
    if not len(columns):
        return None, None, None

    topology = traj.topology
    sampled_time = time_ns[::step][: bits.shape[0]]
    frame_ns = float(np.median(np.diff(sampled_time))) if len(sampled_time) > 1 else 0.0

    # The scored window, in the sampled frame index rather than the full one.
    window_start = min(int(np.ceil(start / step)), max(0, bits.shape[0] - 1))
    window = bits[window_start:]
    if not len(window):
        return None, None, None

    at_t0 = bits[0]
    occupancy = window.mean(axis=0)
    notes: list[str] = [
        "interactions are typed from heavy-atom geometry (polar contact, ionic "
        "and hydrophobic) rather than from a chemistry toolkit, so nothing here "
        "is claimed as a pi-stack and the type vocabulary is coarser than a "
        "ProLIF fingerprint's"
    ]
    if step > 1:
        notes.append(f"the fingerprint is sampled every {step} frames")

    record: dict[str, Any] = {
        "complex_id": cid,
        # DeepFlex's scorer skips the interaction sub-score unless the record
        # carries this, and the spec does not mention it.
        "status": "ok",
        "n_t0_interactions": int(at_t0.sum()),
    }

    if at_t0.any():
        retained = window[:, at_t0].sum(axis=1) / at_t0.sum()
        record["pose_memory_mean"] = float(retained.mean())
        intersection = (window & at_t0).sum(axis=1).astype(float)
        union = (window | at_t0).sum(axis=1).astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            tanimoto = np.where(union > 0, intersection / union, np.nan)
        if np.isfinite(tanimoto).any():
            record["tanimoto_t0_mean"] = float(np.nanmean(tanimoto))
        # Present late and absent at frame 0.
        record["n_interactions_gained"] = int((window.any(axis=0) & ~at_t0).sum())
    else:
        notes.append(
            "no typed interaction at frame 0, so the pose memory and the Tanimoto "
            "to frame 0 are undefined and omitted rather than written as zero")

    key = occupancy >= KEY_INTERACTION_MIN_OCCUPANCY
    record["n_key_interactions"] = int(key.sum())
    if key.any():
        record["key_interaction_occupancy_mean"] = float(occupancy[key].mean())

    hbond = np.array([itype == "HBond" for _, itype in columns], dtype=bool)
    if hbond.any():
        record["hbond_occupancy_mean"] = float(occupancy[hbond].mean())

    if record["n_t0_interactions"] < 4:
        notes.append(
            f"only {record['n_t0_interactions']} typed interaction(s) at frame 0: "
            f"a fingerprint that coarse is not trusted on its own and DeepFlex "
            f"flags the complex")
    record["notes"] = notes

    # ---- the occupancy table, section 5.5 ----------------------------------
    summary: list[dict[str, Any]] = []
    for index, (residue_index, itype) in enumerate(columns):
        # A contact seen only in the discarded head, and not at frame 0, says
        # nothing about the scored window: it would appear as a row at 0%
        # occupancy, which reads as a contact that never formed. It is still in
        # the barcode, which covers the whole trajectory, and that is where a
        # reader should see it.
        if occupancy[index] <= 0.0 and not at_t0[index]:
            continue
        residue = topology.residue(int(residue_index))
        summary.append({
            "residue": f"{residue.name}{residue.resSeq}",
            "resid": int(residue.resSeq),
            "itype": itype,
            "occupancy": float(occupancy[index]),
            "mean_lifetime_ns": _mean_lifetime(window[:, index], frame_ns),
            # 1 or 0 from a measurement, never a guess: a 0 claims the contact
            # formed later and a 1 claims it was there from the start.
            "present_at_t0": int(bool(at_t0[index])),
        })

    # ---- the barcode, section 5.3 ------------------------------------------
    # Busiest first, because the tail is contacts present in a handful of frames
    # and twenty-eight rows is already a tall panel.
    full_occupancy = bits.mean(axis=0)
    order = np.argsort(-full_occupancy)[:fxa2.MAX_BARCODE_ROWS]
    order = order[full_occupancy[order] > 0]
    column_step = max(1, int(np.ceil(bits.shape[0] / fxa2.MAX_BARCODE_COLS)))
    shown = bits[::column_step][:, order]
    labels = []
    for index in order:
        residue_index, itype = columns[int(index)]
        residue = topology.residue(int(residue_index))
        labels.append(f"{residue.name}{residue.resSeq} {itype}")

    barcode = {
        "z": shown.T.astype(np.uint8),
        "labels": np.array(labels, dtype=object),
        "time_ns": sampled_time[::column_step].astype(np.float64),
        "n_contacts_total": np.int32(len(columns)),
        # Frames of the original trajectory per barcode column.
        "step": np.int32(step * column_step),
    }
    return record, summary, barcode


# ===========================================================================
#  Convergence and pose clusters, sections 4 and 6.4
# ===========================================================================


def jensen_shannon(a: np.ndarray, b: np.ndarray, edges: np.ndarray) -> float:
    """JS divergence between two samples on a shared bin grid. 0 is identical.

    Base-2 logarithms, so the divergence is bounded at 1 and the spec's
    thresholds (full marks at 0.05, zero at 0.30) are on that scale. The
    *divergence*, not the distance: the distance is its square root and would
    read about 0.24 where this reads 0.06, which is the difference between full
    marks and most of the way to none.
    """
    pa = np.histogram(a, bins=edges)[0].astype(float)
    pb = np.histogram(b, bins=edges)[0].astype(float)
    if pa.sum() == 0 or pb.sum() == 0:
        return float("nan")
    pa /= pa.sum()
    pb /= pb.sum()
    mean = 0.5 * (pa + pb)

    def kl(p: np.ndarray, q: np.ndarray) -> float:
        mask = p > 0
        return float(np.sum(p[mask] * np.log2(p[mask] / q[mask])))

    return 0.5 * kl(pa, mean) + 0.5 * kl(pb, mean)


def block_divergence(window: np.ndarray) -> tuple[float | None, str | None]:
    """Mean pairwise divergence between four equal blocks of the scored window.

    With one trajectory per complex this stands in for replica agreement and is
    weaker evidence than one: four blocks of the same run share their starting
    conformation, so agreement between them says the run is stationary rather
    than that it found the same answer twice. `n_strides: 1` in the manifest is
    what tells DeepFlex to say so on the page.

    The bin count follows the window length rather than sitting at twenty, and
    that matters more than it looks. A histogram divergence between two finite
    samples is biased upward by roughly the bin count over the sample size, so at
    a fixed twenty bins a stationary 400-frame window reads 0.086 against a
    threshold whose full marks are 0.05: a perfectly converged short run would be
    scored as a diverging one.

    Holding samples-per-bin constant instead keeps the bias constant, and the
    divisor is chosen so that a long run lands on twenty bins, which is what
    DeepFlex uses and therefore what the spec's thresholds (full marks at 0.05,
    zero at 0.30) were calibrated against. Measured on stationary noise, that
    keeps the divergence under 0.02 from 60 frames to 1,800.

    Below about 40 frames the bias is unavoidable and the note says so: the
    number then reads worse than a longer run of the same trajectory would, which
    penalises rather than flatters a short run.
    """
    blocks = np.array_split(window, CONV_BLOCKS)
    shortest = min(len(block) for block in blocks)
    if shortest < MIN_FRAMES_PER_BLOCK:
        return None, (
            f"the scored window is {len(window)} frame(s), which leaves "
            f"{shortest} per block against the {MIN_FRAMES_PER_BLOCK} needed, so "
            f"the convergence sub-score is left pending rather than computed on "
            f"histograms too sparse to compare")

    bins = int(np.clip(round(len(window) / CONV_FRAMES_PER_BIN), 3, CONV_BINS))
    # One bin grid for all four, from the whole window. Per-block grids would
    # compare each distribution against a different ruler.
    edges = np.histogram_bin_edges(window, bins=bins)
    divergences = [
        jensen_shannon(blocks[i], blocks[j], edges)
        for i in range(CONV_BLOCKS) for j in range(i + 1, CONV_BLOCKS)
    ]
    divergences = [d for d in divergences if np.isfinite(d)]
    if not divergences:
        return None, "the four blocks were empty, so the divergence is undefined"

    # Only where the bias is actually material. Holding samples-per-bin constant
    # keeps it negligible down to about sixty frames; below that the bin count is
    # already at its floor and the number genuinely reads worse than a longer run
    # of the same trajectory would.
    note = None
    per_bin = (len(window) / CONV_BLOCKS) / bins
    if per_bin < 10:
        note = (
            f"the block divergence is over {bins} bins with about {per_bin:.0f} "
            f"frame(s) per bin per block, because the scored window is only "
            f"{len(window)} frames. The finite-sample bias inflates it at that "
            f"sampling, so it reads worse than a longer run of the same "
            f"trajectory would")
    return float(np.mean(divergences)), note


def pose_clusters(ligand_xyz: np.ndarray, pose_rmsd: np.ndarray,
                  start: int) -> dict[str, Any] | None:
    """Cluster the scored window's ligand poses, and say whether the pose changed.

    Average-linkage hierarchical clustering on pairwise ligand heavy-atom RMSD,
    cut at 2 A. Returns the cluster records, the dominant population and the
    changed-pose flag, or None when there is too little of the window to cluster.
    """
    from scipy.cluster.hierarchy import fcluster, linkage  # noqa: PLC0415
    from scipy.spatial.distance import squareform  # noqa: PLC0415

    n_frames = len(pose_rmsd)
    scored = np.arange(n_frames)[start:]
    stride = max(1, int(np.ceil(len(scored) / MAX_CLUSTER_FRAMES)))
    frames = scored[::stride]
    if len(frames) < 4:
        return None

    poses = ligand_xyz[frames]
    difference = poses[:, None, :, :] - poses[None, :, :, :]
    matrix = np.sqrt((difference ** 2).sum(axis=3).mean(axis=2))
    np.fill_diagonal(matrix, 0.0)

    labels = fcluster(linkage(squareform(matrix, checks=False), method="average"),
                      t=POSE_CLUSTER_CUTOFF_A, criterion="distance")
    counts = np.bincount(labels)[1:]
    populations = counts / counts.sum()

    records: list[dict[str, Any]] = []
    for label in range(1, len(counts) + 1):
        members = np.flatnonzero(labels == label)
        if not len(members):
            continue
        # The centroid is the member with the smallest mean distance to the
        # others, which is a real frame and can therefore be written out as a
        # structure. An averaged set of coordinates could not.
        within = matrix[np.ix_(members, members)]
        centroid = int(frames[members[int(within.mean(axis=1).argmin())]])
        records.append({
            "cluster_id": int(label),
            "population": float(counts[label - 1] / counts.sum()),
            # An index into the FULL trajectory, not into the strided subsample,
            # so the frame it names can be loaded directly.
            "centroid_frame": centroid,
            "rmsd_to_t0_a": float(pose_rmsd[centroid]),
            "n_frames": int(counts[label - 1]) * stride,
        })
    records.sort(key=lambda entry: -entry["population"])

    # The dominant cluster of each half, rather than the first and last frame.
    # Comparing single frames fires on ordinary thermal motion; section 4 means a
    # ligand that moves to a new pose and then HOLDS it, so the second half's
    # winner has to actually dominate its half.
    half = len(labels) // 2
    changed = False
    if half >= 2:
        first = int(np.bincount(labels[:half]).argmax())
        tail = np.bincount(labels[half:])
        second = int(tail.argmax())
        changed = bool(first != second and tail[second] / tail.sum() >= 0.5)

    return {
        "n_pose_clusters": int(len(counts)),
        "dominant_cluster_population": float(populations.max()),
        "changed_pose_cluster": int(changed),
        "pose_clusters": records,
        "centroid_frames": [record["centroid_frame"] for record in records],
    }


# ===========================================================================
#  Ligand strain and torsions, sections 4 and 5.4
# ===========================================================================


def ligand_mol(traj, atoms: np.ndarray):
    """An RDKit molecule for the ligand, with bond orders, or None.

    The bond orders are **inferred from the trajectory's own coordinates**
    rather than read from the chemical definition the bundle fetched. Three
    reasons, in order of weight:

    * It works for both source kinds. A version 1 `.fxa` carries no SDF at all,
      and a strain that only appears for one kind of source is a column that
      silently means different things.
    * No chemistry file is opened at any point, so there is no path by which one
      could end up in the package. Section 8 is absolute about that and the
      cheapest way to honour a rule is to have no code that could break it.
    * Matching an SDF's atom order onto an MD topology's is its own failure mode,
      and a mismatch produces a plausible wrong number rather than an error.

    It needs hydrogens: the valences are what make the inference possible, and
    without them MMFF cannot type the molecule either. A payload that dropped
    them returns None and the sub-score stays pending, which is honest.
    """
    from rdkit import Chem, RDLogger  # noqa: PLC0415
    from rdkit.Chem import rdDetermineBonds  # noqa: PLC0415

    RDLogger.DisableLog("rdApp.*")

    hydrogens = sum(1 for i in atoms
                    if traj.topology.atom(int(i)).element.symbol == "H")
    if hydrogens == 0:
        return None

    block = _pdb_text(traj[0].atom_slice(atoms))
    base = Chem.MolFromPDBBlock(block, removeHs=False, proximityBonding=True,
                                sanitize=False)
    if base is None or base.GetNumAtoms() != len(atoms):
        return None

    for charge in CANDIDATE_CHARGES:
        candidate = Chem.Mol(base)
        try:
            rdDetermineBonds.DetermineBondOrders(candidate, charge=charge)
            Chem.SanitizeMol(candidate)
        except Exception:  # noqa: BLE001 - a failed assignment is tried again
            continue
        return candidate
    return None


def rotatable_torsions(mol) -> list[tuple[int, int, int, int]]:
    """One dihedral quadruple per rotatable bond.

    The heaviest neighbour is taken on each side rather than the first, so the
    dihedral tracks the substituent that actually moves rather than a hydrogen
    that spins freely and would register a flip on every frame.
    """
    from rdkit import Chem  # noqa: PLC0415

    pattern = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
    torsions: list[tuple[int, int, int, int]] = []
    for j, k in mol.GetSubstructMatches(pattern):
        def heaviest(centre: int, other: int) -> int | None:
            best, best_mass = None, -1.0
            for neighbour in mol.GetAtomWithIdx(centre).GetNeighbors():
                index = neighbour.GetIdx()
                if index == other:
                    continue
                if neighbour.GetMass() > best_mass:
                    best, best_mass = index, neighbour.GetMass()
            return best

        first, last = heaviest(j, k), heaviest(k, j)
        if first is not None and last is not None:
            torsions.append((first, j, k, last))
    return torsions


def assign_wells(angles: np.ndarray) -> np.ndarray:
    """Rotamer well per frame, with hysteresis so boundary rattle is not a flip."""
    low, high = WELL_EDGES
    wells = np.zeros(angles.size, dtype=np.int8)
    # Frame 0 with no hysteresis: there is nothing to be sticky about yet.
    current = 0 if angles[0] < low else (1 if angles[0] < high else 2)
    wells[0] = current
    for index in range(1, angles.size):
        angle = angles[index]
        candidate = 0 if angle < low else (1 if angle < high else 2)
        if candidate != current:
            boundary = low if {current, candidate} == {0, 1} else high
            if abs(angle - boundary) >= HYSTERESIS_DEG:
                current = candidate
        wells[index] = current
    return wells


def strain_records(mol, traj, atoms: np.ndarray, start: int, time_ns: np.ndarray
                   ) -> tuple[dict[str, Any], dict[str, np.ndarray], list[str]]:
    """MMFF94 torsional strain in excess of frame 0, and the torsion series.

    **The excess, not the absolute energy.** These coordinates come from an
    OpenFF or GAFF force field and are being scored by MMFF94, and a 298 K
    snapshot of a drug-sized ligand carries tens of kcal/mol in its internal
    modes. Both offsets are constant for a given ligand and cancel when frame 0
    is subtracted; an absolute MMFF number is not comparable between ligands and
    would be scored as if it were.

    Per frame, the strain is the cost of *holding the observed torsions*:

    1. minimise with every rotatable torsion restrained to its observed value,
       which relaxes the bonds and angles the MD constrained differently;
    2. minimise again with nothing restrained, which falls into the nearest free
       minimum because the conformer has already relaxed;
    3. take the difference.

    Done that way the bond and angle terms largely cancel inside each frame as
    well as between frames, so what is left is torsional.
    """
    from rdkit.Chem import AllChem, rdMolTransforms  # noqa: PLC0415
    from rdkit.Geometry import Point3D  # noqa: PLC0415

    if mol is None:
        return {}, {}, [
            "the ligand's bond orders could not be inferred from the "
            "trajectory (most often because the payload carries no hydrogens), "
            "so the strain and torsion measurements are omitted and the "
            "conformation sub-score stays pending"
        ]

    properties = AllChem.MMFFGetMoleculeProperties(mol)
    if properties is None:
        return {}, {}, [
            "MMFF94 has no parameters for this ligand, so the strain and torsion "
            "measurements are omitted and the conformation sub-score stays pending"
        ]

    torsions = rotatable_torsions(mol)
    conformer = mol.GetConformer()

    # Frame 0 prepended, so the reference is measured by exactly the same code
    # path as every scored frame. A frame-0 value computed differently would
    # reintroduce the offset that subtracting it exists to cancel.
    scored = np.arange(traj.n_frames)[start:]
    stride = max(1, int(np.ceil(len(scored) / MAX_STRAIN_FRAMES)))
    frames = np.concatenate([[0], scored[::stride]])
    coordinates = traj.xyz[frames] * 10.0          # nanometres to angstrom

    strain = np.zeros(len(frames), dtype=np.float32)
    angles = np.zeros((len(frames), len(torsions)), dtype=np.float32)

    for row in range(len(frames)):
        for atom, position in enumerate(coordinates[row][atoms]):
            conformer.SetAtomPosition(
                atom, Point3D(float(position[0]), float(position[1]),
                              float(position[2])))

        # Read the dihedrals before minimising anything: the minimisation moves
        # the conformer, so a reading taken afterwards is of the relaxed pose.
        for index, (i, j, k, l) in enumerate(torsions):
            angles[row, index] = rdMolTransforms.GetDihedralDeg(conformer, i, j, k, l)

        held = AllChem.MMFFGetMoleculeForceField(mol, properties)
        free = AllChem.MMFFGetMoleculeForceField(mol, properties)
        if held is None or free is None:
            continue
        for i, j, k, l in torsions:
            current = rdMolTransforms.GetDihedralDeg(conformer, i, j, k, l)
            held.MMFFAddTorsionConstraint(i, j, k, l, False, current, current,
                                         TORSION_RESTRAINT_K)
        held.Minimize(maxIts=MAX_STRAIN_ITERATIONS)
        held_energy = held.CalcEnergy()
        free.Minimize(maxIts=MAX_STRAIN_ITERATIONS)
        strain[row] = float(max(0.0, held_energy - free.CalcEnergy()))

    # Row 0 is frame 0 and is not part of the scored window.
    strain_t0 = float(strain[0])
    angles_t0 = angles[0].copy()
    strain, angles = strain[1:], angles[1:]
    frames = frames[1:]
    if not len(strain):
        return {}, {}, ["the scored window held no frames to measure strain over"]

    flips = [int(np.count_nonzero(np.diff(assign_wells(angles[:, index]))))
             for index in range(len(torsions))]
    span_ns = float(time_ns[frames[-1]] - time_ns[frames[0]]) if len(frames) > 1 \
        else 0.0

    metrics: dict[str, Any] = {
        # The excess over frame 0, which is the number that gets scored.
        "strain_kcal": float(np.median(strain) - strain_t0),
        # Before the subtraction, for audit, and the value subtracted.
        "strain_absolute_kcal": float(np.median(strain)),
        "strain_t0_kcal": strain_t0,
        "n_rotatable_bonds": len(torsions),
        "n_torsion_flips": int(sum(flips)),
    }
    if span_ns > 0 and torsions:
        metrics["flips_per_ns_per_torsion"] = float(
            sum(flips) / (span_ns * len(torsions)))
    elif torsions:
        metrics["flips_per_ns_per_torsion"] = 0.0

    arrays = {
        "frame": frames.astype(np.int32),
        "time_ns": time_ns[frames].astype(np.float32),
        "strain_kcal": strain.astype(np.float32),
        "strain_t0_kcal": np.float32(strain_t0),
        "strain_excess_kcal": (strain - strain_t0).astype(np.float32),
        "torsion_deg": angles.astype(np.float32),
        "torsion_deg_t0": angles_t0.astype(np.float32),
        "torsion_atoms": np.array(torsions, dtype=np.int32).reshape(-1, 4),
        "flips_per_torsion": np.array(flips, dtype=np.int32),
    }

    notes: list[str] = [
        "the ligand's bond orders were inferred from the trajectory's own "
        "coordinates and hydrogens rather than read from a chemical definition, "
        "so the MMFF94 atom typing rests on that inference"
    ]
    if not torsions:
        notes.append(
            "the ligand has no rotatable bond, so the torsion-flip term is "
            "vacuous and the conformation sub-score rests on the strain alone")
    return metrics, arrays, notes


# ===========================================================================
#  Receptor panels, section 5.2
# ===========================================================================


def _summarise(values: np.ndarray | None) -> dict[str, float | None]:
    """The five numbers the layer manifest carries for a series."""
    if values is None or values.size == 0:
        return {"mean": None, "sd": None, "min": None, "max": None, "last": None}
    return {
        "mean": float(np.mean(values)), "sd": float(np.std(values)),
        "min": float(np.min(values)), "max": float(np.max(values)),
        "last": float(values[-1]),
    }


def receptor_panels(traj, fitted, receptor_ca: np.ndarray, start: int,
                    time_ns: np.ndarray, rmsd_ca: np.ndarray
                    ) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """The receptor's own flexibility: series, per-residue RMSF and two matrices.

    `fitted` is the trajectory already superposed on the receptor Calpha, so the
    RMSF and the correlations describe internal motion rather than the whole
    molecule tumbling.
    """
    import mdtraj as md  # noqa: PLC0415

    topology = traj.topology
    protein = topology.select("protein")

    arrays: dict[str, np.ndarray] = {
        "time_ns": time_ns.astype(np.float32),
        "rmsd_ca_a": rmsd_ca.astype(np.float32),
    }
    summary: dict[str, Any] = {"rmsd_ca": _summarise(rmsd_ca[start:])}

    if len(protein):
        radgyr = md.compute_rg(traj.atom_slice(protein)) * 10.0
        arrays["radgyr_a"] = radgyr.astype(np.float32)
        summary["radgyr"] = _summarise(radgyr[start:])

    # Per-residue Calpha fluctuation about the mean position over the scored
    # window. `rmsf_resid` and `rmsf_a` must be the same length and the same
    # order, and the report marks pocket residues on this plot by looking their
    # numbers up in pocket/<cid>.json, so the numbering has to agree with it:
    # resSeq in both, ascending.
    if len(receptor_ca) >= 3:
        coords = fitted.xyz[start:][:, receptor_ca, :] * 10.0
        mean_position = coords.mean(axis=0)
        fluctuation = np.sqrt(((coords - mean_position) ** 2).sum(axis=2).mean(axis=0))
        resids = np.array([topology.atom(int(i)).residue.resSeq for i in receptor_ca],
                          dtype=np.int32)
        order = np.argsort(resids, kind="stable")
        arrays["rmsf_resid"] = resids[order]
        arrays["rmsf_a"] = fluctuation[order].astype(np.float32)
        summary["rmsf"] = _summarise(fluctuation)
        summary["n_rmsf_residues"] = int(len(resids))

        # Frame against frame Calpha RMSD, reduced before it is written.
        stride = max(1, int(np.ceil(traj.n_frames / MAX_PANEL_SIDE)))
        sample = fitted[::stride]
        matrix = np.empty((sample.n_frames, sample.n_frames), dtype=np.float32)
        for frame in range(sample.n_frames):
            matrix[frame] = md.rmsd(sample, sample, frame,
                                    atom_indices=receptor_ca) * 10.0
        arrays["rmsd_2d"] = matrix

        # Residue cross-correlation, from the same fitted Calpha displacements.
        # Residues are decimated rather than the matrix being cropped, so the
        # panel still spans the whole chain when there are more than 500 of them.
        block = max(1, int(np.ceil(len(receptor_ca) / MAX_PANEL_SIDE)))
        selected = np.arange(len(receptor_ca))[::block]
        displacement = (coords[:, selected, :] - mean_position[selected]).reshape(
            coords.shape[0], -1, 3)
        dot = np.einsum("fid,fjd->ij", displacement, displacement) / coords.shape[0]
        variance = np.sqrt(np.diag(dot))
        with np.errstate(divide="ignore", invalid="ignore"):
            correlation = dot / np.outer(variance, variance)
        arrays["cross_correlation"] = np.nan_to_num(
            np.clip(correlation, -1.0, 1.0), nan=0.0).astype(np.float32)
        arrays["dccm_resid"] = resids[order][selected]

        # The first two principal components of the fitted Calpha coordinates,
        # over the whole trajectory so the discarded head is visible on the plot
        # with everything else.
        whole = fitted.xyz[:, receptor_ca, :].reshape(traj.n_frames, -1) * 10.0
        centred = whole - whole.mean(axis=0)
        if traj.n_frames > 2:
            # SVD rather than an eigendecomposition of the covariance: the
            # covariance of 3N coordinates is 3N by 3N, which for a 200-residue
            # receptor is a 600 by 600 matrix built to extract two vectors.
            _, _, components = np.linalg.svd(centred, full_matrices=False)
            projection = centred @ components[:2].T
            arrays["pc1_package"] = projection[:, 0].astype(np.float32)
            arrays["pc2_package"] = projection[:, 1].astype(np.float32)

    return arrays, summary


# ===========================================================================
#  The 2D depiction, section 6.7
# ===========================================================================

# Low to high mobility. Two hues rather than a rainbow, because the quantity is
# one-dimensional and a reader has to be able to tell which end is which without
# a key; and because a red-green ramp is the one choice that fails for a common
# form of colour blindness.
MOBILITY_LOW = (0.18, 0.44, 0.71)       # blue
MOBILITY_HIGH = (0.84, 0.37, 0.09)      # orange
DEPICTION_SIZE = (340, 260)


def depiction(mol, atoms: np.ndarray, heavy: np.ndarray,
              ligand_rmsf: np.ndarray) -> tuple[str | None, bool]:
    """A 2D diagram of the ligand, shaded by per-atom mobility where possible.

    Returns the SVG and whether it carries the shading. The second value matters:
    the mobility legend is a claim about colours, and `depictions/painted.json`
    lists only the diagrams that actually have them so the legend is not drawn
    beside a plain diagram.
    """
    if mol is None:
        return None, False

    from rdkit import Chem  # noqa: PLC0415
    from rdkit.Chem import AllChem  # noqa: PLC0415
    from rdkit.Chem.Draw import rdMolDraw2D  # noqa: PLC0415

    # The fluctuation is per heavy atom and the molecule carries its hydrogens,
    # so it is attached as an atom property and travels through RemoveHs with
    # the atom it belongs to. Mapping by position afterwards would be mapping
    # onto a renumbered molecule.
    position_of = {int(index): order for order, index in enumerate(atoms)}
    painted = False
    if len(ligand_rmsf) == len(heavy):
        low, high = float(ligand_rmsf.min()), float(ligand_rmsf.max())
        span = high - low
        for value, index in zip(ligand_rmsf, heavy):
            order = position_of.get(int(index))
            if order is None or order >= mol.GetNumAtoms():
                continue
            scaled = 0.5 if span <= 0 else (float(value) - low) / span
            mol.GetAtomWithIdx(order).SetDoubleProp("mobility", scaled)
        painted = True

    try:
        flat = Chem.RemoveHs(Chem.Mol(mol))
        AllChem.Compute2DCoords(flat)
    except Exception:  # noqa: BLE001 - a diagram is polish, not a measurement
        return None, False

    highlights: list[int] = []
    colours: dict[int, tuple[float, float, float]] = {}
    for atom in flat.GetAtoms():
        if not atom.HasProp("mobility"):
            continue
        scaled = atom.GetDoubleProp("mobility")
        highlights.append(atom.GetIdx())
        colours[atom.GetIdx()] = tuple(
            low_channel + scaled * (high_channel - low_channel)
            for low_channel, high_channel in zip(MOBILITY_LOW, MOBILITY_HIGH))
    if not highlights:
        painted = False

    drawer = rdMolDraw2D.MolDraw2DSVG(*DEPICTION_SIZE)
    drawer.drawOptions().addStereoAnnotation = False
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, flat, highlightAtoms=highlights or None,
        highlightAtomColors=colours or None)
    drawer.FinishDrawing()
    return drawer.GetDrawingText(), painted


# ===========================================================================
#  The whole measurement for one complex
# ===========================================================================


def measure_complex(traj, *, complex_id: str, receptor: str, ligand_id: str,
                    ligand_resnames: list[str],
                    frame_spacing_ps: float | None = None) -> Measured:
    """Measure one receptor-ligand complex.

    Raises nothing: a complex that cannot be measured comes back as a failure
    record, because the caller is assembling a package and one bad trajectory
    must not cost the other nineteen.
    """
    try:
        return _measure(traj, complex_id, receptor, ligand_id, ligand_resnames,
                        frame_spacing_ps)
    except MeasureError as exc:
        return failed(complex_id, receptor, ligand_id, str(exc))


def _measure(traj, cid: str, receptor: str, ligand: str, resnames: list[str],
             frame_spacing_ps: float | None) -> Measured:
    atoms, heavy, fragmented = ligand_atoms(traj.topology, resnames)
    pocket = pocket_residues(traj, heavy)

    topology = traj.topology
    ligand_residue = topology.atom(int(heavy[0])).residue.index

    # The three atom sets the measurements fit on. The pocket Calpha set is the
    # one that matters most: it is the superposition frame for the ligand RMSD,
    # and it is the filtered pocket rather than the crop shell for that reason.
    pocket_ca = _atoms_of(topology, pocket["indices_filtered"], name="CA")
    pocket_heavy = _atoms_of(topology, pocket["indices_filtered"], heavy=True)
    receptor_ca = topology.select("protein and name CA")
    if len(pocket_ca) < 3:
        raise MeasureError("pocket_ca_too_few")

    # The crop: the ligand plus every whole residue in the 10 A shell, with the
    # solvent removed. Selected by residue index rather than by residue number,
    # because residue numbers are not unique once insertion codes have been
    # dropped and selecting on one pulls in whichever other residue shares it.
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

    time_ns, time_note = time_axis(traj, frame_spacing_ps)
    if time_note:
        notes.append(time_note)

    # ---- the periodic-boundary gate, before anything is superposed ----------
    #
    # A ligand wrapped across a boundary sits a box length from its receptor in
    # the stored coordinates, and then every distance that spans the two is wrong
    # by that much. `pbc_ok: false` is not a failure, but it makes every distance
    # in the file suspect and the report says so, which is why it is measured
    # rather than defaulted to true.
    raw_com = traj.xyz[:, heavy, :].mean(axis=1) * 10.0
    jumps = np.linalg.norm(np.diff(raw_com, axis=0), axis=1) if traj.n_frames > 1 \
        else np.zeros(0)
    n_jumps = int(np.count_nonzero(jumps > PBC_JUMP_A))
    pbc_ok = n_jumps == 0
    if not pbc_ok:
        notes.append(
            f"{n_jumps} ligand centroid jump(s) over {PBC_JUMP_A:.0f} A between "
            f"consecutive frames (largest {jumps.max():.1f} A): the coordinates "
            f"are wrapped across a periodic boundary, so every distance here is "
            f"suspect"
        )

    # ---- contacts, also before any superposition ---------------------------
    #
    # This has to precede the fits. `superpose` rotates the coordinates and
    # leaves the unitcell vectors untouched, so a periodic distance taken
    # afterwards applies stale box vectors to rotated coordinates.
    contacts = residue_contacts(traj, ligand_residue, pocket["indices_raw"])
    n_native = int(contacts[0])
    threshold = max(1, int(np.ceil(n_native * NATIVE_CONTACT_FRACTION)))
    in_pocket = contacts >= threshold

    # ---- the three superpositions ------------------------------------------
    pose_fit = _fit(traj, pocket_ca)
    pose_rmsd = _drift_rmsd(pose_fit, heavy)
    inplace_rmsd = _drift_rmsd(_fit(traj, heavy), heavy)
    receptor_fit = _fit(traj, receptor_ca) if len(receptor_ca) >= 3 else pose_fit
    receptor_rmsd = _drift_rmsd(receptor_fit, receptor_ca) \
        if len(receptor_ca) >= 3 else np.zeros(traj.n_frames)

    # Centroid displacement, measured from the ligand's OWN frame-0 centroid and
    # not from the pocket's. The pocket centroid is the average of the residues
    # lining the site and sits several angstrom from the ligand centre even in a
    # perfectly held pose, so measuring against it adds a constant offset and
    # the 5 A hard gate then fires on ligands that have not moved. The offset is
    # still reported, as a diagnostic rather than as a drift measure.
    ligand_xyz = pose_fit.xyz[:, heavy, :] * 10.0
    com = ligand_xyz.mean(axis=1)
    com_drift = np.linalg.norm(com - com[0], axis=1)
    pocket_centre = (pose_fit.xyz[0, pocket_heavy, :] * 10.0).mean(axis=0) \
        if len(pocket_heavy) else com[0]
    pocket_offset = np.linalg.norm(com - pocket_centre, axis=1)

    # ---- the scored window -------------------------------------------------
    start, method = detect_equilibration(pose_rmsd)
    scored = slice(start, None)
    window = pose_rmsd[scored]
    window_time = time_ns[scored]

    # Per-atom ligand fluctuation about its mean position over the window, after
    # the pocket fit, which is the ligand flexing rather than the ligand moving.
    mean_position = ligand_xyz[scored].mean(axis=0)
    ligand_rmsf = np.sqrt(
        ((ligand_xyz[scored] - mean_position) ** 2).sum(axis=2).mean(axis=0))

    metrics: dict[str, Any] = {
        "complex_id": cid,
        "receptor": receptor,
        "ligand_id": ligand,
        "status": "ok",
        "n_frames": int(traj.n_frames),
        # DeepFlex's scorer reads this directly rather than with a default, and a
        # FlexAppeal study has no reference compound in it.
        "is_reference": False,

        # The scored window. Every quantity below it is computed from this frame
        # onward: a strain or an energy measured over a different window from the
        # RMSD is not the quantity that belongs beside that grade.
        "equil_start_frame": int(start),
        "equil_method": method,

        # S_pose
        "ligand_rmsd_median": float(np.median(window)),
        "ligand_rmsd_p90": float(np.percentile(window, 90)),
        "ligand_rmsd_max": float(window.max()),
        "inplace_rmsd_median": float(np.median(inplace_rmsd[scored])),

        # S_occ
        "n_native_contacts": n_native,
        "in_pocket_fraction": float(in_pocket[scored].mean()),
        # A transition from inside the pocket to outside it, counted over the
        # scored window and not the whole run. Section 4 says "within the
        # window", and every other scored quantity here is cut the same way: an
        # exit in the equilibration head is a frame the file says was discarded,
        # so counting it would put a departure on the page that nothing scored.
        "n_exit_events": int(np.count_nonzero(
            np.diff(in_pocket[scored].astype(np.int8)) == -1)),
        "com_drift_mean_a": float(com_drift[scored].mean()),
        "com_drift_max_a": float(com_drift[scored].max()),
        "pocket_offset_mean_a": float(pocket_offset[scored].mean()),

        # Quality flags
        "pbc_ok": bool(pbc_ok),
        "extra_ligand_fragments": bool(fragmented),
        "receptor_rmsd_mean": float(receptor_rmsd[scored].mean()),
        "ligand_rmsf_mean": float(ligand_rmsf.mean()),
        "ligand_rmsf_max": float(ligand_rmsf.max()),
    }

    slope = _slope_per_ns(window, window_time)
    if slope is not None:
        metrics["ligand_rmsd_drift_slope"] = slope

    if traj.unitcell_volumes is not None:
        # nm^3 to A^3.
        metrics["box_volume_a3"] = float(np.mean(traj.unitcell_volumes) * 1000.0)

    if n_native < 4:
        notes.append(
            f"only {n_native} receptor residue(s) in contact with the ligand at "
            f"frame 0, so the occupancy threshold is {threshold} residue(s) and "
            f"in_pocket_fraction is a coarse measure here"
        )

    # ---- convergence and the pose clusters ---------------------------------
    divergence, divergence_note = block_divergence(window)
    if divergence is not None:
        metrics["conv_js_mean"] = divergence
    if divergence_note:
        notes.append(divergence_note)

    clusters = pose_clusters(ligand_xyz, pose_rmsd, start)
    if clusters is not None:
        clusters.pop("centroid_frames")
        metrics.update(clusters)
        if clusters["changed_pose_cluster"]:
            notes.append(
                "the dominant pose cluster of the second half of the scored window "
                "differs from the first half's and the new pose holds: DeepFlex "
                "caps the grade at B and flags it as an alternative pose")
    else:
        notes.append(
            "too little of the scored window to cluster the poses, so the "
            "convergence sub-score is left pending")

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
        # The trajectory as one multi-model PDB, cropped to the pocket and
        # superposed on it, so the receptor holds still and what moves is the
        # ligand. Written from the pocket-fitted copy for that reason. Without
        # this member the viewer shows a single frame and the trajectory control
        # is absent.
        structures[f"{cid}_scrub.pdb"] = scrub_pdb(pose_fit, crop_atoms, resnames)

    # One frame per pose-cluster centroid, numbered to match `cluster_id`. The
    # pose-cluster panel links each row to its centroid, and a row without a
    # file is listed without a link rather than dropped.
    for record in metrics.get("pose_clusters") or []:
        frame = int(record["centroid_frame"])
        structures[f"{cid}_cluster{record['cluster_id']}.pdb"] = pocket_pdb(
            pose_fit, frame, crop_atoms, resnames)

    # ---- the receptor's own panels -----------------------------------------
    panels, panels_summary = receptor_panels(
        traj, receptor_fit, receptor_ca, start, time_ns, receptor_rmsd)
    for key, series_key in (("rmsd_ca", "rmsd_ca"), ("radgyr", "radgyr"),
                            ("rmsf", "rmsf_ca")):
        block = panels_summary.get(key) or {}
        if block.get("mean") is not None:
            metrics[f"{series_key}_mean"] = block["mean"]
            metrics[f"{series_key}_max"] = block["max"]
            if series_key == "radgyr":
                metrics["radgyr_sd"] = block["sd"]
    if panels_summary.get("n_rmsf_residues"):
        metrics["n_rmsf_residues"] = panels_summary["n_rmsf_residues"]

    # ---- ligand strain -----------------------------------------------------
    mol = ligand_mol(traj, atoms)
    strain_metrics, conformer, strain_notes = strain_records(
        mol, traj, atoms, start, time_ns)
    metrics.update(strain_metrics)
    notes.extend(strain_notes)
    # The same four numbers again, under the spelling DeepFlex's layer manifest
    # uses. `n_rotatable` rather than `n_rotatable_bonds` is that manifest's own
    # name for it; the spec uses the column's name, so the mapping is explicit
    # here rather than left to a coincidence.
    conformer_line = {
        "n_rotatable": strain_metrics["n_rotatable_bonds"],
        **{key: value for key, value in strain_metrics.items()
           if key != "n_rotatable_bonds"},
    } if strain_metrics else {}

    svg, svg_painted = depiction(mol, atoms, heavy, ligand_rmsf)

    interactions, summary, barcode = interaction_records(
        traj, cid, heavy, pocket["indices_raw"], start, time_ns)
    if interactions is None:
        notes.append(
            "no typed interaction was found in any frame, so the interaction "
            "sub-score is left pending rather than scored on nothing")

    # The whole trajectory, not the scored window: DeepFlex draws the discarded
    # head in a different colour, which is how a reader sees what was thrown
    # away. It decimates to about 400 points for the wire itself.
    series = {
        "time_ns": time_ns.astype(np.float64),
        "pose_rmsd_a": pose_rmsd.astype(np.float32),
        "inplace_rmsd_a": inplace_rmsd.astype(np.float32),
        "receptor_rmsd_a": receptor_rmsd.astype(np.float32),
        "com_drift_a": com_drift.astype(np.float32),
        "pocket_offset_a": pocket_offset.astype(np.float32),
        "contacts": contacts.astype(np.int16),
        "in_pocket": in_pocket.astype(np.int8),
        # Over the ligand's heavy atoms, which is what the fluctuation was
        # measured on: a hydrogen's RMSF is its parent's plus a libration.
        "ligand_rmsf_a": ligand_rmsf.astype(np.float32),
        "equil_start_frame": np.array([start], dtype=np.int64),
        "native_contact_threshold": np.array([threshold], dtype=np.int64),
    }

    if notes:
        metrics["notes"] = notes

    return Measured(
        complex_id=cid, receptor=receptor, ligand_id=ligand,
        metrics=metrics, pocket=pocket_record, structures=structures,
        series=series, interactions=interactions, summary=summary,
        barcode=barcode, conformer=conformer or None,
        conformer_summary=conformer_line or None,
        panels=panels or None, panels_summary=panels_summary or None,
        depiction_svg=svg, depiction_painted=svg_painted, notes=notes,
    )
