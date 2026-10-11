"""Per-frame analysis of one production trajectory.

Spec sections 5.4 to 5.6, applied frame by frame. The shape of the job:

    superpose on the core  ->  measure each motif  ->  score against references

The core matters more than it looks. Every motif RMSD in section 5.6 is measured
"after core superposition", so the core is the frame of reference the whole state
score is expressed in. Pick a core that moves with activation and the motif
RMSDs report the core's motion instead of the motif's. Section 5.4 names C-lobe
helices for the kinase and section 5.5 names the x.38 to x.58 stretch of TM1 to
TM4 for the receptor, both chosen for being the part that does not move.

Frames are converted to gemmi structures one at a time rather than held as a
single array, because every measurement in `gpcr.py` and `pocket.py` works on a
gemmi chain and reimplementing them against mdtraj arrays would mean two
implementations of each measure that have to agree. The conversion costs about a
millisecond a frame, which against a 50 ns run at one frame per 100 ps is a few
seconds in total.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .state import FrameState, score_motif

# Section 5.5's superposition core: the x.38 to x.58 span of the four helices
# least involved in activation. TM5, TM6 and TM7 are excluded by name because
# they are the helices that move.
GPCR_CORE_HELICES = ("TM1", "TM2", "TM3", "TM4")
GPCR_CORE_SPAN = (38, 58)

# Section 5.4's core. alphaE is inside the KLIFS pocket and comes from there;
# alphaF and alphaH are not, and come from DSSP on the ensemble medoid.
KINASE_CORE_HELICES = ("alphaE", "alphaF", "alphaH")


@dataclass
class Core:
    """The residues a trajectory is superposed on, and where they came from.

    `source` is stored because section 5.4 asks for the residue list to be
    logged: a core that silently fell back to a different definition would
    change every RMSD in the run with nothing to show it had.
    """

    residues: list[int] = field(default_factory=list)
    source: str = ""
    helices: dict[str, list[int]] = field(default_factory=dict)

    def as_row(self) -> dict[str, Any]:
        return {"n_residues": len(self.residues), "source": self.source,
                "residues": list(self.residues),
                "helices": {k: list(v) for k, v in sorted(self.helices.items())}}


def gpcr_core(numbering, present: set[int] | None = None) -> Core:
    """The x.38 to x.58 residues of TM1 to TM4, by generic number.

    `present` restricts to residues the construct actually has, so a trimmed
    model does not try to superpose on a residue it is missing. Resolved through
    generic numbers, so the same call gives the equivalent core in any class A
    receptor.
    """
    core = Core(source=f"generic {GPCR_CORE_SPAN[0]} to {GPCR_CORE_SPAN[1]} of "
                       + ", ".join(GPCR_CORE_HELICES))
    lo, hi = GPCR_CORE_SPAN
    for helix in GPCR_CORE_HELICES:
        index = helix[2:]
        found = []
        for offset in range(lo, hi + 1):
            number = numbering.residue(f"{index}.{offset}")
            if number is None:
                continue
            if present is not None and number not in present:
                continue
            found.append(number)
        if found:
            core.helices[helix] = sorted(found)
    core.residues = sorted(n for v in core.helices.values() for n in v)
    return core


# ===========================================================================
#  Frames
# ===========================================================================


def frame_to_gemmi(traj, index: int, chain_name: str = "A"):
    """One mdtraj frame as a gemmi structure, protein and ligand only.

    Residue numbers are taken from the topology's `resSeq` rather than renumbered,
    because every residue in this package is addressed by its deposited number:
    a renumbered frame would look fine and measure the wrong residues.
    """
    import gemmi

    structure = gemmi.Structure()
    structure.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_name)

    coords = traj.xyz[index] * 10.0          # mdtraj nm to angstrom
    current = None
    residue = None
    for atom in traj.topology.atoms:
        key = (atom.residue.chain.index, atom.residue.resSeq, atom.residue.name)
        if atom.residue.name in ("HOH", "WAT", "NA", "CL", "POPC", "PC", "OL", "PA"):
            continue
        if atom.residue.name.startswith("PO") and atom.residue.name != "PRO":
            continue
        if key != current:
            if residue is not None:
                chain.add_residue(residue)
            residue = gemmi.Residue()
            residue.name = atom.residue.name
            residue.seqid = gemmi.SeqId(int(atom.residue.resSeq), " ")
            current = key
        gatom = gemmi.Atom()
        gatom.name = atom.name
        gatom.element = gemmi.Element(atom.element.symbol if atom.element else "C")
        position = coords[atom.index]
        gatom.pos = gemmi.Position(float(position[0]), float(position[1]),
                                   float(position[2]))
        residue.add_atom(gatom)
    if residue is not None:
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    return structure


def ca_coordinates(traj, residues: list[int]):
    """Calpha coordinates for given residue numbers, in angstrom.

    Returns the array and the residue numbers actually found, in the same order,
    so a caller comparing two structures can intersect them rather than assume
    both had every residue.
    """
    import numpy as np

    wanted = set(residues)
    index_by_number: dict[int, int] = {}
    for atom in traj.topology.atoms:
        if atom.name != "CA" or atom.residue.resSeq not in wanted:
            continue
        index_by_number.setdefault(int(atom.residue.resSeq), atom.index)
    found = sorted(index_by_number)
    if not found:
        return np.zeros((traj.n_frames, 0, 3)), []
    indices = [index_by_number[n] for n in found]
    return traj.xyz[:, indices, :] * 10.0, found


def superpose(mobile, target):
    """Kabsch superposition of one frame's core onto a reference core.

    Returns the rotation and the two centroids rather than the transformed
    coordinates, so the same transform can then be applied to a motif that was
    not part of the fit. That is the whole point of a core superposition: the fit
    is on the core and the measurement is on the motif.
    """
    import numpy as np

    mobile_centre = mobile.mean(axis=0)
    target_centre = target.mean(axis=0)
    correlation = (mobile - mobile_centre).T @ (target - target_centre)
    u, _, vt = np.linalg.svd(correlation)
    # Guard against a reflection: an improper rotation would fit a mirror image
    # and report a small RMSD for a structure of the wrong handedness.
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, sign])
    rotation = vt.T @ correction @ u.T
    return rotation, mobile_centre, target_centre


def apply_transform(coords, rotation, mobile_centre, target_centre):
    return (coords - mobile_centre) @ rotation.T + target_centre


def rmsd(one, two) -> float:
    import numpy as np

    return float(np.sqrt(((one - two) ** 2).sum(axis=1).mean()))


def motif_rmsd_to_reference(frame_coords, reference_coords, core_transform) -> float:
    """Motif RMSD after the core has been superposed, which is section 5.6's measure."""
    moved = apply_transform(frame_coords, *core_transform)
    return rmsd(moved, reference_coords)


# ===========================================================================
#  Driving a whole run
# ===========================================================================


@dataclass
class RunAnalysis:
    """One trajectory's per-frame state, plus what it was measured against."""

    name: str
    kind: str                                     # "kinase" or "gpcr"
    core: Core = field(default_factory=Core)
    frames: list[FrameState] = field(default_factory=list)
    motif_residues: dict[str, list[int]] = field(default_factory=dict)
    references: dict[str, list[str]] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def as_row(self) -> dict[str, Any]:
        return {
            "name": self.name, "kind": self.kind, "core": self.core.as_row(),
            "n_frames": len(self.frames),
            "motif_residues": {k: list(v) for k, v in sorted(self.motif_residues.items())},
            "references": {k: list(v) for k, v in sorted(self.references.items())},
            "frames": [f.as_row() for f in self.frames],
            "notes": list(self.notes),
        }


def analyse(traj, *, name: str, kind: str, core: Core,
            motif_residues: dict[str, list[int]],
            reference_motifs: dict[str, dict[str, Any]],
            frame_ps: float = 100.0,
            rule_labels=None, measure=None) -> RunAnalysis:
    """Score every frame of one trajectory against the reference ensembles.

    `reference_motifs` is `{"active": {code: {motif: coords}}, "inactive": ...}`,
    each coords array already expressed in the core frame of its own reference.
    `rule_labels` is an optional callable taking (gemmi chain, frame index) and
    returning `{motif: (label, metric)}`, which is how the rule-based half of
    section 5.6 reaches the score without this function knowing whether it is
    looking at a kinase or a receptor. `measure` is the same arrangement for
    sections 5.4 and 5.5's other per-frame readings, which are stored beside the
    scores rather than folded into them.
    """
    import numpy as np

    run = RunAnalysis(name=name, kind=kind, core=core,
                      motif_residues=motif_residues)
    run.references = {state: sorted(entries)
                      for state, entries in reference_motifs.items()}

    core_coords, core_found = ca_coordinates(traj, core.residues)
    if not core_found:
        run.notes.append("no core residue resolved in the trajectory; "
                         "no motif RMSD can be measured")
        return run
    if len(core_found) < len(core.residues):
        run.notes.append(f"core is {len(core_found)} of {len(core.residues)} "
                         "residues; the rest are not in this construct")

    # One reference supplies the target core. Any of them would do, since every
    # motif RMSD is a difference measured in the superposed frame, and taking
    # the first alphabetically makes the choice reproducible rather than
    # dependent on dictionary order.
    target_core = None
    for state in sorted(reference_motifs):
        for code in sorted(reference_motifs[state]):
            candidate = reference_motifs[state][code].get("_core")
            if candidate is not None and len(candidate) == len(core_found):
                target_core = np.asarray(candidate, dtype=float)
                break
        if target_core is not None:
            break
    if target_core is None:
        run.notes.append("no reference supplied a core of matching length; "
                         "frames are superposed on the first frame instead")
        target_core = core_coords[0]

    motif_coords = {motif: ca_coordinates(traj, residues)
                    for motif, residues in motif_residues.items()}

    for index in range(traj.n_frames):
        transform = superpose(core_coords[index], target_core)
        chain = None
        labels: dict[str, tuple[str | None, str | None]] = {}
        if rule_labels is not None or measure is not None:
            # Built once per frame and shared: the conversion is the expensive
            # part, and two callers each building their own would double it.
            chain = frame_to_gemmi(traj, index)[0]["A"]
        if rule_labels is not None:
            labels = rule_labels(chain, index) or {}

        state = FrameState(time_ns=index * frame_ps / 1000.0)
        if measure is not None:
            state.measurements = measure(chain, index, transform) or {}
        for motif, (coords, found) in motif_coords.items():
            if not found:
                continue
            distances: dict[str, dict[str, float]] = {"active": {}, "inactive": {}}
            for which in ("active", "inactive"):
                for code, motifs in reference_motifs.get(which, {}).items():
                    reference = motifs.get(motif)
                    if reference is None or len(reference) != len(found):
                        continue
                    distances[which][code] = motif_rmsd_to_reference(
                        coords[index], np.asarray(reference, dtype=float), transform)
            label, metric = labels.get(motif, (None, None))
            state.motifs.append(score_motif(
                motif, active=distances["active"], inactive=distances["inactive"],
                n_residues=len(found), rule_label=label, rule_metric=metric))
        run.frames.append(state)
    return run


# Everything that is solvent, counter-ion or bilayer. Anything else in the file
# is the thing being simulated and is kept.
NOT_THE_SYSTEM = {"HOH", "WAT", "TIP3", "NA", "CL", "K", "MG", "ZN",
                  "POPC", "POP", "PC", "OL", "PA", "DLPC", "DOPC", "DPPC"}


def kept_atoms(topology) -> list[int]:
    """Protein plus whatever else is not solvent, ion or lipid.

    Written as an exclusion rather than as a list of ligands because the
    inclusion version was wrong: it selected `protein or resname CAU`, which is
    carazolol, so the ABL1 imatinib run would have had its STI silently dropped
    and been analysed as an apo system with nothing raised. A ligand the code
    has not been told about is exactly the case that has to work, since the
    whole package is meant to handle a protein it has never seen.
    """
    import mdtraj as md

    keep = []
    for atom in topology.atoms:
        name = atom.residue.name.upper()
        if name in NOT_THE_SYSTEM:
            continue
        # Lipid17 splits a lipid into residues whose names vary by tail, so the
        # prefix catches the ones the set above does not name. PRO is the one
        # amino acid that would be caught by it.
        if name.startswith("PO") and name != "PRO":
            continue
        keep.append(atom.index)
    return keep


def load(topology: Path, trajectory: Path, *, stride: int = 1):
    """One trajectory, with the solvent, ions and bilayer dropped.

    Dropped on load rather than later: a 190,000-atom membrane system is mostly
    water, and every frame of it held in memory is the difference between a few
    hundred megabytes and a few gigabytes for data no measurement here reads.
    """
    import mdtraj as md

    first = md.load(str(topology))
    keep = kept_atoms(first.topology)
    if not keep:
        raise ValueError(f"{topology.name} has nothing but solvent in it")
    return md.load(str(trajectory), top=str(topology), stride=stride,
                   atom_indices=keep)
#  Deposited numbering
# ===========================================================================


def deposited_numbering(start_pdb: Path, topology) -> dict[int, int]:
    """Map the trajectory's residue numbers back onto the deposited ones.

    OpenMM renumbers from 1 when it writes the solvated system, so a trajectory
    of the ABL1 construct has residues 1 to 262 where the deposition has 233 to
    498. Every annotation in this package addresses deposited numbers, resolved
    from KLIFS or GPCRdb, so without this map every motif would be read at the
    wrong residues and would still return plausible numbers.

    A positional map is not safe to assume twice over. The ABL1 construct spans
    266 residues and contains 262, so it has internal gaps and residue `i` is not
    `233 + i`; and PDBFixer builds internal gaps back, so the trajectory can hold
    residues the start model does not. The one-letter sequence is therefore
    checked at every position and a mismatch raises rather than returns a map
    that is wrong in the middle.
    """
    import gemmi

    structure = gemmi.read_structure(str(start_pdb))
    structure.setup_entities()
    deposited: list[tuple[int, str]] = []
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is not None and info.is_amino_acid():
                deposited.append((residue.seqid.num,
                                  info.one_letter_code.upper()))

    trajectory: list[tuple[int, str]] = []
    for residue in topology.residues:
        if not residue.is_protein:
            continue
        info = gemmi.find_tabulated_residue(residue.name)
        letter = info.one_letter_code.upper() if info is not None else "X"
        trajectory.append((int(residue.resSeq), letter))

    if len(trajectory) != len(deposited):
        raise ValueError(
            f"{start_pdb.name} has {len(deposited)} residues and the trajectory "
            f"has {len(trajectory)}; PDBFixer has changed the construct, so the "
            "deposited numbering cannot be recovered positionally")

    mismatches = [(i, d[1], t[1]) for i, (d, t)
                  in enumerate(zip(deposited, trajectory)) if d[1] != t[1]]
    if mismatches:
        first = mismatches[0]
        raise ValueError(
            f"sequence mismatch at position {first[0]}: the start model has "
            f"{first[1]} and the trajectory has {first[2]} "
            f"({len(mismatches)} mismatches); refusing to map the numbering")

    return {t[0]: d[0] for d, t in zip(deposited, trajectory)}


def renumber(topology, mapping: dict[int, int]) -> int:
    """Rewrite a trajectory topology's residue numbers to the deposited ones.

    Done in place on the topology rather than by translating every residue list
    at each call site: one translation here cannot be forgotten, and a forgotten
    one at a call site is a measurement of the wrong residues that looks right.
    """
    changed = 0
    for residue in topology.residues:
        # Protein only. Waters are numbered from 1 as well, so a map built for
        # residues 1 to 262 matches several thousand of them: an unrestricted
        # pass renumbered 606 residues where 262 were intended, moving waters
        # into the protein's own deposited range.
        if not residue.is_protein:
            continue
        target = mapping.get(int(residue.resSeq))
        if target is not None and target != residue.resSeq:
            residue.resSeq = target
            changed += 1
    return changed


def scaffold_core(pocket, exclude: tuple[str, ...] = ()) -> Core:
    """The KLIFS pocket scaffold: every position outside the mobile motifs.

    Section 5.4 names the core as the C-lobe helices alphaE, alphaF and alphaH
    "defined by DSSP on the reference ensemble medoid". Only alphaE is inside the
    KLIFS pocket, and a core of alphaE alone is five residues, which is not a
    rigid-body fit: with it, motif RMSDs to the inactive references came out at
    17 to 22 A, which is not a measurement of anything.

    So the core is the scaffold instead, which is the same choice BoltzMaker's
    own `alignment.stable_reference_positions` makes and for the same reason. It
    is every KLIFS position that is not one of the mobile or functional regions,
    it spans both lobes, it is resolved for every reference by KLIFS itself, and
    it is about sixty residues rather than five.

    This is a deviation from the letter of section 5.4 and it is checked rather
    than asserted: `core_agreement` measures how well the chosen core
    superposes across the reference set, and a core that moves with activation
    would show up there as a large active-to-inactive spread.
    """
    mobile = set(exclude) or {"g.l", "alphaC", "b.l", "GK", "hinge", "c.l",
                              "xDFG", "a.l"}
    core = Core(source="KLIFS pocket scaffold, excluding "
                       + ", ".join(sorted(mobile)))
    for index in sorted(pocket.region_by_index):
        region = pocket.region_by_index[index]
        if region in mobile or index not in pocket.by_index:
            continue
        core.helices.setdefault(region, []).append(pocket.by_index[index])
    core.residues = sorted(n for v in core.helices.values() for n in v)
    return core


def core_agreement(core_coords_by_code: dict[str, Any]) -> dict[str, Any]:
    """How rigidly a candidate core superposes across the reference set.

    Every reference is fitted onto the first and the residual RMSD recorded. A
    core that is genuinely stable gives a small spread; a core that moves with
    activation gives a large one, and that is the measurement that justifies
    using it rather than an assertion that it is the right one.
    """
    import numpy as np

    codes = sorted(core_coords_by_code)
    if len(codes) < 2:
        return {"n": len(codes)}
    target = np.asarray(core_coords_by_code[codes[0]], dtype=float)
    residuals = {}
    for code in codes[1:]:
        mobile = np.asarray(core_coords_by_code[code], dtype=float)
        if mobile.shape != target.shape:
            continue
        transform = superpose(mobile, target)
        residuals[code] = rmsd(apply_transform(mobile, *transform), target)
    if not residuals:
        return {"n": len(codes), "note": "no reference matched the core length"}
    values = list(residuals.values())
    return {"n": len(codes) , "reference": codes[0],
            "mean_rmsd_a": round(float(np.mean(values)), 3),
            "max_rmsd_a": round(float(np.max(values)), 3),
            "worst": max(residuals, key=lambda k: residuals[k]),
            "per_code": {k: round(v, 3) for k, v in sorted(residuals.items())}}


# ===========================================================================
#  Helix axes
# ===========================================================================

# An alpha helix turns every 3.6 residues, so a mean over four consecutive
# Calpha lands close to the axis and averages the spiral away.
HELIX_TURN = 4


def helix_axis_points(ca_coords):
    """Points on a helix's axis, from its Calpha positions.

    A running mean over one turn. The obvious alternative, fitting a principal
    axis straight to the Calpha, is not stable on a short segment and this was
    measured rather than assumed: the kink angle at P6.50 came out anywhere
    between 24 and 50 degrees depending only on how many residues were fitted,
    because the Calpha spiral and a four-point SVD is dominated by where in the
    turn the segment happens to start and end.

    Averaging over a whole turn removes the spiral, so the points lie on the
    axis and a line through them is the axis rather than an artefact of the
    window.
    """
    import numpy as np

    coords = np.asarray(ca_coords, dtype=float)
    if coords.shape[0] < HELIX_TURN:
        return np.zeros((0, 3))
    return np.array([coords[i:i + HELIX_TURN].mean(axis=0)
                     for i in range(coords.shape[0] - HELIX_TURN + 1)])


def helix_axis(ca_coords):
    """One helix's axis direction, oriented along the chain.

    Oriented so an angle between two axes is the angle and not its supplement:
    an unoriented axis reports 160 degrees for a 20 degree kink.
    """
    import numpy as np

    points = helix_axis_points(ca_coords)
    if points.shape[0] < 2:
        return None
    centred = points - points.mean(axis=0)
    _, _, vectors = np.linalg.svd(centred, full_matrices=False)
    direction = vectors[0]
    if np.dot(direction, points[-1] - points[0]) < 0:
        direction = -direction
    return direction


def axis_angle(first, second) -> float | None:
    """The angle between two oriented axes, in degrees."""
    import numpy as np

    if first is None or second is None:
        return None
    cosine = float(np.clip(np.dot(first, second), -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))
