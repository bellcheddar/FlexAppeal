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


def kinase_core(pocket, dssp_helices: dict[str, list[int]] | None = None) -> Core:
    """The C-lobe helices alphaE, alphaF and alphaH.

    Only alphaE is inside the KLIFS pocket, so it comes from the pocket map and
    the other two come from `dssp_helices`, which a caller derives by running
    DSSP on the ensemble medoid. Without them the core is alphaE alone, and
    `source` says so rather than pretending to three helices.
    """
    core = Core()
    alpha_e = pocket.region("alphaE")
    if alpha_e:
        core.helices["alphaE"] = sorted(alpha_e)
    for name in ("alphaF", "alphaH"):
        found = (dssp_helices or {}).get(name)
        if found:
            core.helices[name] = sorted(found)
    core.residues = sorted(n for v in core.helices.values() for n in v)
    core.source = ("KLIFS alphaE plus DSSP alphaF and alphaH"
                   if len(core.helices) == 3
                   else "KLIFS alphaE only; alphaF and alphaH were not resolved")
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
            rule_labels=None) -> RunAnalysis:
    """Score every frame of one trajectory against the reference ensembles.

    `reference_motifs` is `{"active": {code: {motif: coords}}, "inactive": ...}`,
    each coords array already expressed in the core frame of its own reference.
    `rule_labels` is an optional callable taking (gemmi chain, frame index) and
    returning `{motif: (label, metric)}`, which is how the rule-based half of
    section 5.6 reaches the score without this function knowing whether it is
    looking at a kinase or a receptor.
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
        if rule_labels is not None:
            chain = frame_to_gemmi(traj, index)[0]["A"]
            labels = rule_labels(chain, index) or {}

        state = FrameState(time_ns=index * frame_ps / 1000.0)
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


def load(topology: Path, trajectory: Path, *, stride: int = 1):
    """One trajectory, protein and ligand only, with the solvent dropped.

    The solvent is dropped on load rather than later: a 190,000-atom membrane
    system is mostly water, and every frame of it held in memory is the
    difference between a few hundred megabytes and a few gigabytes for data
    that no measurement here looks at.
    """
    import mdtraj as md

    first = md.load(str(topology))
    keep = first.topology.select("protein or resname CAU")
    return md.load(str(trajectory), top=str(topology), stride=stride,
                   atom_indices=keep)
