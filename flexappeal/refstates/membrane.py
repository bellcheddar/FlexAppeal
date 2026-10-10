"""Did the bilayer behave? Spec section 5.8.

A membrane simulation can look healthy in every protein measurement while the
bilayer underneath it has thinned, frozen or let water into its core, and every
state call made on that trajectory would inherit the problem without showing it.
So the bilayer is measured on its own terms, per 10 ns chunk, and shown beside
the protein panels rather than assumed.

The expected ranges are POPC at 303 K with Lipid17, and they are the spec's:

    area per lipid      60 to 68 A^2 expected, outside 55 to 75 is G6
    P-P thickness       36 to 41 A expected
    core water          a count, as a sanity metric rather than a threshold
    receptor tilt       against the bilayer normal

Area per lipid is the box cross-section divided by the lipids in one leaflet, so
it needs the leaflet split rather than the total: dividing by all the lipids
would report half the real area and fail a check that should pass.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Section 5.8. The expected window is what POPC should give; the G6 window is
# how far out it has to be before the run is not usable.
APL_EXPECTED_A2 = (60.0, 68.0)
APL_ACCEPTABLE_A2 = (55.0, 75.0)
THICKNESS_EXPECTED_A = (36.0, 41.0)
CHUNK_NS = 10.0

# A water is "in the core" if its oxygen sits within this of the bilayer midplane.
# Half the expected hydrophobic thickness, so it counts waters that are genuinely
# between the acyl chains rather than ones hydrating a headgroup.
CORE_HALF_WIDTH_A = 8.0

LIPID_RESIDUES = ("POPC", "POP", "PC", "OL", "PA")
PHOSPHORUS = "P"


@dataclass
class MembraneChunk:
    """One 10 ns window of bilayer health."""

    start_ns: float
    end_ns: float
    n_frames: int = 0
    apl_a2: float | None = None
    apl_sd_a2: float | None = None
    thickness_a: float | None = None
    thickness_sd_a: float | None = None
    box_z_a: float | None = None
    core_waters: float | None = None
    tilt_deg: float | None = None
    n_lipids: int = 0
    leaflet_counts: tuple[int, int] = (0, 0)
    notes: list[str] = field(default_factory=list)

    @property
    def apl_verdict(self) -> str | None:
        if self.apl_a2 is None:
            return None
        if APL_EXPECTED_A2[0] <= self.apl_a2 <= APL_EXPECTED_A2[1]:
            return "expected"
        if APL_ACCEPTABLE_A2[0] <= self.apl_a2 <= APL_ACCEPTABLE_A2[1]:
            return "outside expected"
        return "G6"

    @property
    def thickness_verdict(self) -> str | None:
        if self.thickness_a is None:
            return None
        lo, hi = THICKNESS_EXPECTED_A
        return "expected" if lo <= self.thickness_a <= hi else "outside expected"

    def as_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {"start_ns": round(self.start_ns, 2),
                               "end_ns": round(self.end_ns, 2),
                               "n_frames": self.n_frames,
                               "n_lipids": self.n_lipids,
                               "leaflet_counts": list(self.leaflet_counts)}
        for key, value in (("apl_a2", self.apl_a2), ("apl_sd_a2", self.apl_sd_a2),
                           ("thickness_a", self.thickness_a),
                           ("thickness_sd_a", self.thickness_sd_a),
                           ("box_z_a", self.box_z_a),
                           ("core_waters", self.core_waters),
                           ("tilt_deg", self.tilt_deg)):
            if value is not None:
                row[key] = round(value, 2)
        for key, value in (("apl_verdict", self.apl_verdict),
                           ("thickness_verdict", self.thickness_verdict)):
            if value is not None:
                row[key] = value
        if self.notes:
            row["notes"] = list(self.notes)
        return row


def phosphorus_indices(topology) -> list[int]:
    """Lipid phosphorus atoms, which mark the two headgroup planes.

    Selected by element within a lipid residue rather than by atom name, because
    Lipid17 splits POPC into PC, OL and PA residues and the headgroup name
    differs between that and a single-residue POPC.
    """
    found = []
    for atom in topology.atoms:
        if atom.element is None or atom.element.symbol != PHOSPHORUS:
            continue
        if atom.residue.name in LIPID_RESIDUES or atom.residue.name.startswith("PO"):
            found.append(atom.index)
    return found


def split_leaflets(z_values) -> tuple[list[int], list[int]]:
    """Which phosphorus sits in which leaflet, split at their own midplane.

    Split at the median rather than at z = 0: the bilayer drifts in z over a run
    under a membrane barostat, and a split at the origin would put a whole
    leaflet on the wrong side once it had moved, halving the apparent area per
    lipid in exactly the chunk where the box changed.
    """
    import numpy as np

    midplane = float(np.median(z_values))
    upper = [i for i, z in enumerate(z_values) if z >= midplane]
    lower = [i for i, z in enumerate(z_values) if z < midplane]
    return upper, lower


def receptor_tilt(ca_coords) -> float | None:
    """Angle between the receptor's principal axis and the bilayer normal.

    The principal axis of the TM Calpha cloud is the helix-bundle axis for a
    7TM receptor, and the normal is z by construction, so the tilt is the angle
    between that axis and z. Folded to 0 to 90 degrees because the axis has no
    direction: a bundle tilted 20 degrees and one tilted 160 are the same tilt.
    """
    import numpy as np

    coords = np.asarray(ca_coords, dtype=float)
    if coords.shape[0] < 3:
        return None
    centred = coords - coords.mean(axis=0)
    _, _, vectors = np.linalg.svd(centred, full_matrices=False)
    axis = vectors[0]
    cosine = abs(float(axis[2])) / float(np.linalg.norm(axis))
    return float(np.degrees(np.arccos(min(1.0, cosine))))


def analyse_chunk(traj, *, start_ns: float, end_ns: float,
                  tm_ca_indices: list[int] | None = None) -> MembraneChunk:
    """Bilayer health for one chunk of trajectory.

    `traj` is an mdtraj Trajectory already sliced to the chunk, in nanometres as
    mdtraj stores things; everything reported here is in angstrom, which is the
    unit the rest of this build uses.
    """
    import numpy as np

    chunk = MembraneChunk(start_ns=start_ns, end_ns=end_ns,
                          n_frames=int(traj.n_frames))
    phosphorus = phosphorus_indices(traj.topology)
    chunk.n_lipids = len(phosphorus)
    if not phosphorus:
        chunk.notes.append("no lipid phosphorus found; not a membrane system")
        return chunk
    if traj.unitcell_lengths is None:
        chunk.notes.append("trajectory carries no box; area per lipid needs one")
        return chunk

    phos_z = traj.xyz[:, phosphorus, 2] * 10.0          # nm to angstrom
    box = traj.unitcell_lengths * 10.0

    # Per frame, because the leaflet split and the box both move during a chunk.
    apl, thickness = [], []
    upper_n, lower_n = [], []
    for frame in range(traj.n_frames):
        upper, lower = split_leaflets(phos_z[frame])
        if not upper or not lower:
            continue
        upper_n.append(len(upper))
        lower_n.append(len(lower))
        area = float(box[frame, 0] * box[frame, 1])
        # One leaflet's worth, and the mean of the two rather than one of them,
        # so an uneven split does not bias the number towards whichever leaflet
        # happened to be chosen.
        apl.append(area / (0.5 * (len(upper) + len(lower))))
        thickness.append(float(np.mean(phos_z[frame][upper])
                               - np.mean(phos_z[frame][lower])))

    if apl:
        chunk.apl_a2 = float(np.mean(apl))
        chunk.apl_sd_a2 = float(np.std(apl))
        chunk.thickness_a = float(np.mean(thickness))
        chunk.thickness_sd_a = float(np.std(thickness))
        chunk.box_z_a = float(np.mean(box[:, 2]))
        chunk.leaflet_counts = (int(round(np.mean(upper_n))),
                                int(round(np.mean(lower_n))))
        if abs(chunk.leaflet_counts[0] - chunk.leaflet_counts[1]) > 2:
            chunk.notes.append(
                f"leaflets differ by {abs(chunk.leaflet_counts[0] - chunk.leaflet_counts[1])} "
                "lipids, which biases area per lipid")

    water_o = traj.topology.select("water and name O")
    if len(water_o):
        midplane = np.median(phos_z, axis=1)
        water_z = traj.xyz[:, water_o, 2] * 10.0
        inside = np.abs(water_z - midplane[:, None]) < CORE_HALF_WIDTH_A
        chunk.core_waters = float(inside.sum(axis=1).mean())

    if tm_ca_indices:
        tilts = [receptor_tilt(traj.xyz[frame, tm_ca_indices, :] * 10.0)
                 for frame in range(traj.n_frames)]
        measured = [t for t in tilts if t is not None]
        if measured:
            chunk.tilt_deg = float(np.mean(measured))
    return chunk


def gate_six(chunks: list[MembraneChunk]) -> dict[str, Any]:
    """Whether the bilayer result is usable, in G6's terms.

    One chunk out of range is not a failed run: the first chunk after
    equilibration often is, and the spec's gate is about a result being
    invalidated rather than about a single window. So G6 is raised only if the
    production chunks are out of the acceptable window on average.
    """
    import numpy as np

    scored = [c for c in chunks if c.apl_a2 is not None]
    if not scored:
        return {"open": False, "reason": "no chunk produced an area per lipid"}
    mean_apl = float(np.mean([c.apl_a2 for c in scored]))
    bad = [c for c in scored if c.apl_verdict == "G6"]
    opens = not (APL_ACCEPTABLE_A2[0] <= mean_apl <= APL_ACCEPTABLE_A2[1])
    return {
        "open": opens,
        "mean_apl_a2": round(mean_apl, 2),
        "acceptable_a2": list(APL_ACCEPTABLE_A2),
        "n_chunks": len(scored),
        "n_chunks_out_of_range": len(bad),
        "reason": (f"mean area per lipid {mean_apl:.1f} A^2 is outside "
                   f"{APL_ACCEPTABLE_A2[0]} to {APL_ACCEPTABLE_A2[1]}")
                  if opens else "mean area per lipid is within the acceptable window",
    }
