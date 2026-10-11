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


def row_verdicts(apl_a2: float | None,
                 thickness_a: float | None) -> dict[str, str | None]:
    """The per-chunk verdicts for a chunk read back out of the store.

    Section 7 gives the `membrane` table six columns and no verdict, so the
    windows a chunk was judged against are not saved with it. They do not need
    to be: a verdict is a pure function of the number and the window. The
    panel recomputes it here, through `MembraneChunk` itself, so that the
    windows have one definition and the explorer cannot drift from what P6
    checks.
    """
    chunk = MembraneChunk(start_ns=0.0, end_ns=0.0, apl_a2=apl_a2,
                          thickness_a=thickness_a)
    return {"apl": chunk.apl_verdict, "thickness": chunk.thickness_verdict}


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


def centre_bilayer(z_values, box_z: float):
    """Shift z so the bilayer sits in the middle of the box, not across its edge.

    A trajectory stores wrapped coordinates, and under a membrane barostat the
    bilayer drifts until it straddles the periodic boundary. When it does, the
    phosphorus z distribution is bimodal at *both* ends of the box with nothing
    in the middle, and every measurement taken from it is wrong in a way that
    looks plausible:

        measured      thickness 78.0 A, 8,622 waters in the hydrophobic core
        actually      37.6 A, which is inside the expected window, and none

    The 78 A is the long way round the box between two leaflets that are 37.6 A
    apart, and the core water count is bulk water sitting at a midplane that is
    not where the midplane is.

    The centre is a circular mean rather than an arithmetic one, because the
    arithmetic mean of two clusters at opposite ends of a periodic axis lands in
    the empty middle, which is exactly the value that caused the problem.

    Returns the shifted z values and the shift applied, so a caller can put
    other atoms into the same frame.
    """
    import numpy as np

    angles = 2.0 * np.pi * np.asarray(z_values, dtype=float) / box_z
    centre = box_z * np.arctan2(np.sin(angles).mean(),
                                np.cos(angles).mean()) / (2.0 * np.pi)
    shift = box_z / 2.0 - centre
    return np.mod(np.asarray(z_values, dtype=float) + shift, box_z), float(shift)


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


# A rod's principal axis is only its axis if the cloud is actually rod shaped.
# Below this ratio of the first to the second singular value the direction is
# set by noise rather than by the helix, and averaging such an axis into a
# bundle makes the bundle worse rather than better.
MIN_AXIS_ANISOTROPY = 2.0
MIN_HELIX_CA = 8
MIN_BUNDLE_HELICES = 3


def helix_axis(ca_coords):
    """The unit axis of one helix, pointing towards +z, or None.

    One helix is a rod, so its Calpha cloud has a well-determined principal
    axis. Returned pointing towards +z rather than as found, because
    consecutive TM helices run in opposite directions and a mean taken over
    raw axes would cancel instead of averaging.
    """
    import numpy as np

    coords = np.asarray(ca_coords, dtype=float)
    if coords.shape[0] < MIN_HELIX_CA:
        return None
    centred = coords - coords.mean(axis=0)
    _, values, vectors = np.linalg.svd(centred, full_matrices=False)
    # A second singular value of zero is a perfectly straight rod, which is the
    # best axis there is rather than the worst. Only a non-zero one can be too
    # close to the first, so the ratio is the test and zero is not a rejection.
    if values[1] > 0 and values[0] / values[1] < MIN_AXIS_ANISOTROPY:
        return None
    axis = vectors[0] / np.linalg.norm(vectors[0])
    return -axis if axis[2] < 0 else axis


def axis_tilt_deg(ca_coords) -> float | None:
    """Angle between one rod's principal axis and the bilayer normal.

    The normal is z by construction, so this is the angle to z, folded to 0 to
    90 degrees because an axis has no direction: a rod tilted 20 degrees and
    one tilted 160 are the same tilt.

    For a single helix this is the tilt. For a whole receptor it is not, which
    is what `bundle_tilt` is for.
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


def bundle_tilt(ca_coords, residue_indices,
                *, min_helices: int = MIN_BUNDLE_HELICES) -> float | None:
    """Receptor tilt, as the mean direction of the transmembrane helix axes.

    Not the principal axis of the whole Calpha cloud, which is what this used
    to be and which does not work. A 7TM bundle inside a 37 A bilayer is about
    as wide as it is tall: measured on beta2AR the membrane-spanning Calpha
    span 40 by 45 by 37 A with singular values 159, 139, 109, a first-to-second
    ratio of 1.15. There is no principal direction in a cloud that round, so
    the axis it returns is noise, and it reported 25 degrees of tilt for a
    receptor sitting almost square in the membrane.

    Section 4.3's fallback orientation already gives the right definition, the
    mean direction of the seven TM helix axes, and each helix on its own is a
    rod whose axis is well determined: the same structure gives ratios of 4.8
    to 7.8 per helix and a bundle tilt of 6.7 degrees.

    `residue_indices` runs parallel to `ca_coords` and is what separates the
    helices: a helix crossing the bilayer is a contiguous run of residues, so
    the runs are the helices and no segment table is needed.
    """
    import numpy as np

    coords = np.asarray(ca_coords, dtype=float)
    residues = [int(r) for r in residue_indices]
    if coords.shape[0] != len(residues):
        raise ValueError(
            f"{coords.shape[0]} coordinates against {len(residues)} residue "
            f"indices; they must run in parallel")

    runs: list[list[int]] = []
    for position, residue in enumerate(residues):
        if runs and residue == residues[position - 1] + 1:
            runs[-1].append(position)
        else:
            runs.append([position])

    axes = [axis for axis in (helix_axis(coords[run]) for run in runs)
            if axis is not None]
    if len(axes) < min_helices:
        return None
    mean = np.mean(axes, axis=0)
    norm = float(np.linalg.norm(mean))
    if norm == 0.0:
        return None
    cosine = abs(float(mean[2])) / norm
    return float(np.degrees(np.arccos(min(1.0, cosine))))


def slab_ca_indices(topology, protein_z, midplane: float,
                    half_thickness: float) -> list[int]:
    """The protein Calpha that sit between the two phosphate planes.

    `bundle_tilt` wants the transmembrane bundle and only that. Every protein
    Calpha would do instead if a receptor were only its bundle, but it is not:
    beta2AR's ICL3 and its termini stick well out of the bilayer, and a
    principal axis fitted through them is the axis of the whole molecule rather
    than of the seven helices.

    The bilayer locates the bundle better than a residue table does. A Calpha
    inside the headgroup planes is membrane-spanning by construction, so the
    set follows the system that was built rather than a published segment range
    that may not match the construct, and it needs no catalogue lookup in a
    function that is already streaming a trajectory.

    `protein_z` is expected in angstrom and in the same frame the bilayer was
    centred into, which is why the midplane and the half thickness are passed
    in rather than measured again here.
    """
    import numpy as np

    offsets = np.abs(np.asarray(protein_z, dtype=float) - midplane)
    return [int(i) for i, offset in enumerate(offsets) if offset <= half_thickness]


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

    raw_phos_z = traj.xyz[:, phosphorus, 2] * 10.0       # nm to angstrom
    box = traj.unitcell_lengths * 10.0

    # Centred per frame before anything is measured. Without this a bilayer that
    # has drifted across the periodic boundary reports a thickness that is the
    # long way round the box.
    phos_z = np.empty_like(raw_phos_z)
    shifts = np.empty(traj.n_frames)
    for frame in range(traj.n_frames):
        phos_z[frame], shifts[frame] = centre_bilayer(
            raw_phos_z[frame], float(box[frame, 2]))

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
        gap = float(np.mean(phos_z[frame][upper])
                    - np.mean(phos_z[frame][lower]))
        height = float(box[frame, 2])
        thickness.append(min(abs(gap), height - abs(gap)))

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
        # The same shift the bilayer got, so the water is in the frame the
        # midplane is measured in. Counting it against an uncentred midplane
        # reported 8,622 waters in a core that had none.
        water_z = traj.xyz[:, water_o, 2] * 10.0
        shifted = np.mod(water_z + shifts[:, None], box[:, 2][:, None])
        midplane = np.median(phos_z, axis=1)
        separation = np.abs(shifted - midplane[:, None])
        # Periodic in z: a water just past the far edge is near the midplane
        # the short way round.
        separation = np.minimum(separation, box[:, 2][:, None] - separation)
        chunk.core_waters = float((separation < CORE_HALF_WIDTH_A).sum(axis=1).mean())

    if tm_ca_indices is None and chunk.thickness_a:
        # Derived here rather than required from the caller. This argument was
        # optional and nothing passed it, so every chunk ever written had no
        # tilt at all and the panel's column was a dash on a bilayer that was
        # otherwise healthy. A metric the caller has to remember to switch on
        # is a metric that gets forgotten, so the default is to measure it.
        ca = traj.topology.select("protein and name CA")
        if len(ca):
            first = np.mod(traj.xyz[0, ca, 2] * 10.0 + shifts[0],
                           float(box[0, 2]))
            tm_ca_indices = [int(ca[i]) for i in slab_ca_indices(
                traj.topology, first, float(np.median(phos_z[0])),
                0.5 * chunk.thickness_a)]
            if len(tm_ca_indices) < 3:
                chunk.notes.append(
                    "fewer than three Calpha inside the bilayer, so no tilt")

    if tm_ca_indices:
        residues = [traj.topology.atom(i).residue.index for i in tm_ca_indices]
        tilts = [bundle_tilt(traj.xyz[frame, tm_ca_indices, :] * 10.0, residues)
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
