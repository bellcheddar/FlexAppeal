"""Kinase state calls, wrapping the vendored classifier rather than editing it.

One function here, and it exists because the vendored absolute threshold
misclassified a third of ABL1's own reference structures.

`sse_comparison.metrics.classify_dfg_state` measures the DFG-Phe ring against two
anchors, the αC-Glu+4 Cα and the catalytic Lys Cα, and then asks whether each
distance falls the correct side of 11.0 Å. The measurement is right and is reused
verbatim. The threshold is what fails:

    DFG-in  references: ring 6.6 to 6.8 A from alphaC, 12.9 to 14.0 from Lys
    DFG-out references: ring 13.1 to 16.3 A from alphaC, 9.1 to 11.3 from Lys

Six DFG-out structures, 2HYY and 3PYY among them, sit at 11.1 to 11.3 A from the
Lys: a tenth of an angstrom past a line drawn at 11.0, and so neither "in" nor
"out" but "other". KLIFS calls all six DFG-out, and they plainly are.

The classes separate by **8.5 A** on the difference between the two distances
(DFG-in clusters at -6.1 to -7.9, DFG-out at +2.4 to +5.8) and there is nothing in
between. So the state is read from which anchor the ring is nearer, which needs no
threshold at the place the data is dense. Agreement with KLIFS goes from 11/17 to
17/17.

A margin is kept, but it is placed where the data is empty rather than where it is
dense. Within `AMBIGUOUS_MARGIN_A` of equidistant the answer is "other", which
matters for trajectory frames rather than references: a DFG mid-flip genuinely is
neither state, and a bare sign test would hide the transition by snapping every
frame to one side.

Logged as **D010** and as a candidate upstream change for BoltzMaker: the same
threshold is in its own `compare-sse` and will mislabel the same structures there.
"""

from __future__ import annotations

from dataclasses import dataclass

# Half-width of the band around equidistant where neither state is claimed. The
# observed class gap is 8.5 A, so this is an order of magnitude inside the empty
# region: it cannot reach either cluster, and it still lets a flipping torsion
# report honestly.
AMBIGUOUS_MARGIN_A = 1.0


@dataclass
class DFGCall:
    """A DFG state, with the two distances it was read from.

    Both distances travel with the call so a reader can see the geometry rather
    than trust the label, and so the absolute-threshold answer stays recoverable
    from stored data without re-running anything.
    """

    state: str | None
    d_alphac_a: float | None
    d_lys_a: float | None

    @property
    def separation_a(self) -> float | None:
        if self.d_alphac_a is None or self.d_lys_a is None:
            return None
        return self.d_alphac_a - self.d_lys_a

    @property
    def absolute_threshold_state(self) -> str | None:
        """What the vendored published rule would have said, for comparison."""
        if self.d_alphac_a is None or self.d_lys_a is None:
            return None
        cutoff = DFG_PUBLISHED_A
        if self.d_alphac_a <= cutoff < self.d_lys_a:
            return "in"
        if self.d_lys_a <= cutoff < self.d_alphac_a:
            return "out"
        return "other"


def dfg_distances(frame, structure, fam_to_pos, motifs_by_name) -> tuple:
    """The two anchor distances, measured by the vendored code.

    Nothing is reimplemented here: `_atom_pos` and `_residue_at` are
    `sse_comparison`'s own, and they carry its handling of missing atoms and of
    the ring atom falling back through CZ, CE1, CD1.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vendor" / "boltzmaker"))
    from sse_comparison import metrics

    dfg = motifs_by_name.get("DFG")
    lys = motifs_by_name.get("catalytic_Lys")
    glu = motifs_by_name.get("alphaC_Glu")
    if not (dfg and lys and glu):
        return None, None

    # KLIFS pocket positions 81, 82, 83 are Asp, Phe, Gly; the ring is the middle.
    phe = dfg.residues[1] if len(dfg.residues) > 1 else dfg.residues[0]
    ring = metrics._atom_pos(
        metrics._residue_at(structure, fam_to_pos, phe), "CZ", "CE1", "CD1")
    near = metrics._atom_pos(
        metrics._residue_at(structure, fam_to_pos, glu.residues[0] + 4), "CA")
    far = metrics._atom_pos(
        metrics._residue_at(structure, fam_to_pos, lys.residues[0]), "CA")
    if ring is None or near is None or far is None:
        return None, None
    return ring.dist(near), ring.dist(far)


def classify_dfg(frame, motifs_by_name, *, half: str = "apo") -> DFGCall:
    """DFG-in or DFG-out from which anchor the ring is nearer.

    `half` picks which side of a comparison frame to read. A reference is built
    against itself so both halves are the same structure; a trajectory frame is
    built against its own topology the same way.
    """
    structure = frame.apo if half == "apo" else frame.holo
    mapping = frame.fam_to_apo if half == "apo" else frame.fam_to_holo
    d_alphac, d_lys = dfg_distances(frame, structure, mapping, motifs_by_name)
    if d_alphac is None or d_lys is None:
        return DFGCall(None, None, None)

    separation = d_alphac - d_lys
    if abs(separation) < AMBIGUOUS_MARGIN_A:
        state = "other"
    else:
        state = "out" if separation > 0 else "in"
    return DFGCall(state, d_alphac, d_lys)


# ===========================================================================
#  The rest of section 5.4
# ===========================================================================

# Section 5.4. A salt bridge is a side-chain contact and this is the distance
# the kinase literature uses for the alphaC-in call.
ALPHAC_SALT_BRIDGE_A = 4.0
# The published DFG two-distance cut, kept for comparison rather than used:
# D010 replaced it with the nearer-anchor rule because it called six of ABL1's
# own DFG-out references neither state. Named because it was written out as
# 11.0 here and again in the state map's `published_threshold`, and a figure
# that draws a line somebody else's rule is at should draw the line that rule
# actually is.
DFG_PUBLISHED_A = 11.0
# Section 5.4 names 30 degrees as BoltzMaker's own phi/psi flag threshold.
PHI_PSI_FLAG_DEG = 30.0


def _residue(chain, number: int):
    for residue in chain:
        if residue.seqid.num == number:
            return residue
    return None


def _atoms(residue, names: tuple[str, ...]):
    if residue is None:
        return []
    return [a for a in residue if a.name in names]


def alphac_salt_bridge(chain, pocket) -> dict:
    """Catalytic Lys NZ to the alphaC Glu carboxylate, and the state it implies.

    Measured to the nearer carboxylate oxygen rather than to a named one: the
    group is symmetric, so naming OE1 makes the distance depend on which way the
    side chain was built rather than on where it is.

    `unresolved` rather than `out` when an atom is missing, because an
    unmeasured salt bridge is not a broken one. That distinction is the whole
    reason this returns a dict rather than a float.
    """
    catalytic = pocket.region("c.l")
    lys_positions = [pocket.at(17)]
    glu_positions = [pocket.at(24)]
    lys = next((_residue(chain, n) for n in lys_positions if n), None)
    glu = next((_residue(chain, n) for n in glu_positions if n), None)

    nz = _atoms(lys, ("NZ",))
    carboxylate = _atoms(glu, ("OE1", "OE2"))
    # Named as well as numbered, the way the receptor's measurements are. The
    # numbers were already here and resolved from KLIFS positions 17 and 24,
    # never typed, and the page had nothing to show but "alphaC salt bridge".
    # A reader could not tell which two residues it was between.
    row: dict = {"lys": pocket.at(17), "glu": pocket.at(24)}
    if row["lys"] is not None and row["glu"] is not None:
        row["pair"] = f"Lys{row['lys']} to Glu{row['glu']}"
    if not nz or not carboxylate:
        row["state"] = "unresolved"
        return row
    distance = min(a.pos.dist(b.pos) for a in nz for b in carboxylate)
    row["distance_a"] = round(distance, 2)
    row["state"] = "in" if distance < ALPHAC_SALT_BRIDGE_A else "out"
    return row


def _ca_array(chain, numbers: list[int]):
    import numpy as np

    found = []
    for number in numbers:
        residue = _residue(chain, number)
        atom = residue.find_atom("CA", "*") if residue is not None else None
        if atom is not None:
            found.append([atom.pos.x, atom.pos.y, atom.pos.z])
    return np.asarray(found, dtype=float)


def _axis(coords):
    """One helix's axis, turn-averaged rather than fitted to raw Calpha.

    A principal axis straight from the Calpha is not stable on a short helix:
    see `traj.helix_axis_points` for the measurement that showed it.
    """
    from .traj import helix_axis

    return helix_axis(coords)


def core_fit(chain, reference_chain, core_numbers: list[int]):
    """The superposition that puts one chain in another's frame, on the core.

    Fitted here rather than taken from the caller, and that distinction cost a
    debugging round: the per-frame transform the scoring uses maps a frame into
    the *reference ensemble's* frame, while a comparison against the start model
    needs the start model's frame. Passing the first where the second was wanted
    gave a 78 angstrom centroid shift and an 82 angstrom connector RMSD, which
    were at least obviously wrong rather than plausibly so.
    """
    from .traj import superpose

    mine = _ca_array(chain, core_numbers)
    theirs = _ca_array(reference_chain, core_numbers)
    if mine.shape[0] < 3 or mine.shape != theirs.shape:
        return None
    return superpose(mine, theirs)


def helix_displacement(chain, reference_chain, numbers: list[int],
                       core_transform=None) -> dict:
    """Centroid shift and axis rotation of one helix against a reference.

    Section 5.4's alphaC displacement, and section 5.5's TM6 rotation, are the
    same measurement on different helices, so one function does both.

    `core_transform` must be a fit of this chain onto `reference_chain` on the
    structural core, which `core_fit` above produces. Without it this reports
    the whole protein drifting rather than the helix moving, which is the error
    the core exists to prevent.
    """
    import numpy as np

    from .traj import apply_transform

    mine = _ca_array(chain, numbers)
    theirs = _ca_array(reference_chain, numbers)
    if mine.shape[0] < 3 or mine.shape != theirs.shape:
        return {"status": "not resolved in both"}
    if core_transform is not None:
        mine = apply_transform(mine, *core_transform)

    shift = float(np.linalg.norm(mine.mean(axis=0) - theirs.mean(axis=0)))
    one, two = _axis(mine), _axis(theirs)
    row: dict = {"status": "ok", "n_residues": int(mine.shape[0]),
                 "centroid_shift_a": round(shift, 2)}
    if one is not None and two is not None:
        cosine = abs(float(np.dot(one, two)))
        # Folded to 0 to 90: a principal axis has no direction, so a helix
        # rotated 170 degrees and one rotated 10 are the same rotation.
        row["axis_rotation_deg"] = round(
            float(np.degrees(np.arccos(min(1.0, cosine)))), 1)
    return row


def p_loop(chain, pocket) -> dict:
    """The glycine-rich loop, by the fallback section 5.4 allows.

    Section 5.4 asks for the KLIFS G-rich distance, angle and rotation "if the
    definitions in the KLIFS literature can be verified at build time;
    otherwise P-loop tip Calpha to hinge distance plus P-loop Calpha RMSD (log
    which)". The KLIFS API exposes the pocket positions and their region names
    but not those three geometric definitions, and nothing in the responses
    pins down which residues the angle and rotation are measured between. So
    this is the fallback, and `definition` says so in the output rather than in
    a comment nobody reads.
    """
    import numpy as np

    g_rich = pocket.region("g.l")
    hinge = pocket.region("hinge")
    row: dict = {"definition": "fallback: tip to hinge distance",
                 "n_g_rich": len(g_rich), "n_hinge": len(hinge)}
    if not g_rich or not hinge:
        row["status"] = "not resolved"
        return row

    # The tip is the middle of the loop, which is the part that closes over the
    # nucleotide; the ends are anchored in the beta sheet and barely move.
    tip = g_rich[len(g_rich) // 2]
    tip_ca = _ca_array(chain, [tip])
    hinge_ca = _ca_array(chain, hinge)
    if tip_ca.shape[0] == 0 or hinge_ca.shape[0] == 0:
        row["status"] = "not resolved"
        return row
    row["status"] = "ok"
    row["tip_residue"] = tip
    row["tip_to_hinge_a"] = round(
        float(np.linalg.norm(tip_ca[0] - hinge_ca.mean(axis=0))), 2)
    return row


def backbone_flags(chain, reference_chain, numbers: list[int],
                   threshold_deg: float = PHI_PSI_FLAG_DEG) -> dict:
    """Residues whose backbone phi or psi has moved more than the threshold.

    Section 5.4's hinge measure, at BoltzMaker's own 30 degree threshold.
    Reported as which residues moved rather than as a count, because a hinge
    where one residue flipped and one where three did are different events and a
    count cannot tell them apart.
    """
    import gemmi

    def dihedrals(target, number: int):
        previous, current, following = (_residue(target, number - 1),
                                        _residue(target, number),
                                        _residue(target, number + 1))
        if current is None:
            return None, None
        def pos(residue, name):
            if residue is None:
                return None
            atom = residue.find_atom(name, "*")
            return atom.pos if atom is not None else None

        n, ca, c = pos(current, "N"), pos(current, "CA"), pos(current, "C")
        previous_c, next_n = pos(previous, "C"), pos(following, "N")
        phi = (gemmi.calculate_dihedral(previous_c, n, ca, c) * 180.0 / np.pi
               if None not in (previous_c, n, ca, c) else None)
        psi = (gemmi.calculate_dihedral(n, ca, c, next_n) * 180.0 / np.pi
               if None not in (n, ca, c, next_n) else None)
        return phi, psi

    import numpy as np

    def delta(one, two):
        if one is None or two is None:
            return None
        # Circular: 179 and -179 are two degrees apart, not 358.
        gap = abs(one - two) % 360.0
        return min(gap, 360.0 - gap)

    moved, measured = [], 0
    for number in numbers:
        mine = dihedrals(chain, number)
        theirs = dihedrals(reference_chain, number)
        for label, index in (("phi", 0), ("psi", 1)):
            change = delta(mine[index], theirs[index])
            if change is None:
                continue
            measured += 1
            if change > threshold_deg:
                moved.append({"residue": number, "angle": label,
                              "change_deg": round(change, 1)})
    return {"threshold_deg": threshold_deg, "n_measured": measured,
            "n_flagged": len(moved), "flagged": moved[:12]}


def dfg_plane(chain, pocket) -> dict:
    """The two Dunbrack distances, from a chain and its KLIFS pocket map.

    The same geometry `dfg_distances` reads through the vendored comparison
    frame, computed directly here because a trajectory frame has no comparison
    frame to read it through: the state map needs these per frame to draw the
    trajectory's path across the plane its references are plotted on.

    The ring atom falls back through CZ, CE1, CD1 exactly as the vendored
    `_atom_pos` does, so the two routes agree on a structure where both can run.
    """
    xdfg = pocket.region("xDFG")
    catalytic = pocket.at(17)
    glu = pocket.at(24)
    if len(xdfg) < 3 or catalytic is None or glu is None:
        return {"status": "positions not resolved"}

    phe = _residue(chain, xdfg[2])
    ring = None
    for name in ("CZ", "CE1", "CD1"):
        found = _atoms(phe, (name,))
        if found:
            ring = found[0]
            break
    # The alphaC anchor is the glutamate plus four, which is how the criterion
    # is defined rather than the glutamate itself.
    near = _atoms(_residue(chain, glu + 4), ("CA",))
    far = _atoms(_residue(chain, catalytic), ("CA",))
    if ring is None or not near or not far:
        return {"status": "atoms missing"}

    d_alphac = ring.pos.dist(near[0].pos)
    d_lys = ring.pos.dist(far[0].pos)
    separation = d_alphac - d_lys
    return {
        "status": "ok",
        "d_alphac_a": round(d_alphac, 2),
        "d_lys_a": round(d_lys, 2),
        "separation_a": round(separation, 2),
        # D010's rule: which anchor the ring is nearer, with the margin placed
        # where the reference classes leave an 8.5 A gap.
        "state": ("other" if abs(separation) < AMBIGUOUS_MARGIN_A
                  else "out" if separation > 0 else "in"),
    }
