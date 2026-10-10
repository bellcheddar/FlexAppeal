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
        """What the vendored 11.0 A rule would have said, for comparison."""
        if self.d_alphac_a is None or self.d_lys_a is None:
            return None
        cutoff = 11.0
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
