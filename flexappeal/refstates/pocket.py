"""KLIFS's own 85-position pocket map, per structure.

Spec section 5.4 asks for every kinase residue to be "resolved from KLIFS
positions" with "KLIFS-style region names in the output", and the vendored
annotator does not go that far: it returns six motifs and lumps the remaining
sixty positions into one `pocket_scaffold`. That is enough to superpose on and
not enough to name a P-loop, a beta-4 strand or an activation loop.

`interactions_match_residues` gives the whole map directly, as KLIFS labels it:

    index 1..85   KLIFS_position   Xray_position
        4         g.l.4            249

So the region names here are KLIFS's, the boundaries are KLIFS's, and the residue
numbers are the deposited ones. Nothing in this module types a residue number or
a region boundary, which is the property section 5.4 is really asking for: the
same code names the regions of a kinase it has never seen.

The nineteen regions, in sequence order, are I, g.l, II, III, alphaC, b.l, IV, V,
GK, hinge, linker, alphaD, alphaE, VI, c.l, VII, VIII, xDFG and a.l. They are
discovered from the response rather than listed, so a change at KLIFS shows up as
a new region name rather than as a silently dropped one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from .sources import KLIFS_ROOT, Catalogue, SourceError

MATCH_RESIDUES = f"{KLIFS_ROOT}/interactions_match_residues"
STRUCTURE_LIST = f"{KLIFS_ROOT}/structure_list"

# KLIFS writes a position as `<region>.<index>`, and for the two single-residue
# regions as `GK.45` and `VIII.79` in the same form. A label that does not parse
# is kept whole rather than guessed at.
_LABEL = re.compile(r"^(?P<region>.+)\.(?P<index>\d+)$")

# Section 5.4's R-spine, named the way the kinase literature names it. Each entry
# is resolved below from the pocket map rather than from a residue number.
R_SPINE_ORDER = ("beta4", "alphaC_plus4", "DFG_Phe", "HRD_His")

# Section 5.4 puts this at 6.0 A. Measured on ABL1's own reference set, the
# DFG-Phe to HRD-His centroid distance separates the states completely:
#
#     active    5.33 to 5.89 A, median 5.70   (12 of 12 inside 6.0)
#     inactive  9.53 to 11.08 A, median 9.86  (0 of 45 inside 6.0)
#
# So 6.0 A is correct for crystal structures, and it sits a tenth of an angstrom
# above the top of the active range, which leaves no room for a side chain that
# is breathing. Trajectory frames of the apo DFG-in system read 6.03 to 6.69 A,
# plainly on the active side of a 3.6 A gap, and a 6.0 A line called all of them
# broken.
#
# The threshold is therefore placed in the empty region between the clusters
# rather than at the edge of one, which is the same correction D010 made to the
# DFG criterion and for the same reason. A crystal structure is one minimised
# conformation; an MD frame at 300 K is not.
R_SPINE_CONTACT_A = 7.7
R_SPINE_CONTACT_CRYSTAL_A = 6.0

# The C-lobe helices section 5.4 names as the superposition core. Only alphaE of
# the three is inside the KLIFS pocket, so the other two are found by DSSP on the
# medoid; see `core.py`.
KINASE_CORE_REGIONS = ("alphaE",)

# KLIFS spells alphaC and alphaD with the Greek letter. Normalised so a caller
# can write `regions["alphaC"]` without carrying the encoding around.
_GREEK = {"αC": "alphaC", "αD": "alphaD", "αE": "alphaE"}


@dataclass
class PocketMap:
    """One structure's KLIFS pocket, addressable by position or by region."""

    structure_id: int
    pdb: str = ""
    chain: str = ""
    by_index: dict[int, int] = field(default_factory=dict)
    label_by_index: dict[int, str] = field(default_factory=dict)
    region_by_index: dict[int, str] = field(default_factory=dict)
    unresolved: list[int] = field(default_factory=list)

    def region(self, name: str) -> list[int]:
        """The deposited residue numbers of one region, in sequence order."""
        return [self.by_index[i] for i in sorted(self.region_by_index)
                if self.region_by_index[i] == name and i in self.by_index]

    def regions(self) -> dict[str, list[int]]:
        out: dict[str, list[int]] = {}
        for index in sorted(self.region_by_index):
            out.setdefault(self.region_by_index[index], [])
            if index in self.by_index:
                out[self.region_by_index[index]].append(self.by_index[index])
        return out

    def at(self, index: int) -> int | None:
        return self.by_index.get(index)

    def label(self, residue: int) -> str | None:
        """`g.l.4` for a residue number, for labelling a plot or a table."""
        for index, number in self.by_index.items():
            if number == residue:
                return self.label_by_index.get(index)
        return None


def pocket_map(catalogue: Catalogue, structure_id: int) -> PocketMap:
    """The 85-position map for one KLIFS structure.

    A position KLIFS could not resolve comes back with an `Xray_position` of `-1`
    or a blank, and lands in `unresolved` rather than being written as a residue
    number. Section 0.1's rule about never substituting a stand-in applies to
    residue numbers as much as to scores: a `-1` here would be looked up as a
    real residue somewhere downstream and would silently measure nothing.
    """
    rows = catalogue.get(MATCH_RESIDUES, {"structure_ID": structure_id})
    if not isinstance(rows, list) or not rows:
        raise SourceError(f"KLIFS returned no pocket map for structure {structure_id}")

    found = PocketMap(structure_id=int(structure_id))
    for row in rows:
        index = int(row["index"])
        label = str(row.get("KLIFS_position", "")).strip()
        found.label_by_index[index] = label
        match = _LABEL.match(label)
        region = match.group("region") if match else label
        found.region_by_index[index] = _GREEK.get(region, region)

        raw = str(row.get("Xray_position", "")).strip()
        # KLIFS marks a gap as -1, and an insertion code arrives attached to the
        # number ("317A"), which is a residue this build cannot address by
        # seqid alone.
        if not raw or raw == "-1" or not raw.lstrip("-").isdigit():
            found.unresolved.append(index)
            continue
        number = int(raw)
        if number < 0:
            found.unresolved.append(index)
            continue
        found.by_index[index] = number
    return found


def structure_ids(catalogue: Catalogue, kinase_id: int) -> dict[str, int]:
    """PDB code to KLIFS structure id, so a pocket map can be fetched by code.

    Where a code has several KLIFS structures, as a multi-chain entry does, the
    first is taken and the choice is visible in the returned map rather than
    buried: the chain is part of `PocketMap` for exactly this reason.
    """
    rows = catalogue.get(f"{KLIFS_ROOT}/structures_list", {"kinase_ID": kinase_id})
    out: dict[str, int] = {}
    for row in rows or []:
        code = str(row.get("pdb", "")).upper()
        if code:
            out.setdefault(code, int(row["structure_ID"]))
    return out


# ===========================================================================
#  The regions section 5.4 needs, resolved from the map
# ===========================================================================


def activation_loop(pocket: PocketMap, chain=None,
                    last_residue: int | None = None) -> list[int]:
    """The activation loop, from the DFG aspartate to the APE motif.

    KLIFS's `a.l` region is two residues long: the pocket alignment stops where
    the loop leaves the pocket, so taking `a.l` as the activation loop gives five
    residues for a loop that is about twenty, and every A-loop RMSD measured on
    it would be an RMSD of the DFG motif under another name.

    So KLIFS gives the start, as the aspartate of `xDFG` past the X, and the end
    is the APE motif, which is where the activation loop is defined to end. APE
    is found in the chain's own sequence downstream of DFG rather than typed,
    and `ANE`/`APQ`-style variants are not searched for: if APE is absent the
    span falls back to `a.l` and the caller can see the loop is short.
    """
    xdfg = pocket.region("xDFG")
    if not xdfg:
        return []
    start = xdfg[1] if len(xdfg) > 1 else xdfg[0]      # the D of DFG, past the X

    end = None
    if chain is not None:
        end = _ape_residue(chain, start)
    if end is None:
        tail = pocket.region("a.l")
        end = max(tail) if tail else start
    if last_residue is not None:
        end = min(end, last_residue)
    return list(range(start, end + 1)) if end >= start else []


# The activation loop runs from DFG to the APE motif, and nothing is searched
# beyond this far downstream: in every protein kinase the loop is well inside it,
# and an unbounded search would find an unrelated APE in the C-lobe.
APE_SEARCH_LIMIT = 40


def _ape_residue(chain, start: int) -> int | None:
    """The glutamate of the first APE triplet downstream of `start`.

    The loop ends at the E, which is the residue the span should include.
    """
    import gemmi

    letters: dict[int, str] = {}
    for residue in chain:
        info = gemmi.find_tabulated_residue(residue.name)
        if info is not None and info.is_amino_acid():
            letters[residue.seqid.num] = info.one_letter_code.upper()

    for number in range(start, start + APE_SEARCH_LIMIT):
        window = "".join(letters.get(number + offset, "-") for offset in range(3))
        if window == "APE":
            return number + 2
    return None


def r_spine(pocket: PocketMap, chain=None) -> dict[str, int]:
    """The four regulatory-spine residues, resolved from the pocket map.

    Three of the four are positional and come straight from KLIFS: the HRD
    histidine opens `c.l`, the DFG phenylalanine is the middle of `xDFG` past the
    X, and the alphaC anchor is the alphaC glutamate plus four.

    The beta-4 leucine is the one that is not positional. It sits in KLIFS region
    `IV`, which is four residues long, and which of them carries the spine is a
    matter of the fold rather than of the alignment. With a `chain` it is read as
    the region IV residue whose side chain is nearest the alphaC anchor, which is
    what being on the spine means. Without one the region is returned unresolved
    rather than guessed, so a caller that has no coordinates gets three residues
    and knows it.
    """
    spine: dict[str, int] = {}
    catalytic = pocket.region("c.l")
    if catalytic:
        spine["HRD_His"] = catalytic[0]
    xdfg = pocket.region("xDFG")
    if len(xdfg) > 2:
        spine["DFG_Phe"] = xdfg[2]
    alpha_c = pocket.region("alphaC")
    glu = alpha_c[4] if len(alpha_c) > 8 else None
    if glu is not None:
        spine["alphaC_plus4"] = glu + 4

    beta4 = pocket.region("IV")
    if beta4 and chain is not None and "alphaC_plus4" in spine:
        anchor = _side_chain_centroid(chain, spine["alphaC_plus4"])
        if anchor is not None:
            scored = []
            for number in beta4:
                centroid = _side_chain_centroid(chain, number)
                if centroid is not None:
                    scored.append((anchor.dist(centroid), number))
            if scored:
                spine["beta4"] = min(scored)[1]
    return spine


def _side_chain_centroid(chain, number: int):
    """The mean position of a residue's side-chain heavy atoms.

    Backbone excluded because the spine is a side-chain stack: measuring it on
    Calpha would report the helix packing rather than the spine, and the two come
    apart exactly in the structures this is meant to tell apart.
    """
    import gemmi

    for residue in chain:
        if residue.seqid.num != number:
            continue
        atoms = [a for a in residue
                 if a.name not in ("N", "CA", "C", "O", "OXT")
                 and a.element != gemmi.Element("H")]
        if not atoms:
            # Glycine, or a side chain that was not resolved. Calpha is the
            # honest fallback and the caller is told by the empty side chain.
            atoms = [a for a in residue if a.name == "CA"]
        if not atoms:
            return None
        total = gemmi.Position(0, 0, 0)
        for atom in atoms:
            total += gemmi.Position(atom.pos.x, atom.pos.y, atom.pos.z)
        return gemmi.Position(total.x / len(atoms), total.y / len(atoms),
                              total.z / len(atoms))
    return None


def r_spine_state(chain, spine: dict[str, int]) -> dict[str, Any]:
    """`assembled` if every consecutive spine pair is within R_SPINE_CONTACT_A.

    Consecutive along the spine, not along the sequence: the order is beta4,
    alphaC+4, DFG-Phe, HRD-His, which is how the four side chains stack in space.
    A missing residue makes the state `unresolved` rather than `broken`, because
    an unmeasured pair is not an open one.

    The threshold is 7.7 A and not the published 6.0: D015 moved it into the
    empty gap between the crystal clusters, because 6.0 sits a tenth of an
    angstrom above the top of the active range and left no room for a side chain
    at 300 K. `tools/build_references.py` re-derives the separation it was fitted
    to, so the number can be checked rather than taken on trust.
    """
    centroids = {name: _side_chain_centroid(chain, number)
                 for name, number in spine.items()}
    pairs: list[dict[str, Any]] = []
    complete = True
    for first, second in zip(R_SPINE_ORDER, R_SPINE_ORDER[1:]):
        one, two = centroids.get(first), centroids.get(second)
        if one is None or two is None:
            complete = False
            pairs.append({"pair": f"{first}-{second}"})
            continue
        distance = one.dist(two)
        pairs.append({"pair": f"{first}-{second}", "distance_a": round(distance, 2),
                      "contact": distance < R_SPINE_CONTACT_A})

    if not complete:
        state = "unresolved"
    elif all(p.get("contact") for p in pairs):
        state = "assembled"
    else:
        state = "broken"
    return {"state": state, "pairs": pairs,
            "threshold_a": R_SPINE_CONTACT_A,
            "residues": {k: v for k, v in sorted(spine.items())}}
