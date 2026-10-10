"""Class A GPCR state calls, resolved through GPCRdb generic numbering.

Every residue here is addressed by its generic number and looked up in the
catalogue's own residue table. Nothing is typed: `R3.50` becomes R131 in β2AR
because GPCRdb says so, and would become a different number in any other
receptor, which is the point of the scheme and the reason a module written this
way works for a receptor it has never seen.

The verification that this is wired up correctly is cheap and worth stating: the
spec's own worked example is "R3.50 (R131)", and the mapping built from the live
table returns exactly that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Spec section 5.5. Each motif is named by the generic numbers it is read from,
# so the definitions survive a change of receptor.
TM3_TM6_OPEN = ("3.50", "6.34")      # the intracellular cleft
TM3_TM6_ALT = ("3.50", "6.30")       # reported beside it
IONIC_LOCK = ("3.50", "6.30")        # R3.50 to E6.30
NPXXY = ("7.49", "7.53")
NPXXY_TO_TM5 = ("5.58", "7.53")      # Y5.58 OH to Y7.53 OH
PIF = ("3.40", "5.50", "6.44")
TOGGLE = ("6.48",)
SODIUM_POCKET = ("2.50", "3.39", "7.49")
KINK_PROLINES = ("6.50", "7.50")


@dataclass
class GenericNumbering:
    """One receptor's map between generic numbers and its own residue numbers."""

    entry_name: str
    by_generic: dict[str, int] = field(default_factory=dict)
    by_residue: dict[int, str] = field(default_factory=dict)
    residue_letter: dict[int, str] = field(default_factory=dict)
    segment: dict[int, str] = field(default_factory=dict)

    def residue(self, generic: str) -> int | None:
        return self.by_generic.get(generic)

    def label(self, generic: str) -> str:
        """`R3.50 (R131)`, the form section 5.5 asks every GPCR label to take."""
        number = self.by_generic.get(generic)
        if number is None:
            return generic
        letter = self.residue_letter.get(number, "")
        return f"{letter}{generic} ({letter}{number})"

    def segment_residues(self, name: str) -> list[int]:
        return sorted(n for n, seg in self.segment.items() if seg == name)


def build_numbering(rows: list[dict[str, Any]],
                    entry_name: str = "adrb2_human",
                    scheme: str = "BW") -> GenericNumbering:
    """The map, from GPCRdb's extended residue table.

    Ballesteros-Weinstein by default because that is what the labels use, but the
    scheme is a parameter: GPCRdb's own structure-based numbering is in
    `display_generic_number` and is the better choice for anything that has to
    align across receptors with different insertions.
    """
    numbering = GenericNumbering(entry_name=entry_name)
    for row in rows:
        number = row.get("sequence_number")
        if number is None:
            continue
        numbering.residue_letter[number] = row.get("amino_acid", "")
        numbering.segment[number] = row.get("protein_segment", "")
        for alt in (row.get("alternative_generic_numbers") or []):
            if alt.get("scheme") == scheme and alt.get("label"):
                label = str(alt["label"])
                # First wins: a generic number maps to one residue, and a later
                # duplicate would mean the table is not what this assumes.
                numbering.by_generic.setdefault(label, number)
                numbering.by_residue.setdefault(number, label)
    return numbering


# ===========================================================================
#  Per-structure measurement
# ===========================================================================


def ca_distance(chain, numbering: GenericNumbering, first: str, second: str
                ) -> float | None:
    """Cα to Cα between two generic positions, or None if either is unresolved.

    The residue identity is checked against the catalogue rather than trusted.
    Most β2AR entries use the native numbering and 3.50 really is residue 131,
    but engineered constructs renumber, and a distance measured between the
    wrong two residues is a number that looks fine.
    """
    wanted = {}
    for generic in (first, second):
        number = numbering.residue(generic)
        if number is None:
            return None
        wanted[number] = numbering.residue_letter.get(number, "")

    import gemmi

    found: dict[int, Any] = {}
    for residue in chain:
        if residue.seqid.num not in wanted:
            continue
        info = gemmi.find_tabulated_residue(residue.name)
        letter = info.one_letter_code.upper() if info else ""
        if letter != wanted[residue.seqid.num].upper():
            return None                      # renumbered: refuse rather than guess
        atom = residue.find_atom("CA", "*")
        if atom is None:
            return None
        found[residue.seqid.num] = atom.pos
    if len(found) != 2:
        return None
    (one, two) = list(found.values())
    return one.dist(two)


def tm3_tm6_distance(chain, numbering: GenericNumbering) -> float | None:
    """The intracellular cleft opening: the headline GPCR activation measure.

    Larger on activation. TM6 swings outward to open the G-protein site, and this
    distance is what the Dror 2011 state maps are drawn against.
    """
    return ca_distance(chain, numbering, *TM3_TM6_OPEN)


def fit_threshold(active: list[float], inactive: list[float]) -> dict[str, Any]:
    """The cleft distance that best separates the two states, leave-one-out.

    Section 9's P3. The threshold is fitted on the reference set and scored with
    each entry left out in turn, so the number reported is not the number the cut
    was chosen to flatter.
    """
    values = [(v, "active") for v in active] + [(v, "inactive") for v in inactive]
    if not active or not inactive:
        return {"threshold": None, "accuracy": None, "n": len(values)}

    def best_cut(sample: list[tuple[float, str]]) -> float:
        points = sorted({v for v, _ in sample})
        candidates = [(a + b) / 2 for a, b in zip(points, points[1:])] or points
        scored = []
        for cut in candidates:
            correct = sum(1 for v, state in sample
                          if (v >= cut) == (state == "active"))
            scored.append((correct, cut))
        return max(scored)[1]

    correct = 0
    for index in range(len(values)):
        held_out = values[index]
        rest = values[:index] + values[index + 1:]
        cut = best_cut(rest)
        if (held_out[0] >= cut) == (held_out[1] == "active"):
            correct += 1

    return {
        "threshold": round(best_cut(values), 2),
        "accuracy": correct / len(values),
        "correct": correct,
        "n": len(values),
        "n_active": len(active),
        "n_inactive": len(inactive),
        "active_median": round(sorted(active)[len(active) // 2], 2),
        "inactive_median": round(sorted(inactive)[len(inactive) // 2], 2),
    }
