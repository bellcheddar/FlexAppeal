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


# ===========================================================================
#  Micro-switches (spec section 5.5)
# ===========================================================================

# Section 5.5's thresholds, each from the measure it belongs to rather than
# shared: the ionic lock is a salt bridge and uses a salt-bridge distance, the
# sodium pocket is a coordination shell and uses a coordination distance.
IONIC_LOCK_A = 4.0
SODIUM_COORDINATION_A = 3.5


def _residue_by_number(chain, number: int):
    for residue in chain:
        if residue.seqid.num == number:
            return residue
    return None


def _heavy_atoms(residue, names: tuple[str, ...] = ()):
    """Named heavy atoms, or every heavy atom if no names are given."""
    import gemmi

    hydrogen = gemmi.Element("H")
    return [a for a in residue
            if a.element != hydrogen and (not names or a.name in names)]


def min_heavy_distance(chain, numbering: GenericNumbering, first: str, second: str,
                       *, first_atoms: tuple[str, ...] = (),
                       second_atoms: tuple[str, ...] = ()) -> float | None:
    """The closest approach between two residues' heavy atoms.

    A salt bridge is a side-chain contact, so it is measured between side-chain
    atoms and not between Calpha: the ionic lock can break while the two Calpha
    barely move, which is the whole reason it is reported separately from the
    TM3-TM6 distance.
    """
    numbers = [numbering.residue(first), numbering.residue(second)]
    if any(n is None for n in numbers):
        return None
    residues = [_residue_by_number(chain, n) for n in numbers]
    if any(r is None for r in residues):
        return None
    one = _heavy_atoms(residues[0], first_atoms)
    two = _heavy_atoms(residues[1], second_atoms)
    if not one or not two:
        return None
    return min(a.pos.dist(b.pos) for a in one for b in two)


def ionic_lock(chain, numbering: GenericNumbering) -> dict[str, Any]:
    """R3.50 to E6.30, the classic class A lock.

    Intact in the inactive receptor, so `intact` maps to inactive in the state
    vocabulary. Measured between the guanidinium and carboxylate heavy atoms
    rather than between the named CZ and CD of section 5.5, because both groups
    are symmetric and naming one atom of a symmetric pair makes the distance
    depend on which way the side chain was built.
    """
    distance = min_heavy_distance(chain, numbering, *IONIC_LOCK,
                                  first_atoms=("NE", "CZ", "NH1", "NH2"),
                                  second_atoms=("CD", "OE1", "OE2"))
    row: dict[str, Any] = {"pair": f"{numbering.label(IONIC_LOCK[0])}"
                                   f" - {numbering.label(IONIC_LOCK[1])}"}
    if distance is None:
        row["state"] = "unresolved"
        return row
    row["distance_a"] = round(distance, 2)
    row["state"] = "intact" if distance < IONIC_LOCK_A else "broken"
    return row


def chi_angles(chain, numbering: GenericNumbering, generic: str) -> dict[str, Any]:
    """chi1 and chi2 for one generic position, for the rotamer switches.

    The toggle tryptophan and Y7.53 are rotamer switches: the backbone hardly
    moves and the side chain flips, so a Calpha measure reports nothing and the
    dihedral reports everything.
    """
    import gemmi

    number = numbering.residue(generic)
    residue = _residue_by_number(chain, number) if number is not None else None
    row: dict[str, Any] = {"position": numbering.label(generic)}
    if residue is None:
        return row

    # chi1 is N-CA-CB-<first side-chain atom>, chi2 one bond further out. The
    # branch atom differs by residue type, so it is taken from the residue rather
    # than assumed: CG for most, OG/OG1/SG for the small polar ones.
    chi1_third = next((n for n in ("CG", "CG1", "OG", "OG1", "SG", "CG2")
                       if residue.find_atom(n, "*") is not None), None)
    chi2_third = next((n for n in ("CD", "CD1", "OD1", "ND1", "SD", "CD2")
                       if residue.find_atom(n, "*") is not None), None)

    def dihedral(*names: str) -> float | None:
        atoms = [residue.find_atom(n, "*") for n in names]
        if any(a is None for a in atoms):
            return None
        return gemmi.calculate_dihedral(*[a.pos for a in atoms]) * 180.0 / 3.141592653589793

    if chi1_third:
        value = dihedral("N", "CA", "CB", chi1_third)
        if value is not None:
            row["chi1_deg"] = round(value, 1)
    if chi1_third and chi2_third:
        value = dihedral("CA", "CB", chi1_third, chi2_third)
        if value is not None:
            row["chi2_deg"] = round(value, 1)
    return row


def sodium_pocket(chain, numbering: GenericNumbering,
                  ions: list | None = None) -> dict[str, Any]:
    """D2.50 to S3.39 and to N7.49, plus whether a sodium is coordinated.

    A sodium in the pocket favours the inactive state, and its absence is not
    evidence either way: the ion is often simply not modelled, and in a
    simulation it has to diffuse in. So `sodium_present` is reported only when
    there are ion positions to look at, and omitted otherwise rather than
    written as false.
    """
    row: dict[str, Any] = {}
    for label, pair in (("d_s3x39_a", ("2.50", "3.39")),
                        ("d_n7x49_a", ("2.50", "7.49"))):
        distance = min_heavy_distance(chain, numbering, *pair)
        if distance is not None:
            row[label] = round(distance, 2)

    if ions is None:
        return row
    number = numbering.residue("2.50")
    residue = _residue_by_number(chain, number) if number is not None else None
    if residue is None:
        return row
    carboxylate = _heavy_atoms(residue, ("OD1", "OD2"))
    if not carboxylate:
        return row
    nearest = min((min(a.pos.dist(p) for a in carboxylate) for p in ions),
                  default=None)
    if nearest is not None:
        row["nearest_sodium_a"] = round(nearest, 2)
        row["sodium_present"] = nearest < SODIUM_COORDINATION_A
    return row


def micro_switches(chain, numbering: GenericNumbering,
                   ions: list | None = None) -> dict[str, Any]:
    """Every section 5.5 measure that one structure can give on its own.

    The reference-dependent ones (TM6 rotation, the NPxxY and connector RMSDs to
    the medoids) are not here: they need the reference ensemble and belong with
    the state score in `state.py`. What is here is what a single frame answers by
    itself, which is also what can be checked against a crystal structure.
    """
    cleft = tm3_tm6_distance(chain, numbering)
    alternate = ca_distance(chain, numbering, *TM3_TM6_ALT)
    row: dict[str, Any] = {
        "ionic_lock": ionic_lock(chain, numbering),
        "toggle_w6x48": chi_angles(chain, numbering, TOGGLE[0]),
        "npxxy_y7x53": chi_angles(chain, numbering, NPXXY[1]),
        "sodium_pocket": sodium_pocket(chain, numbering, ions),
    }
    if cleft is not None:
        row["tm3_tm6_a"] = round(cleft, 2)
        row["tm3_tm6_pair"] = (f"{numbering.label(TM3_TM6_OPEN[0])}"
                               f" - {numbering.label(TM3_TM6_OPEN[1])}")
    if alternate is not None:
        row["tm3_tm6_alt_a"] = round(alternate, 2)

    # Y5.58 OH to Y7.53 OH: the two tyrosines come together on activation, and
    # this is the one NPxxY measure that needs no reference.
    bridge = min_heavy_distance(chain, numbering, *NPXXY_TO_TM5,
                                first_atoms=("OH",), second_atoms=("OH",))
    if bridge is not None:
        row["y5x58_y7x53_oh_a"] = round(bridge, 2)
    return row


# ===========================================================================
#  The rest of section 5.5
# ===========================================================================

# Section 5.5's TM6 rotation is of the *intracellular half*, because that is the
# end that swings out to open the G-protein site; the extracellular half barely
# moves and including it would average the signal away.
TM6_INTRACELLULAR = ("6.28", "6.40")
# Kink angles are measured over several window sizes rather than one, and the
# median is reported with the spread beside it.
#
# One window was not defensible. With a principal axis fitted straight to the
# Calpha the angle at P6.50 came out anywhere from 24 to 50 degrees depending
# only on the window, which is a parameter choice masquerading as a
# measurement. Turn-averaging the axis fixed most of it, and the rest is that
# a short window simply has too few axis points to fit a line through: below
# seven residues there are fewer than four. So the window is swept and the
# answer carries its own uncertainty.
KINK_WINDOWS = (7, 8, 9, 10, 11, 12)
# Fewer axis points than this and the line through them is noise.
KINK_MIN_AXIS_POINTS = 4


def _ca_array(chain, numbers: list[int]):
    import numpy as np

    found = []
    for number in numbers:
        residue = _residue_by_number(chain, number)
        atom = residue.find_atom("CA", "*") if residue is not None else None
        if atom is not None:
            found.append([atom.pos.x, atom.pos.y, atom.pos.z])
    return np.asarray(found, dtype=float)


def generic_span(numbering: GenericNumbering, helix: str,
                 first: str, last: str) -> list[int]:
    """The residue numbers between two generic positions on one helix."""
    index = helix[2:] if helix.startswith("TM") else helix
    lo = int(float(first) * 100) % 100
    hi = int(float(last) * 100) % 100
    found = []
    for offset in range(min(lo, hi), max(lo, hi) + 1):
        number = numbering.residue(f"{index}.{offset}")
        if number is not None:
            found.append(number)
    return sorted(found)


def tm6_rotation(chain, reference_chain, numbering: GenericNumbering,
                 core_transform=None) -> dict:
    """Axis rotation and centroid shift of TM6's intracellular half.

    The same measurement as section 5.4's alphaC displacement, so it reuses that
    implementation rather than carrying a second copy: one helix against the
    same helix in a reference, after the core superposition.
    """
    from .kinase import helix_displacement

    numbers = generic_span(numbering, "TM6", *TM6_INTRACELLULAR)
    if not numbers:
        return {"status": "TM6 span not resolved"}
    found = helix_displacement(chain, reference_chain, numbers, core_transform)
    found["span"] = [min(numbers), max(numbers)]
    found["generic"] = list(TM6_INTRACELLULAR)
    return found


def connector_rmsd(chain, reference_chain, numbering: GenericNumbering,
                   core_transform=None) -> dict:
    """The PIF connector: I3.40, P5.50 and F6.44 side chains against a reference.

    Side chains rather than Calpha, because the connector rearranges by
    repacking three side chains against each other while their backbones stay
    put. Measured on every heavy atom the three residues share with the
    reference, so a residue modelled with a different number of atoms is
    dropped rather than compared on a subset.
    """
    import numpy as np

    from .traj import apply_transform, rmsd

    numbers = [numbering.residue(g) for g in PIF]
    if any(n is None for n in numbers):
        return {"status": "PIF positions not resolved"}

    mine, theirs, used = [], [], []
    for generic, number in zip(PIF, numbers):
        one = _residue_by_number(chain, number)
        two = _residue_by_number(reference_chain, number)
        if one is None or two is None:
            continue
        names = ({a.name for a in _heavy_atoms(one)}
                 & {a.name for a in _heavy_atoms(two)})
        # Backbone excluded: the point of the connector is the side chains.
        names -= {"N", "CA", "C", "O", "OXT"}
        if not names:
            continue
        for name in sorted(names):
            first = one.find_atom(name, "*")
            second = two.find_atom(name, "*")
            if first is None or second is None:
                continue
            mine.append([first.pos.x, first.pos.y, first.pos.z])
            theirs.append([second.pos.x, second.pos.y, second.pos.z])
        used.append(f"{generic} ({number})")

    if len(mine) < 3:
        return {"status": "too few shared side-chain atoms"}
    moved = np.asarray(mine, dtype=float)
    if core_transform is not None:
        moved = apply_transform(moved, *core_transform)
    return {"status": "ok", "positions": used, "n_atoms": len(mine),
            "rmsd_a": round(rmsd(moved, np.asarray(theirs, dtype=float)), 2)}


def helix8_angle(chain, numbering: GenericNumbering) -> dict:
    """Helix 8's angle to the membrane plane.

    H8 lies along the membrane rather than through it, so the informative angle
    is to the plane and not to the normal: a well-behaved H8 reads near zero and
    one that has lifted off reads larger. The bilayer normal is z by
    construction, the receptor having been oriented onto OPM, so this is
    90 degrees minus the angle to z.
    """
    import numpy as np

    numbers = numbering.segment_residues("H8")
    coords = _ca_array(chain, numbers)
    if coords.shape[0] < 4:
        return {"status": "H8 not resolved"}
    from .traj import helix_axis

    axis = helix_axis(coords)
    if axis is None:
        return {"status": "H8 axis not determined"}
    to_normal = np.degrees(np.arccos(min(1.0, abs(float(axis[2])))))
    return {"status": "ok", "n_residues": int(coords.shape[0]),
            "angle_to_membrane_plane_deg": round(90.0 - float(to_normal), 1),
            "mean_z_a": round(float(coords[:, 2].mean()), 1)}


def kink_angles(chain, numbering: GenericNumbering) -> dict:
    """Kink angle at each conserved proline (section 5.5's TM kinks).

    The angle between the helix axis fitted either side of the proline, swept
    over several window sizes because a single window is a parameter choice
    rather than a measurement. The median is the figure and the spread says how
    far to trust it: a straight helix reads near zero, and the conserved
    prolines of TM6 and TM7 hold real kinks of which the one at 6.50 changes
    with activation.
    """
    import statistics as st

    from .traj import axis_angle, helix_axis, helix_axis_points

    out: dict = {}
    for generic in KINK_PROLINES:
        number = numbering.residue(generic)
        if number is None:
            out[generic] = {"status": "not resolved"}
            continue

        angles = []
        for half in KINK_WINDOWS:
            before = _ca_array(chain, list(range(number - half, number)))
            after = _ca_array(chain, list(range(number + 1, number + 1 + half)))
            if (helix_axis_points(before).shape[0] < KINK_MIN_AXIS_POINTS
                    or helix_axis_points(after).shape[0] < KINK_MIN_AXIS_POINTS):
                continue
            angle = axis_angle(helix_axis(before), helix_axis(after))
            if angle is not None:
                angles.append(angle)

        if not angles:
            out[generic] = {"status": "too few residues either side"}
            continue
        median = st.median(angles)
        spread = max(angles) - min(angles)
        out[generic] = {
            "status": "ok", "residue": number,
            "label": numbering.label(generic),
            "kink_deg": round(median, 1),
            "spread_deg": round(spread, 1),
            "n_windows": len(angles),
            # A spread of more than half the value means the windows do not
            # agree on the answer, so the figure is reported and not claimed.
            "resolved": spread < 0.5 * median if median > 0 else False,
        }
    return out
