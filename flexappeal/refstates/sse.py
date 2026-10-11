"""Secondary structure per frame, and how far it agrees with each reference set.

Spec section 5.6. Two separate things, and they answer different questions:

    per motif      helix fraction, strand fraction and where the element starts
                   and ends, which is how a motif unwinding shows up
    per residue    the fraction of active references and of inactive references
                   that share this frame's SSE class at this residue

The second is what makes a local unwinding interpretable. A frame whose NPxxY
backbone has lost its helical turn might be drifting towards the active
references, which also lose it, or it might be doing something neither reference
set does, and only a comparison against both tells those apart.

`compute_dssp(simplified=False)` throughout, so the eight DSSP classes survive
into storage. The simplified three-letter form is derived where it is wanted
rather than computed, because going the other way is impossible and the
distinction between a 3-10 helix and an alpha helix is exactly the kind of local
change this is looking for.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

# DSSP's eight classes. Mapped to the simplified three only where a caller asks,
# and the mapping is mdtraj's own so a comparison against its simplified output
# cannot drift from it.
HELIX_CLASSES = ("H", "G", "I")
STRAND_CLASSES = ("E", "B")
# How far apart the two agreements must be before the lean means anything.
# Chosen from the observed spread rather than picked: the discriminating motifs
# separate by 0.07 to 0.15 and the ones that carry no signal separate by 0.02
# to 0.04, so the line goes between them.
DECISIVE_MARGIN = 0.05

SIMPLIFIED = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E",
              "T": "C", "S": "C", " ": "C", "": "C", "C": "C",
              "NA": "NA"}


# Display character for a residue DSSP could not assign. mdtraj writes it as
# the two-character "NA", which is why nothing here joins labels into a string
# before substituting this: one two-character label would shift every index
# after it and quietly misalign the whole motif.
UNASSIGNED = "NA"
UNASSIGNED_CHAR = "-"


@dataclass
class MotifSSE:
    """One motif's secondary structure in one frame.

    `labels` is the measurement, one entry per residue. `string` is for display
    and for an SSE identity comparison, and is only safe because every label is
    rendered to exactly one character.
    """

    motif: str
    labels: list[str] = field(default_factory=list)
    helix_fraction: float | None = None
    strand_fraction: float | None = None
    boundaries: list[dict[str, Any]] = field(default_factory=list)

    @property
    def string(self) -> str:
        """One character per residue, for display and for equality tests."""
        return "".join(UNASSIGNED_CHAR if c == UNASSIGNED else (c if c.strip() else "C")
                       for c in self.labels)

    def as_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {"motif": self.motif, "sse": self.string}
        if self.helix_fraction is not None:
            row["helix_fraction"] = round(self.helix_fraction, 3)
        if self.strand_fraction is not None:
            row["strand_fraction"] = round(self.strand_fraction, 3)
        if self.boundaries:
            row["boundaries"] = self.boundaries
        return row


def simplify(labels: list[str]) -> list[str]:
    return [SIMPLIFIED.get(c, "C") for c in labels]


def fractions(labels: list[str]) -> tuple[float | None, float | None]:
    """Helix and strand fraction, over the residues that were assigned.

    Unassigned residues are excluded from the denominator rather than counted as
    coil: counting them as coil would report a motif as unwinding when what
    actually happened is that a residue went missing.
    """
    usable = [c for c in labels if c != UNASSIGNED]
    if not usable:
        return None, None
    helix = sum(1 for c in usable if c in HELIX_CLASSES)
    strand = sum(1 for c in usable if c in STRAND_CLASSES)
    return helix / len(usable), strand / len(usable)


def boundaries(labels: list[str], residue_numbers: list[int]) -> list[dict[str, Any]]:
    """Where each run of one SSE class starts and ends, in residue numbers.

    Reported as residue numbers rather than offsets into the motif so a boundary
    shift can be read directly against a structure, which is what section 5.6
    means by an SSE boundary position.
    """
    runs: list[dict[str, Any]] = []
    for index, letter in enumerate(labels):
        number = residue_numbers[index] if index < len(residue_numbers) else None
        if runs and runs[-1]["sse"] == letter:
            runs[-1]["end"] = number
            runs[-1]["length"] += 1
        else:
            runs.append({"sse": letter, "start": number, "end": number, "length": 1})
    return [r for r in runs if r["sse"] != UNASSIGNED]


def motif_sse(dssp_row, motif: str, indices: list[int],
              residue_numbers: list[int]) -> MotifSSE:
    """One motif's SSE, sliced out of a whole-chain DSSP row.

    `indices` are positions into the DSSP row, which is per protein residue, and
    not residue numbers: mdtraj indexes its DSSP output by its own residue order
    and a trimmed construct's numbering does not start at one.
    """
    labels = [str(dssp_row[i]) if 0 <= i < len(dssp_row) else UNASSIGNED
              for i in indices]
    helix, strand = fractions(labels)
    return MotifSSE(motif=motif, labels=labels, helix_fraction=helix,
                    strand_fraction=strand,
                    boundaries=boundaries(labels, residue_numbers))


# ===========================================================================
#  Agreement with the reference ensembles
# ===========================================================================


def reference_profile(reference_rows: dict[str, list[str]]) -> dict[int, Counter]:
    """Per residue index, how many references show each SSE class.

    Keyed by index into the aligned residue list, so the caller is responsible
    for having aligned the references to the trajectory's residues. That is
    deliberate: silently aligning here would hide a construct mismatch, and a
    mismatch is the one thing that would make every agreement number wrong while
    all of them still looked plausible.
    """
    profile: dict[int, Counter] = {}
    for row in reference_rows.values():
        for index, letter in enumerate(row):
            if letter == UNASSIGNED:
                continue
            profile.setdefault(index, Counter())[letter] += 1
    return profile


def agreement(frame_row: list[str], active: dict[int, Counter],
              inactive: dict[int, Counter]) -> list[dict[str, Any]]:
    """Per residue, the fraction of each reference set sharing this frame's class.

    A residue the frame could not assign, or that no reference resolved, is
    omitted rather than scored zero: zero agreement and no data are different
    claims, and only one of them is evidence of a conformational change.
    """
    rows: list[dict[str, Any]] = []
    for index, letter in enumerate(frame_row):
        if letter == UNASSIGNED:
            continue
        act, inact = active.get(index), inactive.get(index)
        if not act and not inact:
            continue
        row: dict[str, Any] = {"index": index, "sse": letter}
        if act:
            row["active_fraction"] = round(act.get(letter, 0) / sum(act.values()), 3)
            row["n_active"] = sum(act.values())
        if inact:
            row["inactive_fraction"] = round(
                inact.get(letter, 0) / sum(inact.values()), 3)
            row["n_inactive"] = sum(inact.values())
        rows.append(row)
    return rows


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Mean agreement with each reference set, over the residues that have both.

    Restricted to residues present in both sets so the two numbers are comparable:
    averaging each over whatever residues it happened to resolve would make a
    difference between them that came from coverage look like one that came from
    conformation.
    """
    both = [r for r in rows
            if "active_fraction" in r and "inactive_fraction" in r]
    if not both:
        return {"n_residues": 0}
    active = sum(r["active_fraction"] for r in both) / len(both)
    inactive = sum(r["inactive_fraction"] for r in both) / len(both)
    margin = abs(active - inactive)
    return {"n_residues": len(both),
            "mean_active_agreement": round(active, 3),
            "mean_inactive_agreement": round(inactive, 3),
            "margin": round(margin, 3),
            # Secondary structure is far less discriminating than geometry:
            # most residues of a transmembrane helix are helical in both
            # states, so every agreement here is high and the differences are
            # small. H8 came out 0.864 active against 0.888 inactive, and
            # calling that "inactive" on a 0.024 margin would be reading noise.
            "decisive": margin >= DECISIVE_MARGIN,
            "leans": ("active" if active > inactive else
                      "inactive" if inactive > active else "neither")
                     if margin >= DECISIVE_MARGIN else "neither"}


def dssp_rows(traj, simplified: bool = False):
    """`compute_dssp` on every frame, with the chain restricted to protein.

    mdtraj refuses to run DSSP on a topology containing lipids, so the protein
    is sliced out first. Done here rather than at each call site because every
    membrane caller would otherwise have to remember it, and the failure is a
    confusing one: it reports an unknown residue rather than a selection problem.
    """
    import mdtraj as md

    protein = traj.topology.select("protein")
    sliced = traj.atom_slice(protein) if len(protein) < traj.n_atoms else traj
    return md.compute_dssp(sliced, simplified=simplified)


# ===========================================================================
#  Reference secondary structure
# ===========================================================================


def structure_dssp(path, residues: list[int]) -> list[str]:
    """DSSP labels for chosen residue numbers of one structure file.

    Returned in the order `residues` was given, with `NA` where the residue is
    absent, so the result lines up with the caller's list and a gap stays a gap
    rather than shifting everything after it.

    Keyed on the deposited residue number rather than on position, because the
    references do not share a numbering: ABL1's are split between two schemes
    19 residues apart, and a positional read would compare different residues.
    """
    import mdtraj as md

    try:
        loaded = md.load(str(path))
    except Exception:                                           # noqa: BLE001
        return [UNASSIGNED] * len(residues)

    protein = loaded.topology.select("protein")
    if len(protein) == 0:
        return [UNASSIGNED] * len(residues)
    sliced = loaded.atom_slice(protein)
    try:
        labels = md.compute_dssp(sliced, simplified=False)[0]
    except Exception:                                           # noqa: BLE001
        return [UNASSIGNED] * len(residues)

    by_number: dict[int, str] = {}
    for position, residue in enumerate(sliced.topology.residues):
        if position < len(labels):
            by_number[int(residue.resSeq)] = str(labels[position])
    return [by_number.get(number, UNASSIGNED) for number in residues]


def motif_profiles(references: dict[str, list[dict]], *,
                   motif_residues_for, motifs: list[str],
                   cache_dir=None
                   ) -> dict[str, dict[str, dict[int, Counter]]]:
    """Per state, per motif, how many references show each SSE class per residue.

    This is what section 5.6's comparison needs: "per residue, the fraction of
    active references and inactive references sharing the frame's SSE class".
    Built per motif rather than whole-chain because that is the only level at
    which the references are reliably aligned, each one's motif residues coming
    from its own pocket map or residue table.

    Cached on disk, and that is not an optimisation so much as a requirement:
    running DSSP over 116 receptor references takes minutes, and the streaming
    analysis in section 2.3 re-runs every half hour while production is going.
    The labels are a property of the deposited structures and never change, so
    the first pass pays for it and the rest do not. The key covers the residues
    asked for as well as the codes, so a change to a motif's span invalidates
    it rather than silently reusing labels for different residues.
    """
    import hashlib
    import json as _json
    from pathlib import Path as _Path

    raw: dict[str, dict[str, dict[str, list[str]]]] = {}
    key_parts = []
    for state in sorted(references):
        for row in references[state]:
            code = str(row["pdb"]).upper()
            wanted = motif_residues_for(code) or {}
            for motif in motifs:
                residues = wanted.get(motif)
                if residues:
                    key_parts.append(f"{state}:{code}:{motif}:"
                                     f"{residues[0]}-{residues[-1]}:{len(residues)}")

    cache_path = None
    if cache_dir is not None:
        digest = hashlib.sha256("|".join(sorted(key_parts)).encode()).hexdigest()[:20]
        cache_path = _Path(cache_dir) / f"sse_{digest}.json"
        if cache_path.is_file():
            try:
                raw = _json.loads(cache_path.read_text())
            except _json.JSONDecodeError:
                cache_path.unlink(missing_ok=True)
                raw = {}

    if not raw:
        raw = {}
        for state, rows in references.items():
            per_motif: dict[str, dict[str, list[str]]] = {}
            for row in rows:
                code = str(row["pdb"]).upper()
                wanted = motif_residues_for(code) or {}
                for motif in motifs:
                    residues = wanted.get(motif)
                    if not residues:
                        continue
                    per_motif.setdefault(motif, {})[code] = structure_dssp(
                        row["path"], residues)
            raw[state] = per_motif
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(_json.dumps(raw))

    return {state: {motif: reference_profile(by_code)
                    for motif, by_code in per_motif.items()}
            for state, per_motif in raw.items()}


def motif_agreement(frame_labels: list[str], motif: str,
                    profiles: dict[str, dict[str, dict[int, Counter]]]
                    ) -> dict[str, object]:
    """One motif's SSE agreement with each reference set, for one frame."""
    active = (profiles.get("active") or {}).get(motif) or {}
    inactive = (profiles.get("inactive") or {}).get(motif) or {}
    if not active and not inactive:
        return {}
    rows = agreement(frame_labels, active, inactive)
    found = summarise(rows)
    found["motif"] = motif
    return found
