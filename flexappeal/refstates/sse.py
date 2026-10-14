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

# How many discriminating residues a motif needs before its lean is called.
# One or two is not a measurement: TM5 has two and gave the identical answer
# for a structure built from the active reference and one built from the
# inactive, while TM7 has seven and told them apart correctly in both
# directions. The line goes between those.
MIN_DISCRIMINATING = 3


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


def discriminating_residues(active: dict[int, Counter],
                            inactive: dict[int, Counter]) -> set[int]:
    """Residues where the two reference sets disagree about the dominant class.

    Everywhere else carries no state information, however large the agreement
    difference looks, and this is not a theoretical point. Measured across all
    residues, TM6 gave **0.934 against 0.868 for two different structures**, one
    built from the active reference and one from the inactive: bit-identical,
    because a transmembrane helix is helical in both states so its SSE string is
    the same and the agreement is the same. The lean was reporting that the
    active reference set is more uniformly helical than the inactive one, which
    is a fact about the catalogue rather than about the frame.

    Restricting to residues where the consensus differs is what makes the
    comparison a measurement of the frame. A motif with no such residues cannot
    carry secondary-structure information and says so.
    """
    found = set()
    for index in set(active) & set(inactive):
        one = active[index].most_common(1)
        two = inactive[index].most_common(1)
        if one and two and one[0][0] != two[0][0]:
            found.add(index)
    return found


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


def profiles_from_json(obj) -> dict:
    """Reference SSE profiles, back from a round trip through JSON.

    `motif_profiles` returns {state: {motif: {residue index: Counter}}}, and a
    results file or a bundle member can carry none of that shape: a Counter
    becomes a plain object and an integer key becomes a string. Rebuilt here
    rather than at the call site, because the failure is quiet in exactly the
    wrong way -- `most_common` is missing and the ribbon simply does not draw,
    with the scores around it all correct.
    """
    from collections import Counter as _Counter

    if not obj:
        return {}
    out: dict[str, dict[str, dict[int, _Counter]]] = {}
    for state, motifs in obj.items():
        if not isinstance(motifs, dict):
            continue
        rebuilt: dict[str, dict[int, _Counter]] = {}
        for motif, residues in motifs.items():
            if not isinstance(residues, dict):
                continue
            rebuilt[motif] = {
                int(index): (counts if isinstance(counts, _Counter)
                             else _Counter(counts))
                for index, counts in residues.items()
            }
        out[state] = rebuilt
    return out


def motif_agreement(frame_labels: list[str], motif: str,
                    profiles: dict[str, dict[str, dict[int, Counter]]]
                    ) -> dict[str, object]:
    """One motif's SSE agreement with each reference set, for one frame.

    Scored only over the residues where the reference sets disagree about the
    dominant class, because the rest cannot distinguish the states and including
    them produced the same answer for two different structures.
    """
    active = (profiles.get("active") or {}).get(motif) or {}
    inactive = (profiles.get("inactive") or {}).get(motif) or {}
    if not active and not inactive:
        return {}

    informative = discriminating_residues(active, inactive)
    rows = [r for r in agreement(frame_labels, active, inactive)
            if r["index"] in informative]
    found = summarise(rows)
    found["motif"] = motif
    found["n_discriminating"] = len(informative)
    found["n_compared"] = sum(1 for _ in set(active) & set(inactive))
    if not informative:
        found["leans"] = "neither"
        found["decisive"] = False
        found["note"] = ("the reference sets agree on this motif's secondary "
                         "structure at every residue, so it carries no state "
                         "information")
    elif len(informative) < MIN_DISCRIMINATING:
        found["leans"] = "neither"
        found["decisive"] = False
        found["note"] = (f"only {len(informative)} residue(s) distinguish the "
                         f"reference sets here, too few to read a state from")
    return found


# ---------------------------------------------------------------------------
#  The run's own secondary structure
# ---------------------------------------------------------------------------
#
# Moved here from `tools/analyse_run.py` so a generated bundle can draw the
# DSSP ribbon too. It needs mdtraj and the reference profiles the build-time
# setup already resolved, and nothing else: no catalogue, no network.


def secondary_structure(loaded, motif_residues: dict[str, list[int]],
                        profiles: dict | None = None) -> dict:
    """Per-frame DSSP, and each motif's helix and strand content (section 5.6).

    Stored as the per-motif summary plus the whole-chain per-frame strings. The
    strings are what the DSSP ribbon is drawn from and they are the bulk of it,
    so they go in the run's JSON rather than into SQLite: section 7 puts
    per-residue secondary structure outside the database for exactly this
    reason, being a residues-by-frames matrix rather than a set of rows.

    `simplified=False` so the eight DSSP classes survive. The three-class form
    is derived where it is wanted; going the other way is impossible, and the
    difference between a 3-10 helix and an alpha helix is the kind of local
    change this is looking for.
    """
    try:
        rows = dssp_rows(loaded, simplified=False)
    except Exception as exc:                                    # noqa: BLE001
        log(f"DSSP unavailable: {str(exc)[:80]}")
        return {}
    if rows is None or len(rows) == 0:
        return {}

    protein = [r for r in loaded.topology.residues if r.is_protein]
    numbers = [int(r.resSeq) for r in protein]
    index_of = {number: position for position, number in enumerate(numbers)}

    motifs: dict[str, dict] = {}
    for name, residues in motif_residues.items():
        indices = [index_of[n] for n in residues if n in index_of]
        if not indices:
            continue
        helix, strand = [], []
        for frame in range(len(rows)):
            found = motif_sse(rows[frame], name, indices,
                                  [n for n in residues if n in index_of])
            if found.helix_fraction is not None:
                helix.append(found.helix_fraction)
                strand.append(found.strand_fraction)
        if not helix:
            continue
        motifs[name] = {
            "n_residues": len(indices),
            "helix_fraction_mean": sum(helix) / len(helix),
            "helix_fraction_first": helix[0],
            "helix_fraction_last": helix[-1],
            "strand_fraction_mean": sum(strand) / len(strand),
        }

        # Section 5.6: per residue, the fraction of each reference set sharing
        # this frame's SSE class. Taken at the last frame, which is the one a
        # reader is asking about; the per-frame series is a motif by residue by
        # frame cube and belongs nowhere near a page.
        if profiles:
            found = motif_agreement(
                motif_sse(rows[-1], name, indices,
                              [n for n in residues if n in index_of]).labels,
                name, profiles)
            if found.get("n_residues"):
                motifs[name]["reference_agreement"] = found

    return {
        "n_frames": int(len(rows)),
        "n_residues": len(numbers),
        "residues": numbers,
        # One character per residue per frame, which is only safe because every
        # DSSP label renders to exactly one character; mdtraj writes an
        # unassigned residue as the two-character "NA".
        "frames": ["".join(UNASSIGNED_CHAR if c == UNASSIGNED
                           else (c if str(c).strip() else "C")
                           for c in row) for row in rows],
        "motifs": motifs,
    }
