"""Reference motif coordinates, in the form the state score consumes.

`state.py` scores a frame by comparing each motif against every reference, which
needs each reference's motif as coordinates in the same core frame as the frame
being scored. This module produces that:

    {"active": {"2GQG": {"A-loop": array, "hinge": array, "_core": array}}, ...}

Every residue is addressed by its deposited number, resolved through the same
KLIFS pocket map or GPCRdb residue table the trajectory uses. That is the part
worth being careful about. ABL1's references are split between the ABL1a and
ABL1b numbering, 19 residues apart, and a motif compared across that offset
would be compared against the wrong residues while returning a perfectly
plausible RMSD. Nothing here types a residue number for exactly that reason.

A reference missing a residue of a motif is dropped for that motif and kept for
the others, rather than dropped entirely or padded. Padding would invent
coordinates; dropping the whole reference would throw away the motifs it does
resolve, and the ensembles are small enough that each one matters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CORE_KEY = "_core"


@dataclass
class MotifEnsemble:
    """Motif coordinates for one state's references, plus what was dropped."""

    state: str
    by_code: dict[str, dict[str, Any]] = field(default_factory=dict)
    dropped: list[dict[str, Any]] = field(default_factory=list)

    def as_row(self) -> dict[str, Any]:
        return {
            "state": self.state,
            "codes": sorted(self.by_code),
            "motifs_per_code": {code: sorted(k for k in motifs if k != CORE_KEY)
                                for code, motifs in sorted(self.by_code.items())},
            "dropped": list(self.dropped),
        }


def ca_by_residue(path: Path, chain_id: str = "") -> dict[int, tuple]:
    """Calpha positions keyed by deposited residue number.

    Alternative conformations are removed first. Two conformations of the same
    residue both carry a Calpha, and taking whichever came second would make a
    motif's coordinates depend on file order.
    """
    import gemmi

    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()

    found: dict[int, tuple] = {}
    for chain in structure[0]:
        if chain_id and chain.name != chain_id:
            continue
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is None or not info.is_amino_acid():
                continue
            atom = residue.find_atom("CA", "*")
            if atom is not None:
                found.setdefault(residue.seqid.num,
                                 (atom.pos.x, atom.pos.y, atom.pos.z))
    return found


def motif_coordinates(ca: dict[int, tuple], residues: list[int]):
    """Coordinates for a motif, or None if any residue is missing.

    All or nothing per motif. A motif measured on the subset of residues one
    reference happens to resolve is not comparable with the same motif measured
    on a different subset of another, and the RMSD between them would still be
    a number.
    """
    import numpy as np

    if not residues:
        return None
    try:
        return np.array([ca[n] for n in residues], dtype=float)
    except KeyError:
        return None


# How many of the references a residue has to appear in to stay in a motif.
# There is a real trade-off here and it is not avoidable: requiring every
# reference leaves a 38-residue helix with nothing, and requiring none leaves
# the motif defined by whichever structure resolved least. 0.9 keeps the motif
# nearly whole while discarding the few partial depositions that would otherwise
# set its size for everyone.
MOTIF_COVERAGE = 0.9
MOTIF_MINIMUM_RESIDUES = 3


def common_motif_residues(references: dict[str, list[dict[str, Any]]], *,
                          motif_residues_for: callable,
                          coverage: float = MOTIF_COVERAGE,
                          minimum: int = MOTIF_MINIMUM_RESIDUES
                          ) -> tuple[dict[str, list[int]], dict[str, Any]]:
    """Per motif, the residues enough references resolve to be comparable.

    All-or-nothing per motif is the right rule for comparing one motif across
    two structures, and applied to the residue list a caller *wants* it silently
    drops the large motifs. TM6 is 38 residues through the most mobile helix of
    the receptor and no deposition resolves all of them, so every reference
    failed the check and the motif went unscored; the kinase activation loop
    went the same way. Both are the motif that carries the state, which makes
    them the worst two to lose, and nothing in the output said they were gone.

    Intersecting across every reference is no better: on 116 GPCRdb entries it
    leaves every motif empty, because a handful of partial depositions decide
    the answer for all of them.

    So a residue stays if `coverage` of the references resolve it, and the
    references that do not resolve the surviving set are dropped for that motif
    by `build`. The correspondence stays exact residue to residue, which is what
    makes an RMSD between two structures mean anything, and the cost is visible:
    the report says how many residues each motif kept and how many references
    can still be compared on it.
    """
    import numpy as np

    counts: dict[str, dict[int, int]] = {}
    present_in: dict[str, int] = {}
    wanted: dict[str, list[int]] = {}
    seen = 0
    for rows in references.values():
        for row in rows:
            path = Path(row["path"])
            if not path.is_file():
                continue
            ca = ca_by_residue(path, str(row.get("chain") or ""))
            if not ca:
                continue
            seen += 1
            for motif, residues in (motif_residues_for(
                    str(row["pdb"]).upper()) or {}).items():
                wanted.setdefault(motif, list(residues))
                tally = counts.setdefault(motif, {})
                found = [n for n in residues if n in ca]
                # Coverage is taken over the references that have this motif at
                # all, not over every reference. Helix 8 is absent or disordered
                # in most beta2AR depositions, and measured against all of them
                # it kept no residues; measured against the ones that have it,
                # it is perfectly comparable on a smaller ensemble. The ensemble
                # size is reported so a reader knows which it is.
                if len(found) >= minimum:
                    present_in[motif] = present_in.get(motif, 0) + 1
                    for number in found:
                        tally[number] = tally.get(number, 0) + 1

    if not seen:
        return {}, {"n_references": 0}

    kept: dict[str, list[int]] = {}
    report: dict[str, Any] = {"n_references": seen, "coverage": coverage,
                              "motifs": {}, "dropped": []}
    for motif, tally in counts.items():
        available = present_in.get(motif, 0)
        residues = sorted(n for n, count in tally.items()
                          if count >= coverage * available)
        report["motifs"][motif] = {
            "requested": len(wanted.get(motif, [])),
            "kept": len(residues),
            "references_with_motif": available,
        }
        if len(residues) < minimum or available < 2:
            report["dropped"].append(motif)
            continue
        kept[motif] = residues
    return kept, report


def build(references: dict[str, list[dict[str, Any]]], *,
          motif_residues_for: callable,
          core_residues_for: callable) -> dict[str, MotifEnsemble]:
    """Motif coordinates for every reference, grouped by state.

    `motif_residues_for(code)` and `core_residues_for(code)` return that
    reference's own residue numbers, which is where the per-structure resolution
    happens: the caller passes functions that consult the pocket map or residue
    table for that code rather than one shared residue list.

    `references` maps a state to rows each carrying at least `pdb` and `path`.
    """
    out: dict[str, MotifEnsemble] = {}
    for state, rows in references.items():
        ensemble = MotifEnsemble(state=state)
        for row in rows:
            code = str(row["pdb"]).upper()
            path = Path(row["path"])
            if not path.is_file():
                ensemble.dropped.append({"pdb": code, "reason": "structure not on disk"})
                continue
            ca = ca_by_residue(path, str(row.get("chain") or ""))
            if not ca:
                ensemble.dropped.append({"pdb": code, "reason": "no Calpha read"})
                continue

            core = motif_coordinates(ca, core_residues_for(code))
            if core is None:
                ensemble.dropped.append(
                    {"pdb": code, "reason": "core incomplete in this structure"})
                continue

            motifs: dict[str, Any] = {CORE_KEY: core}
            missing = []
            for motif, residues in (motif_residues_for(code) or {}).items():
                coords = motif_coordinates(ca, residues)
                if coords is None:
                    missing.append(motif)
                    continue
                motifs[motif] = coords
            ensemble.by_code[code] = motifs
            if missing:
                ensemble.dropped.append({"pdb": code, "partial": True,
                                         "motifs_missing": sorted(missing)})
        out[state] = ensemble
    return out


def for_scoring(ensembles: dict[str, MotifEnsemble]) -> dict[str, dict[str, Any]]:
    """The plain nested dict `traj.analyse` takes."""
    return {state: ensemble.by_code for state, ensemble in ensembles.items()}


def align_to_common_frame(ensembles: dict[str, MotifEnsemble]) -> dict[str, Any]:
    """Put every reference's coordinates into one shared frame, on its core.

    Each reference comes from its own mmCIF and so sits in its own
    crystallographic frame, which are unrelated to each other. Scoring a frame
    against them without this step compares the frame against one reference's
    orientation and against every other reference's arbitrary one, which
    produced motif RMSDs of 17 to 22 A to the inactive set while the active
    reference that happened to be the alignment target read 1 to 3 A.

    That failure is worth describing because of how it looked: the scores were
    confidently active for every frame, every motif, with a plausible spread, and
    the cause was that half the comparison was nonsense. Nothing in the output
    said so.

    The target is the first reference, in sorted order, whose core is the most
    common length in the set. Sorted order alone is not enough and the
    difference is not cosmetic: a reference whose core is one residue shorter
    than everyone else's fails the shape check below for every other reference,
    and if that reference happens to sort first the whole ensemble is thrown
    away. On CDK2 it did. The target was 1F5Q with a 41-residue core, 28 of the
    other 29 references had 42, and one reference survived out of twenty-nine.

    Choosing by modal length cannot change a set that is already uniform, which
    every ABL1 and beta2AR ensemble here is, so the production rows are
    unaffected and CDK2 goes from one reference to twenty-eight.

    Intersecting the cores instead would be the more thorough fix and it is not
    available: these are residue numbers, and ABL1's references are split
    between the ABL1a and ABL1b numbering nineteen residues apart, so the
    intersection across its nineteen references is three residues. The
    intersection would have to be taken over KLIFS positions and resolved back
    per structure, which is a change to how cores are built rather than to how
    one is chosen.

    The transform is fitted on the core and applied to every motif, which is the
    same rule the per-frame scoring uses.
    """
    from .traj import apply_transform, rmsd, superpose

    import numpy as np

    codes = [(state, code) for state in sorted(ensembles)
             for code in sorted(ensembles[state].by_code)]
    if not codes:
        return {"n_aligned": 0}

    from collections import Counter

    lengths = Counter(
        len(ensembles[state].by_code[code][CORE_KEY]) for state, code in codes)
    modal = lengths.most_common(1)[0][0]
    target_state, target_code = next(
        (state, code) for state, code in codes
        if len(ensembles[state].by_code[code][CORE_KEY]) == modal)
    target = np.asarray(
        ensembles[target_state].by_code[target_code][CORE_KEY], dtype=float)

    aligned, skipped, residuals = 0, [], {}
    for state, code in codes:
        motifs = ensembles[state].by_code[code]
        core = np.asarray(motifs[CORE_KEY], dtype=float)
        if core.shape != target.shape:
            skipped.append(code)
            continue
        transform = superpose(core, target)
        for name in list(motifs):
            motifs[name] = apply_transform(
                np.asarray(motifs[name], dtype=float), *transform)
        residuals[code] = rmsd(motifs[CORE_KEY], target)
        aligned += 1

    for code in skipped:
        for state in ensembles:
            if code in ensembles[state].by_code:
                del ensembles[state].by_code[code]
                ensembles[state].dropped.append(
                    {"pdb": code, "reason": "core length differs from the "
                                            "alignment target"})
    return {"n_aligned": aligned, "target": target_code,
            "n_dropped": len(skipped),
            "mean_core_rmsd_a": round(float(np.mean(list(residuals.values()))), 3)
                                if residuals else None,
            "max_core_rmsd_a": round(float(np.max(list(residuals.values()))), 3)
                               if residuals else None}


def rmsd_to_nearest(ensembles: dict[str, MotifEnsemble], motif: str,
                    state: str) -> dict[str, float]:
    """Per reference, its motif RMSD to the nearest reference of one state.

    The second axis of section 5.6's GPCR state map, which wants the TM3-TM6
    cleft distance against the NPxxY RMSD to inactive. The cleft distance comes
    from the phase 1 manifest; this is the other half, and it is computed from
    the reference coordinates rather than stored, because they are already
    aligned into one frame and recomputing is cheaper than another manifest
    field to keep in step.

    A reference is not compared with itself, so an inactive structure's distance
    to the inactive set is its distance to the *other* inactive structures. The
    alternative is a column of zeros down one state, which would make the map
    look like two lines rather than two clouds.
    """
    import numpy as np

    from .traj import rmsd

    targets = ensembles.get(state)
    if targets is None:
        return {}

    out: dict[str, float] = {}
    for holder in ensembles.values():
        for code, motifs in holder.by_code.items():
            mine = motifs.get(motif)
            if mine is None:
                continue
            distances = [
                rmsd(np.asarray(mine, dtype=float),
                     np.asarray(theirs[motif], dtype=float))
                for other, theirs in targets.by_code.items()
                if other != code and theirs.get(motif) is not None
                and len(theirs[motif]) == len(mine)]
            if distances:
                out[code] = float(min(distances))
    return out
