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


def common_core(core_by_code: dict[str, list[int]]) -> list[int]:
    """Residues present in every reference's core, as section 5.5 asks.

    The core has to be the same length everywhere or the superposition cannot be
    applied across the set, so it is the intersection and not the union. Taken on
    generic or KLIFS positions by the caller, which is what makes an intersection
    meaningful across structures with different numbering.
    """
    if not core_by_code:
        return []
    shared: set[int] | None = None
    for residues in core_by_code.values():
        current = set(residues)
        shared = current if shared is None else (shared & current)
    return sorted(shared or [])


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

    The target is the first reference in sorted order, so the choice is
    reproducible. The transform is fitted on the core and applied to every motif,
    which is the same rule the per-frame scoring uses.
    """
    from .traj import apply_transform, rmsd, superpose

    import numpy as np

    codes = [(state, code) for state in sorted(ensembles)
             for code in sorted(ensembles[state].by_code)]
    if not codes:
        return {"n_aligned": 0}

    target_state, target_code = codes[0]
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
