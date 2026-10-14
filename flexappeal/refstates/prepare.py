"""Everything the reference-state score needs, resolved before the run starts.

Scoring a frame is numpy: superpose it on a core and measure each motif against
a set of reference coordinates. Finding out which residues those are is not. It
needs KLIFS or GPCRdb over the network, the deposited structures behind them,
and gemmi to read those, and a generated bundle is meant to run on a machine
that has an OpenMM environment and nothing else.

So the two halves run in different places. This module runs where the network
and the caches already are, which is wherever the bundle is built, and returns
a few hundred kilobytes of plain JSON: the motif residues of *this* construct,
the superposition core, and the reference motif coordinates to score against.
The bundle carries that, and `analyse.py` does the numpy.

The size is the reason this works at all. The largest of the four worked
examples is the beta2 adrenergic receptor, 100 reference structures across both
states, 107 residues of Calpha coordinates each: 125 kB as float32. The ligand
chemistry a bundle already carries is the same idea and a similar size.

Nothing here is a hand-typed residue number. The motifs are resolved from KLIFS
pocket positions or GPCRdb generic numbers at build time, against the
construct's own deposited numbering, which is the rule the whole build is held
to.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

# What goes into the bundle, and what the results file then carries.
MEMBER = "refstates/ensemble.json"

# A payload this module cannot produce is omitted rather than guessed at, so
# the version is here to tell a reader of an old bundle what shape to expect.
PAYLOAD_VERSION = 1


class PrepareError(RuntimeError):
    """The annotation cannot be set up for this structure.

    Raised with a sentence a user can act on, because the alternative is a
    bundle that builds happily and fails hours later on their machine.
    """


def resolve(structure_path: Path, cfg: dict[str, Any]) -> dict[str, Any]:
    """The reference-state payload for one configuration, or raise.

    `structure_path` is the deposited structure as uploaded, before solvation,
    because every motif here is addressed by its deposited residue number and
    OpenMM renumbers from 1 the moment it solvates.
    """
    from . import ensembles
    from .setup import GPCR_MOTIFS, KINASE_MOTIFS, gpcr_setup, kinase_setup

    family = str(cfg.get("refstates_family") or "auto")
    if family == "auto":
        family = detect_family(structure_path)
        if family == "unknown":
            raise PrepareError(
                "could not tell whether this is a kinase or a receptor from "
                "its structure. Set the protein family explicitly, or turn the "
                "reference-state scoring off: scoring against the wrong "
                "catalogue's motifs is worse than not scoring.")

    name = str(cfg.get("refstates_protein") or "").strip()
    limit = int(cfg.get("refstates_max_references") or 0) or None

    try:
        if family == "kinase":
            setup = kinase_setup(
                structure_path,
                str(cfg.get("refstates_reference_pdb") or "").strip().upper(),
                protein=name or "ABL1", max_references=limit)
            motifs = KINASE_MOTIFS
        else:
            setup = gpcr_setup(structure_path, protein=name or "adrb2_human",
                               max_references=limit)
            motifs = GPCR_MOTIFS
    except Exception as exc:                        # noqa: BLE001
        # Deliberately broad. The catalogues raise their own exception types,
        # gemmi raises its own, and a user who asked for this in a web form
        # needs one sentence rather than a class name. The original is kept in
        # the message because it is usually the useful part.
        raise PrepareError(
            f"the reference-state setup failed for "
            f"{name or 'the default protein for this family'}: {exc}") from exc

    motif_residues, core, built, _own, report = setup[:5]
    if not motif_residues:
        raise PrepareError(
            f"none of the {family} motifs could be resolved against this "
            f"structure's own residue numbering")

    scoring = ensembles.for_scoring(built)
    counts = {state: len(codes) for state, codes in scoring.items()}
    if not any(counts.values()):
        raise PrepareError(
            f"no reference structures were resolved for "
            f"{name or family}, so there would be nothing to score against")

    return {
        "version": PAYLOAD_VERSION,
        "kind": family,
        "protein": name,
        "motifs": list(motifs),
        "motif_residues": {k: list(v) for k, v in motif_residues.items()},
        "core": {"residues": list(core.residues), "source": core.source,
                 "helices": core.helices},
        "reference_motifs": _plain(scoring),
        "ensembles": {state: ensemble.as_row()
                      for state, ensemble in built.items()},
        "report": _plain(report),
        "counts": counts,
    }


def detect_family(structure_path: Path) -> str:
    """kinase, gpcr or unknown, from the structure alone.

    The same detection the CLI uses, so a bundle and a command line agree about
    what a file is.
    """
    from .cli import detect_family as _detect

    return _detect(structure_path)


def _plain(value: Any) -> Any:
    """Whatever this is, as something `json.dumps` will take.

    The setups return numpy arrays inside dictionaries inside dataclasses, and
    the payload has to survive a round trip through a zip member.
    """
    import numpy as np

    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "as_row"):
        return _plain(value.as_row())
    return str(value)


def dumps(payload: dict[str, Any]) -> bytes:
    """The member's bytes. Compact, because it travels in two files."""
    return json.dumps(payload, separators=(",", ":")).encode("utf-8")
