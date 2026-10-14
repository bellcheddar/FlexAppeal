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
        "reference_structures": _reference_structures(family, name, setup,
                                                      limit),
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


def _resolution(entry: dict) -> dict[str, Any]:
    """The resolution as a number, or not at all.

    KLIFS returns it as a string and GPCRdb as a number, and the table that
    prints it formats with `.2f`. Omitted rather than passed through where it
    cannot be read, which is this build's rule for a value it cannot produce:
    a reference with no resolution is a reference with no resolution, not one
    at zero angstrom.
    """
    raw = entry.get("resolution")
    try:
        return {"resolution": float(raw)}
    except (TypeError, ValueError):
        return {}


def _reference_structures(family: str, name: str, setup,
                          limit: int | None) -> list[dict[str, Any]]:
    """Each reference measured on the plane its published criterion is drawn on.

    This is the state map: the catalogue's own structures as a cloud, with the
    run's path across it. The cloud is a property of the references rather than
    of the trajectory, so it is measured here, once, and carried.

    It is the one part of this payload that needs the vendored annotator, and
    only for a kinase: the two DFG anchor distances are measured by
    `sse_comparison`'s own atom lookups, which carry its handling of a missing
    atom and of the ring falling back through CZ, CE1 and CD1. Wrapping rather
    than reimplementing is the rule, and it costs nothing here because this
    runs where that code is.

    Returns an empty list rather than raising. The other panels are the
    measurement; this one is the context for it, and a run that cannot have it
    is better off without it than refused.
    """
    import tempfile

    from . import refs
    from .setup import GPCRDB_CACHE, KLIFS_CACHE

    try:
        if family == "kinase":
            return _kinase_references(name or "ABL1", KLIFS_CACHE,
                                      Path(tempfile.mkdtemp()), refs, limit)
        return _gpcr_references(name or "adrb2_human", GPCRDB_CACHE, setup,
                                refs, limit)
    except Exception:                               # noqa: BLE001
        return []


def _kinase_references(name: str, cache_dir: Path, work_dir: Path,
                       refs, limit: int | None) -> list[dict[str, Any]]:
    """The DFG plane, from the catalogue's own entries.

    Capped the same way the scoring ensemble is. CDK2 has 337 active entries
    and each one costs a structure read and an annotation, so an uncapped cloud
    would turn a form submission into a quarter of an hour. The cap is per
    state, so the default of 40 draws up to 80 points, which is a cloud.
    """
    import sys

    from . import Catalogue, kinase_id, kinase_references

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vendor"
                           / "boltzmaker"))
    from sse_comparison.annotators.klifs import KLIFSAnnotator

    catalogue = Catalogue(cache_dir)
    ensembles_by_state = kinase_references(catalogue,
                                           kinase_id(catalogue, name=name))
    annotator = KLIFSAnnotator()
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for state, ensemble in ensembles_by_state.items():
        entries = list(getattr(ensemble, "entries", []))
        for entry in (entries[:limit] if limit else entries):
            code = str(entry["pdb"]).upper()
            if code in seen:
                continue
            seen.add(code)
            prepared = refs.prepare_kinase_reference(
                entry, cache_dir=cache_dir, work_dir=work_dir,
                annotator=annotator, kinase_name=name)
            if prepared.dfg_d_alphac is None or prepared.dfg_d_lys is None:
                continue
            rows.append({
                "pdb": prepared.pdb, "chain": prepared.chain,
                "source": "KLIFS", "source_state": state,
                "our_state": prepared.our_dfg or "",
                **_resolution(entry),
                "metrics": {"dfg_d_alphac": prepared.dfg_d_alphac,
                            "dfg_d_lys": prepared.dfg_d_lys,
                            "our_dfg": prepared.our_dfg,
                            "their_dfg": prepared.their_dfg},
            })
    return rows


def _gpcr_references(name: str, cache_dir: Path, setup,
                     refs, limit: int | None) -> list[dict[str, Any]]:
    """The TM3-TM6 and NPxxY plane, from the catalogue's own entries.

    No vendored annotator here: the receptor's axes are a distance between two
    generic-numbered Calphas and a motif RMSD this build already computed while
    assembling the ensemble.
    """
    import gemmi

    from . import Catalogue, receptor_references, receptor_residues
    from .gpcr import build_numbering, tm3_tm6_distance

    npxxy = ((setup[4] or {}).get("map_axis") or {}).get("values") or {}
    catalogue = Catalogue(cache_dir)
    # `receptor_residues` returns the catalogue's rows; the lookup the distance
    # needs is built from them.
    numbering = build_numbering(receptor_residues(catalogue, name))
    rows: list[dict[str, Any]] = []
    for state, ensemble in receptor_references(catalogue, name).items():
        entries = list(getattr(ensemble, "entries", []))
        for entry in (entries[:limit] if limit else entries):
            # GPCRdb calls it `pdb_code` where KLIFS calls it `pdb`, which is
            # why this reads one key and the kinase path above reads another.
            code = str(entry.get("pdb_code") or entry.get("pdb") or "").upper()
            chain = str(entry.get("preferred_chain") or "")[:1]
            if not code:
                continue
            try:
                structure = gemmi.read_structure(
                    str(refs.fetch_mmcif(code, cache_dir)))
                structure.setup_entities()
                structure.remove_alternative_conformations()
            except Exception:                       # noqa: BLE001
                continue
            target = next((c for c in structure[0] if c.name == chain),
                          structure[0][0] if len(structure[0]) else None)
            distance = (tm3_tm6_distance(target, numbering)
                        if target is not None else None)
            if distance is None or code not in npxxy:
                continue
            rows.append({
                "pdb": code, "chain": chain, "source": "GPCRdb",
                "source_state": state, "our_state": "",
                **_resolution(entry),
                "metrics": {"tm3_tm6": distance,
                            "npxxy_to_inactive": npxxy[code]},
            })
    return rows


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
