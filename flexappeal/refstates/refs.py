"""Reference structures, prepared down the same code path a trajectory frame takes.

Section 4.4 of the spec is specific about this and it is the design decision the
positive controls rest on: a reference is a one-frame trajectory. The same
annotator resolves its motifs, the same classifier reads its DFG and αC states,
the same metrics come out. If references went down a shortcut, P1 to P3 would be
measuring the shortcut rather than the code that annotates a run.

So what this module does is narrow: fetch the coordinates, clean them to one
chain, hand them to the vendored `sse_comparison` annotator, and read the labels
back. The labels are then checked against what KLIFS and GPCRdb say themselves,
which is how the controls are evaluated.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RCSB_CIF = "https://files.rcsb.org/download/{pdb}.cif.gz"

# Section 9's agreement thresholds, which are what the controls pass or fail on.
P1_MIN_AGREEMENT = 0.90        # our DFG label against KLIFS's
P2_MIN_AGREEMENT = 0.85        # our alphaC label against KLIFS's


@dataclass
class PreparedReference:
    """One reference structure, with what both we and the catalogue think of it."""

    pdb: str
    chain: str
    path: Path
    sequence: str = ""
    motifs: list = field(default_factory=list)
    our_dfg: str | None = None
    our_alphac: str | None = None
    their_dfg: str | None = None
    their_alphac: str | None = None
    status: str = "ok"
    note: str = ""
    # The geometry the DFG call was read from, carried so a reader can check the
    # label and so the old absolute-threshold answer stays recoverable.
    dfg_d_alphac: float | None = None
    dfg_d_lys: float | None = None
    dfg_absolute_state: str | None = None

    @property
    def dfg_agrees(self) -> bool | None:
        """None where the catalogue has no opinion, which is not a disagreement."""
        if not self.their_dfg or self.their_dfg == "na" or not self.our_dfg:
            return None
        return self.our_dfg == self.their_dfg

    @property
    def alphac_agrees(self) -> bool | None:
        if not self.their_alphac or self.their_alphac == "na" or not self.our_alphac:
            return None
        return self.our_alphac == self.their_alphac


def fetch_mmcif(pdb: str, cache_dir: Path, *, refresh: bool = False) -> Path:
    """One structure from the RCSB, cached gzipped.

    Kept compressed on disk: an mmCIF is mostly text and a reference ensemble is
    a hundred of them, so this is the difference between tens and hundreds of
    megabytes for data that is downloaded once and read a few times.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb.lower()}.cif.gz"
    if path.is_file() and not refresh and path.stat().st_size > 0:
        return path

    import requests
    from tenacity import (retry, retry_if_exception_type, stop_after_attempt,
                          wait_exponential)

    @retry(stop=stop_after_attempt(4),
           wait=wait_exponential(multiplier=1, min=1, max=20),
           retry=retry_if_exception_type((requests.ConnectionError, requests.Timeout)),
           reraise=True)
    def _get() -> bytes:
        response = requests.get(RCSB_CIF.format(pdb=pdb.upper()), timeout=60)
        if response.status_code >= 500:
            raise requests.ConnectionError(f"RCSB returned {response.status_code}")
        response.raise_for_status()
        return response.content

    path.write_bytes(_get())
    return path


def _uncompressed(path: Path, work_dir: Path) -> Path:
    """gemmi reads a plain file; the cache holds a gzipped one."""
    work_dir.mkdir(parents=True, exist_ok=True)
    out = work_dir / path.name.replace(".gz", "")
    if not out.is_file():
        with gzip.open(path, "rb") as src:
            out.write_bytes(src.read())
    return out


def prepare_kinase_reference(entry: dict[str, Any], *, cache_dir: Path,
                             work_dir: Path, annotator,
                             kinase_name: str = "ABL1") -> PreparedReference:
    """Fetch, clean and annotate one KLIFS entry, then read our own labels off it.

    `entry` is a row from KLIFS's `structures_list`, so it carries the
    catalogue's own `DFG` and `aC_helix` calls. Those are recorded beside ours
    rather than used: the whole point of P1 and P2 is that the two are compared.
    """
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vendor" / "boltzmaker"))
    from sse_comparison import metrics, structures
    from sse_comparison.alignment import build_comparison_frame

    pdb = str(entry["pdb"]).upper()
    chain = str(entry.get("chain") or "A")
    reference = PreparedReference(
        pdb=pdb, chain=chain, path=Path(),
        their_dfg=str(entry.get("DFG") or "").lower() or None,
        their_alphac=str(entry.get("aC_helix") or "").lower() or None)

    try:
        reference.path = fetch_mmcif(pdb, cache_dir)
        plain = _uncompressed(reference.path, work_dir)
        loaded = structures.load_structure_for_comparison(
            plain, chain, expected_sequence="", keep_ligand_codes=frozenset())
        reference.sequence = loaded.sequence
        # `name_hint` is not optional in practice: `KLIFSAnnotator.annotate`
        # returns an empty list without one, because it is what the pocket
        # lookup is keyed on. Passing the chain and the path alone gets no
        # motifs and no error, which is how this first ran.
        reference.motifs = annotator.annotate(
            loaded.sequence, pdb_id=pdb, structure_path=plain, chain_id=chain,
            name_hint=kinase_name)
        if not reference.motifs:
            reference.status = "failed:no_motifs"
            reference.note = "the annotator returned no motifs for this chain"
            return reference

        # A reference compared against itself: the mapping is the identity and
        # the classifier reads the geometry, which is exactly what a trajectory
        # frame does against its own topology.
        frame = build_comparison_frame(
            loaded, loaded, loaded.sequence, reference.motifs, "kinase")
        by_name = {m.name: m for m in reference.motifs}
        dfg = by_name.get("DFG")
        lys = by_name.get("catalytic_Lys")
        glu = by_name.get("alphaC_Glu")
        if dfg and lys and glu:
            # KLIFS positions 81, 82, 83 are Asp, Phe, Gly, so the ring the
            # classifier wants is the middle one.
            phe = dfg.residues[1] if len(dfg.residues) > 1 else dfg.residues[0]
            # The classifiers return (apo, holo, changed) because they are built
            # to compare two structures. Here both halves are this reference, so
            # the two labels are the same one and `changed` is always False.
            # The DFG call comes from `flexappeal.refstates.kinase`, not from the
            # vendored classifier: its absolute 11.0 A threshold cuts through a
            # cluster of real DFG-out structures and mislabels six of ABL1's own
            # references as "other" (D010). The distances behind the call are the
            # vendored measurement, unchanged.
            from .kinase import classify_dfg

            call = classify_dfg(frame, by_name)
            reference.our_dfg = call.state
            reference.dfg_d_alphac = call.d_alphac_a
            reference.dfg_d_lys = call.d_lys_a
            reference.dfg_absolute_state = call.absolute_threshold_state
            apo_state, holo_state, _ = metrics.classify_alphac_state(
                frame, glu.residues[0], lys.residues[0])
            reference.our_alphac = apo_state or holo_state
        else:
            reference.status = "failed:missing_motifs"
            reference.note = f"needed DFG, catalytic_Lys and alphaC_Glu; got {sorted(by_name)}"
    except Exception as exc:                                  # noqa: BLE001
        # One reference never blocks the ensemble (spec 0.1 rule 3).
        reference.status = f"failed:{type(exc).__name__}"
        reference.note = str(exc)[:200]
    return reference


def agreement(references: list[PreparedReference], which: str) -> dict[str, Any]:
    """Positive control P1 or P2: our label against the catalogue's.

    Entries the catalogue labels `na` are excluded rather than counted as
    disagreements, which is what section 9 asks for: no opinion is not a wrong
    opinion.
    """
    attr = "dfg_agrees" if which == "dfg" else "alphac_agrees"
    verdicts = [getattr(r, attr) for r in references if r.status == "ok"]
    compared = [v for v in verdicts if v is not None]
    agreed = sum(1 for v in compared if v)
    threshold = P1_MIN_AGREEMENT if which == "dfg" else P2_MIN_AGREEMENT
    return {
        "control": "P1" if which == "dfg" else "P2",
        "label": which,
        "compared": len(compared),
        "agreed": agreed,
        "fraction": (agreed / len(compared)) if compared else None,
        "threshold": threshold,
        "passes": bool(compared) and (agreed / len(compared)) >= threshold,
        "no_opinion": len(verdicts) - len(compared),
        "failed_to_prepare": sum(1 for r in references if r.status != "ok"),
    }


def write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    """One JSON object per line, so a restart can skip what is already done."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, default=str) + "\n")
