#!/usr/bin/env python3
"""Rebuild the reference manifests and the controls that come from them.

Spec phase 1. The manifests in `data/manifests/` feed the reference browser,
the state maps and P3's fitted threshold, and until this script existed nothing
in the repository could regenerate them: they were artefacts of one session,
with the code that produced them reachable only by reading it. That is a
reproducibility gap rather than a tidiness one, because a committed number
nobody can re-derive is a number nobody can check.

    python tools/build_references.py            # rebuild both, report the controls
    python tools/build_references.py --verify   # rebuild nothing, compare instead

`--verify` is the useful mode day to day: it recomputes everything and diffs it
against what is on disk, so a change in a catalogue or in our own code shows up
as a disagreement rather than as a silent drift.

Note that the manifests are gitignored, so they are not committed and a fresh
clone has none. That is why this script matters: without it there is no way to
produce them on another machine at all.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

KLIFS_CACHE = ROOT / "data" / "cache" / "klifs"
GPCRDB_CACHE = ROOT / "data" / "cache" / "gpcrdb"
MANIFESTS = ROOT / "data" / "manifests"

# Section 9's pass marks, so this script reports against them rather than
# leaving a reader to compare the numbers by hand.
P1_MIN = 0.90
P2_MIN = 0.85
P3_MIN = 0.90


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def build_kinase(work: Path, all_chains: bool = False) -> list[dict]:
    """Every ABL1 reference, prepared and labelled by our own code.

    One row per PDB entry, not per chain. KLIFS lists a row per chain and per
    alternate conformation, so ABL1 comes back as 61 rows over 20 depositions:
    2G1T alone appears four times. Those copies sit in one asymmetric unit of
    one crystal and are not independent measurements of whether our DFG call
    agrees with KLIFS's, so counting them would weight a four-chain crystal
    four times and report an n that overstates what was actually checked.

    `all_chains` keeps every chain instead. That is a sensitivity check rather
    than the manifest: it answers whether our labels are stable within a
    crystal, which is a different question from whether they agree with KLIFS.
    """
    import sys as _sys

    from flexappeal.refstates import (Catalogue, kinase_id, kinase_references,
                                      refs)

    _sys.path.insert(0, str(ROOT / "vendor" / "boltzmaker"))
    from sse_comparison.annotators.klifs import KLIFSAnnotator

    catalogue = Catalogue(KLIFS_CACHE)
    kid = kinase_id(catalogue)
    ensembles = kinase_references(catalogue, kid)
    annotator = KLIFSAnnotator()

    rows: list[dict] = []
    seen: set[tuple[str, ...]] = set()
    for which in ("active", "inactive", "src_like"):
        ensemble = ensembles.get(which)
        if ensemble is None:
            continue
        for entry in ensemble.entries:
            code = str(entry["pdb"]).upper()
            chain = str(entry.get("chain") or "A")
            key = (code, chain) if all_chains else (code,)
            if key in seen:
                continue
            seen.add(key)
            prepared = refs.prepare_kinase_reference(
                entry, cache_dir=KLIFS_CACHE, work_dir=work,
                annotator=annotator)
            rows.append({
                "pdb": prepared.pdb, "chain": prepared.chain,
                "status": prepared.status,
                "our_dfg": prepared.our_dfg, "their_dfg": prepared.their_dfg,
                "dfg_d_alphac": prepared.dfg_d_alphac,
                "dfg_d_lys": prepared.dfg_d_lys,
                "dfg_absolute_state": prepared.dfg_absolute_state,
                "our_alphac": prepared.our_alphac,
                "their_alphac": prepared.their_alphac,
                "dfg_agrees": prepared.dfg_agrees,
                "alphac_agrees": prepared.alphac_agrees,
                "note": prepared.note,
            })
            log(f"  {prepared.pdb} {prepared.chain}: {prepared.status}, "
                f"DFG {prepared.our_dfg} against {prepared.their_dfg}")
    return rows


def build_receptor() -> tuple[list[dict], dict]:
    """Every beta2AR reference with its cleft distance, and P3's threshold."""
    import gemmi

    from flexappeal.refstates import (Catalogue, receptor_references,
                                      receptor_residues, refs)
    from flexappeal.refstates.gpcr import (build_numbering, fit_threshold,
                                           tm3_tm6_distance)

    catalogue = Catalogue(GPCRDB_CACHE)
    numbering = build_numbering(receptor_residues(catalogue))
    ensembles = receptor_references(catalogue)

    rows: list[dict] = []
    by_state: dict[str, list[float]] = {"active": [], "inactive": []}
    for which in ("active", "inactive", "intermediate"):
        ensemble = ensembles.get(which)
        if ensemble is None:
            continue
        for entry in ensemble.entries:
            code = str(entry.get("pdb_code") or "").upper()
            if not code:
                continue
            chain = str(entry.get("preferred_chain") or "")[:1]
            try:
                path = refs.fetch_mmcif(code, GPCRDB_CACHE)
                structure = gemmi.read_structure(str(path))
                structure.setup_entities()
                structure.remove_alternative_conformations()
            except Exception as exc:                            # noqa: BLE001
                log(f"  {code}: not read ({str(exc)[:50]})")
                continue

            target = next((c for c in structure[0] if c.name == chain),
                          structure[0][0] if len(structure[0]) else None)
            distance = (tm3_tm6_distance(target, numbering)
                        if target is not None else None)
            rows.append({"pdb": code, "chain": chain,
                         "state": which, "tm3_tm6": distance})
            if distance is not None and which in by_state:
                by_state[which].append(distance)

    control = fit_threshold(by_state["active"], by_state["inactive"])
    return rows, control


def rspine_separation() -> dict:
    """Re-derive the crystal separation D015 fitted the R-spine threshold to.

    `R_SPINE_CONTACT_A` is 7.7 angstrom where the kinase literature uses 6.0,
    and the justification was a measurement: the two crystal clusters leave an
    empty region between them and 6.0 sits at the very top of the active one,
    with no room for a side chain at 300 K. That measurement was quoted in
    DECISIONS and nothing could reproduce it, which is the same gap this script
    exists to close for the manifests.

    The state is `broken` when any consecutive pair exceeds the threshold, so
    the quantity that decides the call is each reference's *largest* pair, and
    that is what the gap is reported on.
    """
    import gemmi

    from flexappeal.refstates import (Catalogue, kinase_id, kinase_references,
                                      pocket, refs)
    from flexappeal.refstates.pocket import (R_SPINE_CONTACT_A,
                                             R_SPINE_CONTACT_CRYSTAL_A,
                                             r_spine, r_spine_state)

    catalogue = Catalogue(KLIFS_CACHE)
    ensembles = kinase_references(catalogue, kinase_id(catalogue))

    widest: dict[str, list[float]] = {}
    states: dict[str, list[str]] = {}
    seen: set[str] = set()
    for which in ("active", "inactive", "src_like"):
        ensemble = ensembles.get(which)
        if ensemble is None:
            continue
        for entry in ensemble.entries:
            code = str(entry["pdb"]).upper()
            if code in seen:
                continue
            seen.add(code)
            chain_id = str(entry.get("chain") or "A")
            try:
                structure = gemmi.read_structure(
                    str(refs.fetch_mmcif(code, KLIFS_CACHE)))
                structure.setup_entities()
                structure.remove_alternative_conformations()
                chain = next((c for c in structure[0] if c.name == chain_id),
                             None)
                own = pocket.pocket_map(catalogue, entry.get("structure_ID"))
                state = r_spine_state(chain, r_spine(own, chain))
            except Exception as exc:                            # noqa: BLE001
                log(f"    {code} {chain_id}: not measured "
                    f"({type(exc).__name__})")
                continue
            pairs = [p["distance_a"] for p in state["pairs"]
                     if "distance_a" in p]
            if not pairs:
                continue
            widest.setdefault(which, []).append(max(pairs))
            states.setdefault(which, []).append(state["state"])

    out = spine_verdict(widest, R_SPINE_CONTACT_A, R_SPINE_CONTACT_CRYSTAL_A)
    for which, values in widest.items():
        out[which] = {"n": len(values), "min": round(min(values), 2),
                      "max": round(max(values), 2),
                      "states": sorted(set(states[which]))}
    return out


def spine_verdict(widest: dict[str, list[float]], threshold: float,
                  published: float) -> dict:
    """The gap arithmetic, separated so it can be tested without the network.

    Kept apart from the measurement above because this is the part that decides
    whether the run passes, and a verdict that silently inverted would read as a
    clean result. The measurement needs 20 structures and a catalogue; this
    needs two lists of numbers.

    The empty region runs from the top of the active cluster to the bottom of
    the inactive one. Any threshold inside it classifies every crystal
    correctly, which is what makes 7.7 a position in a gap rather than a tuned
    number.
    """
    out: dict = {"threshold_a": threshold, "published_a": published}
    active, inactive = widest.get("active"), widest.get("inactive")
    if not active or not inactive:
        return out
    low, high = max(active), min(inactive)
    out["gap_a"] = [round(low, 2), round(high, 2)]
    out["gap_width_a"] = round(high - low, 2)
    out["threshold_in_gap"] = low < threshold < high
    out["published_in_gap"] = low < published < high
    out["published_margin_a"] = round(published - low, 2)
    # `assembled` needs every pair under the threshold, so an active reference
    # is correct when its widest pair is under and an inactive one when its
    # widest pair is not.
    out["correct"] = (all(v < threshold for v in active)
                      and all(v >= threshold for v in inactive))
    out["n_classified"] = len(active) + len(inactive)
    return out


def _load(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def compare(fresh: list[dict], committed: list[dict], key: str = "pdb") -> dict:
    """What changed between a rebuild and the manifest on disk.

    Compared per record rather than as whole files, because the order a
    catalogue returns its rows in is not stable and a diff on the file would
    report a reordering as a change.
    """
    def index(rows):
        out = {}
        for row in rows:
            out.setdefault((row.get(key), row.get("chain")), row)
        return out

    one, two = index(fresh), index(committed)
    added = sorted(k for k in one if k not in two)
    removed = sorted(k for k in two if k not in one)
    changed = []
    for shared in sorted(set(one) & set(two)):
        a, b = one[shared], two[shared]
        fields = sorted(set(a) | set(b))
        diffs = {}
        for field in fields:
            first, second = a.get(field), b.get(field)
            if isinstance(first, float) and isinstance(second, float):
                if abs(first - second) > 0.01:
                    diffs[field] = (second, first)
            elif first != second:
                diffs[field] = (second, first)
        if diffs:
            changed.append({"record": shared, "fields": diffs})
    return {"added": added, "removed": removed, "changed": changed,
            "n_fresh": len(one), "n_committed": len(two)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verify", action="store_true",
                        help="rebuild in memory and diff against the committed "
                             "manifests instead of overwriting them")
    parser.add_argument("--only", choices=("kinase", "receptor"), default=None)
    parser.add_argument("--all-chains", action="store_true",
                        help="keep every KLIFS chain rather than one row per "
                             "PDB entry; a sensitivity check, and it will not "
                             "write the manifest")
    args = parser.parse_args()
    if args.all_chains and not args.verify:
        # The committed manifest is per-entry by construction. Writing the
        # per-chain set over it would change what P1 and P2 are counted from
        # without saying so, which is exactly the silent drift this script
        # exists to catch.
        parser.error("--all-chains is a check, not a build; pass --verify too")

    import tempfile

    from flexappeal.refstates import refs

    MANIFESTS.mkdir(parents=True, exist_ok=True)
    problems = 0

    if args.only != "receptor":
        log("rebuilding the ABL1 references")
        with tempfile.TemporaryDirectory() as work:
            kinase = build_kinase(Path(work), all_chains=args.all_chains)
        usable = [r for r in kinase if r.get("status") == "ok"]
        p1 = refs.agreement([_Row(r) for r in usable], "dfg")
        p2 = refs.agreement([_Row(r) for r in usable], "alphac")
        for control in (p1, p2):
            log(f"  {control['control']} {control['label']}: "
                f"{control['agreed']} of {control['compared']} "
                f"({(control.get('fraction') or 0):.0%}), needs "
                f"{control['threshold']:.0%} -> "
                f"{'pass' if control['passes'] else 'FAIL'}"
                f"  ({control['no_opinion']} with no catalogue opinion, "
                f"{control['failed_to_prepare']} not prepared)")
        if not (p1["passes"] and p2["passes"]):
            problems += 1
            log("  one of P1 or P2 is below its pass mark")

        # D015's threshold, re-derived rather than restated.
        spine = rspine_separation()
        for which in ("active", "inactive", "src_like"):
            row = spine.get(which)
            if row:
                log(f"  R-spine {which}: n={row['n']}, deciding pair "
                    f"{row['min']} to {row['max']} A, {'/'.join(row['states'])}")
        if "gap_a" in spine:
            low, high = spine["gap_a"]
            where = "inside" if spine["threshold_in_gap"] else "OUTSIDE"
            verdict = ("every reference correct"
                       if spine["correct"] else "MISCLASSIFIES a reference")
            log(f"  R-spine empty region {low} to {high} A "
                f"({spine['gap_width_a']} A wide); {spine['threshold_a']} A "
                f"sits {where} it, {verdict} "
                f"({spine['n_classified']} crystals)")
            log(f"  the published {spine['published_a']} A clears the active "
                f"cluster by only {spine['published_margin_a']} A, which is "
                f"why D015 moved it")
            if not (spine.get("threshold_in_gap") and spine.get("correct")):
                problems += 1
                log("  R_SPINE_CONTACT_A no longer sits in the empty region")

        path = MANIFESTS / "refs_abl1.jsonl"
        if args.all_chains:
            # Nothing to diff: the manifest is per-entry and this set is per
            # chain, so every extra chain would read as an addition.
            disagree = [r for r in usable
                        if r.get("dfg_agrees") is False
                        or r.get("alphac_agrees") is False]
            log(f"  all-chains sensitivity: {len(usable)} chains, "
                f"{len(disagree)} disagreeing with KLIFS")
            for row in disagree:
                log(f"    {row['pdb']} {row['chain']}: DFG "
                    f"{row['our_dfg']} against {row['their_dfg']}, aC "
                    f"{row['our_alphac']} against {row['their_alphac']}")
        elif args.verify:
            diff = compare(kinase, _load(path))
            log(f"  verify: {diff['n_fresh']} rebuilt against "
                f"{diff['n_committed']} on disk, {len(diff['changed'])} "
                f"changed, {len(diff['added'])} added, {len(diff['removed'])} removed")
            for row in diff["changed"][:5]:
                log(f"    {row['record']}: {row['fields']}")
            problems += bool(diff["changed"] or diff["added"] or diff["removed"])
        else:
            refs.write_manifest(path, kinase)
            log(f"  wrote {path}")

    if args.only != "kinase":
        log("rebuilding the beta2AR references")
        receptor, p3 = build_receptor()
        log(f"  P3 cleft threshold {p3.get('threshold')} A, leave-one-out "
            f"{p3.get('correct')} of {p3.get('n')} "
            f"({(p3.get('accuracy') or 0):.0%}), needs {P3_MIN:.0%}")
        log(f"  active median {p3.get('active_median')} A, "
            f"inactive median {p3.get('inactive_median')} A")
        if (p3.get("accuracy") or 0) < P3_MIN:
            problems += 1
            log("  P3 is below its pass mark")

        path = MANIFESTS / "refs_b2ar.jsonl"
        if args.verify:
            diff = compare(receptor, _load(path))
            log(f"  verify: {diff['n_fresh']} rebuilt against "
                f"{diff['n_committed']} on disk, {len(diff['changed'])} "
                f"changed, {len(diff['added'])} added, {len(diff['removed'])} removed")
            for row in diff["changed"][:5]:
                log(f"    {row['record']}: {row['fields']}")
            problems += bool(diff["changed"] or diff["added"] or diff["removed"])
        else:
            refs.write_manifest(path, receptor)
            log(f"  wrote {path}")

        # ingest.py carries this threshold as a constant so the explorer does
        # not refit it on every page load. A drift between the two is worth
        # saying out loud rather than leaving to be noticed.
        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "_ingest", ROOT / "tools" / "ingest.py")
        ingest = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ingest)
        stored = ingest.GPCR_CLEFT_THRESHOLD_A
        fitted = p3.get("threshold")
        if fitted is not None and abs(fitted - stored) > 0.01:
            problems += 1
            log(f"  ingest.py uses {stored} A but the fit gives {fitted} A; "
                f"update GPCR_CLEFT_THRESHOLD_A")
        else:
            log(f"  ingest.py's threshold of {stored} A matches the fit")

    return 1 if problems else 0


class _Row:
    """Adapts a manifest dict back to what `refs.agreement` expects.

    `agreement` takes PreparedReference objects and reads two properties off
    them. Rebuilding those objects to re-run a control would mean fetching and
    annotating everything twice, so the dict is wrapped instead.
    """

    def __init__(self, row: dict) -> None:
        self._row = row

    @property
    def status(self):
        return self._row.get("status", "ok")

    @property
    def dfg_agrees(self):
        return self._row.get("dfg_agrees")

    @property
    def alphac_agrees(self):
        return self._row.get("alphac_agrees")


if __name__ == "__main__":
    raise SystemExit(main())
