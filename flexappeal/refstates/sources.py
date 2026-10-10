"""KLIFS and GPCRdb: the two catalogues this build takes its residue identities from.

Nothing here types a residue number. A kinase motif is a KLIFS pocket position
and a receptor motif is a GPCRdb generic number, and both are resolved against
the live catalogue for the structure in hand. That is the whole point of using
them: ABL1 alone has two isoforms nineteen residues apart and every PDB entry
numbers its construct differently, so a hard-coded "Asp381" is right for one file
and silently wrong for the next.

What this module is, and is not. It fetches and filters: kinase and receptor
identity, the reference structure lists, and the GPCRdb residue table. The motif
annotation itself is BoltzMaker's `sse_comparison`, vendored and wrapped rather
than reimplemented (spec section 4.1), and the metrics live beside it.

Both services are somebody else's infrastructure, so every call is cached on disk,
rate limited to two requests a second, and retried on the failures that are worth
retrying. A build that re-runs should hit the network once.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# --- the catalogues ---------------------------------------------------------

KLIFS_ROOT = "https://klifs.net/api"
GPCRDB_ROOT = "https://gpcrdb.org/services"

# Two requests a second each, from spec section 2.2. These are public academic
# services and a build that hammers them is a build that gets blocked.
MIN_INTERVAL_S = 0.5

# --- the filters, from spec sections 4.2 and 4.3 ---------------------------

KINASE_MAX_RESOLUTION_A = 2.8
KINASE_MIN_QUALITY = 6.0
KINASE_MAX_MISSING_POCKET = 5
GPCR_MAX_RESOLUTION_A = 3.6

# ABL1's KLIFS id is asserted rather than hard-coded: the swagger lists 392 and
# section 4.2 says to check the lookup agrees rather than trusting either.
EXPECTED_ABL1_KINASE_ID = 392


class SourceError(RuntimeError):
    """A catalogue could not be reached or returned something unusable."""


@dataclass
class Ensemble:
    """One state's reference structures, with the filter that produced them.

    `rejected` is carried rather than dropped. A reference set is evidence, and
    "forty-five inactive structures" means something different if four hundred
    were considered than if forty-six were.
    """

    name: str
    state: str
    entries: list[dict[str, Any]] = field(default_factory=list)
    considered: int = 0
    rejected: dict[str, int] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.entries)

    @property
    def pdb_codes(self) -> list[str]:
        key = "pdb" if self.entries and "pdb" in self.entries[0] else "pdb_code"
        return sorted({str(e[key]).upper() for e in self.entries})


class Catalogue:
    """A cached, rate-limited reader for one of the two services.

    The cache is keyed on the full request rather than on an endpoint name, so a
    different kinase or a different receptor cannot collide with an earlier
    answer. `refresh` bypasses it, which is what a build does when a catalogue
    has been updated and the stored copy is the thing under suspicion.
    """

    def __init__(self, cache_dir: Path, *, refresh: bool = False,
                 min_interval_s: float = MIN_INTERVAL_S) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.refresh = refresh
        self.min_interval_s = min_interval_s
        self._last_call = 0.0

    # -- plumbing ----------------------------------------------------------

    def _cache_path(self, url: str, params: dict | None) -> Path:
        import hashlib

        key = url + "?" + json.dumps(params or {}, sort_keys=True)
        digest = hashlib.sha256(key.encode()).hexdigest()[:20]
        return self.cache_dir / f"{digest}.json"

    def _wait(self) -> None:
        elapsed = time.monotonic() - self._last_call
        if elapsed < self.min_interval_s:
            time.sleep(self.min_interval_s - elapsed)
        self._last_call = time.monotonic()

    def get(self, url: str, params: dict | None = None) -> Any:
        """One GET, cached on disk and retried on the failures worth retrying."""
        path = self._cache_path(url, params)
        if path.is_file() and not self.refresh:
            try:
                return json.loads(path.read_text())
            except json.JSONDecodeError:
                path.unlink(missing_ok=True)          # a truncated write

        import requests
        from tenacity import (retry, retry_if_exception_type,
                              stop_after_attempt, wait_exponential)

        @retry(stop=stop_after_attempt(4),
               wait=wait_exponential(multiplier=1, min=1, max=20),
               retry=retry_if_exception_type(
                   (requests.ConnectionError, requests.Timeout)),
               reraise=True)
        def _fetch() -> Any:
            self._wait()
            response = requests.get(url, params=params, timeout=45)
            # A 4xx is not retried: a wrong kinase name will be wrong again.
            if response.status_code >= 500:
                raise requests.ConnectionError(f"{url} returned {response.status_code}")
            if not response.ok:
                raise SourceError(f"{url} returned {response.status_code}")
            return response.json()

        payload = _fetch()
        path.write_text(json.dumps(payload))
        return payload


# ===========================================================================
#  KLIFS: kinases
# ===========================================================================


def kinase_id(catalogue: Catalogue, name: str = "ABL1",
              species: str = "Human") -> int:
    """The KLIFS numeric id for a kinase, looked up rather than assumed.

    Section 4.2 is specific about this: the swagger lists ABL1 as 392 and the
    rule is to assert the lookup agrees, not to write 392 down. A catalogue that
    renumbered would otherwise produce a build that fetched the wrong kinase's
    structures and said nothing.
    """
    rows = catalogue.get(f"{KLIFS_ROOT}/kinase_ID",
                         {"kinase_name": name, "species": species})
    if not isinstance(rows, list) or not rows:
        raise SourceError(f"KLIFS returned no kinase for {name}/{species}")
    exact = [r for r in rows if str(r.get("name", "")).upper() == name.upper()]
    chosen = (exact or rows)[0]
    found = int(chosen["kinase_ID"])
    if name.upper() == "ABL1" and found != EXPECTED_ABL1_KINASE_ID:
        raise SourceError(
            f"KLIFS resolved ABL1 to kinase_ID {found}, not the documented "
            f"{EXPECTED_ABL1_KINASE_ID}. The catalogue has changed under this "
            f"build and the reference ensembles would be another kinase's.")
    return found


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def kinase_references(catalogue: Catalogue, kid: int) -> dict[str, Ensemble]:
    """The ABL1 reference ensembles, filtered and split by KLIFS's own labels.

    Three sets, not two. Section 4.2 keeps the Src-like inactive conformation
    (DFG-in with the αC helix out) as its own category rather than folding it
    into the inactive ensemble, because in Abl it is a distinct state and a frame
    sitting in it is not the same finding as a frame sitting in DFG-out.
    """
    rows = catalogue.get(f"{KLIFS_ROOT}/structures_list", {"kinase_ID": kid})
    if not isinstance(rows, list) or not rows:
        raise SourceError(f"KLIFS returned no structures for kinase {kid}")

    kept, rejected = [], {}
    for row in rows:
        reason = None
        if row.get("species") != "Human":
            reason = "not human"
        elif not 0 < _as_float(row.get("resolution"), 99) <= KINASE_MAX_RESOLUTION_A:
            reason = f"resolution over {KINASE_MAX_RESOLUTION_A} A"
        elif _as_float(row.get("quality_score")) < KINASE_MIN_QUALITY:
            reason = f"quality below {KINASE_MIN_QUALITY}"
        elif _as_float(row.get("missing_residues"), 99) > KINASE_MAX_MISSING_POCKET:
            reason = f"over {KINASE_MAX_MISSING_POCKET} missing pocket residues"
        if reason:
            rejected[reason] = rejected.get(reason, 0) + 1
        else:
            kept.append(row)

    def ensemble(label: str, state: str, predicate) -> Ensemble:
        return Ensemble(name=label, state=state,
                        entries=[r for r in kept if predicate(r)],
                        considered=len(rows), rejected=dict(rejected))

    return {
        "active": ensemble(
            "ABL1 active-like", "active",
            lambda r: r.get("DFG") == "in" and r.get("aC_helix") == "in"),
        "inactive": ensemble(
            "ABL1 inactive (DFG-out)", "inactive",
            lambda r: r.get("DFG") == "out"),
        "src_like": ensemble(
            "ABL1 Src-like inactive", "src_like",
            lambda r: r.get("DFG") == "in" and r.get("aC_helix") == "out"),
    }


# ===========================================================================
#  GPCRdb: receptors
# ===========================================================================


def receptor_residues(catalogue: Catalogue,
                      entry_name: str = "adrb2_human") -> list[dict[str, Any]]:
    """The receptor's residue table: sequence number to generic number.

    This is the mapping every GPCR motif in section 5.5 is resolved through.
    `display_generic_number` is GPCRdb's structure-based scheme and
    `alternative_generic_numbers` carries Ballesteros-Weinstein, which is the one
    the labels use.
    """
    rows = catalogue.get(f"{GPCRDB_ROOT}/residues/extended/{entry_name}/")
    if not isinstance(rows, list) or not rows:
        raise SourceError(f"GPCRdb returned no residues for {entry_name}")
    return rows


def receptor_references(catalogue: Catalogue,
                        entry_name: str = "adrb2_human") -> dict[str, Ensemble]:
    """The β2AR reference ensembles, split by GPCRdb's own state call.

    Intermediate is kept as a third set where the catalogue has any: section 4.3
    wants them on the state maps even though they are neither endpoint. β2AR
    currently has none, which is recorded rather than assumed away.
    """
    rows = catalogue.get(f"{GPCRDB_ROOT}/structure/protein/{entry_name}/")
    if not isinstance(rows, list) or not rows:
        raise SourceError(f"GPCRdb returned no structures for {entry_name}")

    kept, rejected = [], {}
    for row in rows:
        if _as_float(row.get("resolution"), 99) > GPCR_MAX_RESOLUTION_A:
            reason = f"resolution over {GPCR_MAX_RESOLUTION_A} A"
            rejected[reason] = rejected.get(reason, 0) + 1
        else:
            kept.append(row)

    def ensemble(label: str, state: str) -> Ensemble:
        return Ensemble(name=label, state=state,
                        entries=[r for r in kept if r.get("state") == state.title()],
                        considered=len(rows), rejected=dict(rejected))

    return {
        "active": ensemble("β2AR active", "active"),
        "inactive": ensemble("β2AR inactive", "inactive"),
        "intermediate": ensemble("β2AR intermediate", "intermediate"),
    }
