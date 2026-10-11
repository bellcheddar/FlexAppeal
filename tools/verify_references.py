#!/usr/bin/env python3
"""Resolve every reference through Crossref and write what it actually returns.

Spec section 12.3: "verify every DOI via Crossref at build time; drop any that
fail". The spec lists the works rather than their DOIs, and is explicit that two
of them are uncertain, Shan 2009's venue and which lipid parameter set applies.
So this searches Crossref by title and records what comes back, including the
venue and year, rather than taking the spec's citation as given.

A reference that does not resolve is written with `resolved: false` and kept in
the file rather than deleted, because the About page has to be able to say that
a citation could not be verified. Dropping it silently would leave the method
uncited and look deliberate.

Run once during the build; the output is committed and the page reads it, so
serving a page never touches the network.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "references.json"

CROSSREF = "https://api.crossref.org/works"
# Crossref asks for a contact in the user agent and gives better service for it.
USER_AGENT = "FlexAppeal/1.0 (https://flexappeal.mdeller.com; mailto:marc@marcdeller.com)"

# What each reference is for, so the About page can group them by the thing they
# support rather than alphabetically.
REFERENCES = [
    {"key": "klifs", "expect_year": 2021, "role": "Kinase pocket definitions and reference states",
     "title": "KLIFS: an overhaul after the first 5 years of supporting kinase research",
     "author": "Kanev"},
    {"key": "gpcrdb", "expect_year": 2025, "role": "Receptor generic numbering and reference states",
     "title": "GPCRdb in 2025: adding odorant receptors, data mapper, structure similarity search and models",
     "author": "Taracena Herrera"},
    {"key": "ballesteros", "expect_year": 1995, "role": "Generic numbering scheme",
     "title": "Integrated methods for the construction of three-dimensional models and computational probing of structure-function relations in G protein-coupled receptors",
     "author": "Ballesteros"},
    {"key": "dunbrack", "expect_year": 2019, "role": "DFG state nomenclature",
     "title": "Defining a new nomenclature for the structures of active and inactive kinases",
     "author": "Modi"},
    {"key": "shan", "expect_year": 2009, "role": "DFG-Asp protonation and the Abl flip",
     "title": "A conserved protonation-dependent switch controls drug binding in the Abl kinase",
     "author": "Shan"},
    {"key": "dror", "expect_year": 2011, "expect_venue": "Proceedings of the National Academy", "role": "Beta2AR activation state maps",
     "title": "Activation mechanism of the \u03b22-adrenergic receptor",
     "author": "Dror", "doi": "10.1073/pnas.1110499108"},
    {"key": "rasmussen", "expect_year": 2011, "expect_venue": "Nature", "role": "Active-state beta2AR reference (3P0G)",
     "title": "Structure of a nanobody-stabilized active state of the \u03b22 adrenoceptor",
     "author": "Rasmussen"},
    {"key": "cherezov", "expect_year": 2007, "expect_venue": "Science", "role": "Inactive-state beta2AR reference (2RH1)",
     "title": "High-Resolution Crystal Structure of an Engineered Human \u03b22-Adrenergic G Protein-Coupled Receptor",
     "author": "Cherezov", "doi": "10.1126/science.1150577"},
    {"key": "openmm", "expect_year": 2024,
     "role": "Simulation engine, and structure preparation through PDBFixer",
     "title": "OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials", "venue_hint": "Journal of Physical Chemistry B",
     "author": "Eastman", "doi": "10.1021/acs.jpcb.3c06662"},
    {"key": "mdtraj", "expect_year": 2015, "role": "Trajectory analysis and DSSP",
     "title": "MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories",
     "author": "McGibbon"},
    {"key": "lipid17", "role": "Lipid parameters",
     "title": "Lipid21: Complex Lipid Membrane Simulations with AMBER",
     "author": "Dickson"},
    {"key": "opm", "expect_year": 2012, "role": "Membrane orientation",
     "title": "OPM database and PPM web server: resources for positioning of proteins in membranes",
     "author": "Lomize"},
    {"key": "pdbfixer", "role": "Structure preparation",
     "title": "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics",
     "author": "Eastman"},
    {"key": "openff", "role": "Ligand parameters",
     "title": "Development and Benchmarking of Open Force Field 2.0.0: The Sage Small Molecule Force Field",
     "author": "Boothroyd"},
    {"key": "gemmi", "expect_year": 2023, "role": "Structure file handling",
     "title": "GEMMI and Servalcat restrain REFMAC5",
     "author": "Wojdyr"},
]


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


# Crossref indexes post-publication commentary as its own record, and a search
# for a well-known paper returns the commentary above the paper. Three of these
# references resolved to a "Faculty Opinions recommendation of ..." record on
# the first attempt, which is a different work with a different DOI.
COMMENTARY_PREFIXES = ("faculty opinions recommendation",
                       "faculty of 1000", "correction to", "erratum")

# Only record types that are the paper itself. Filtering on title alone put four
# confidently wrong citations in the file: ACS deposits its supporting
# information as its own Crossref record with the paper's exact title and a DOI
# ending `.s001`, type `component`, and the PDB deposits each entry as type
# `dataset`, so 3P0G resolved to `10.2210/pdb3p0g/pdb` instead of the Nature
# paper describing it. Both are real DOIs for the wrong thing, which is worse
# than not resolving at all.
ACCEPTED_TYPES = {"journal-article", "book-chapter", "proceedings-article",
                  "posted-content", "book-section"}


def _plain(text: str) -> str:
    """A title with Crossref's markup and whitespace removed.

    Crossref returns some titles with inline markup, so the Dror paper arrives
    as "Activation mechanism of the <i>\u03b2</i><sub>2</sub>-adrenergic
    receptor" and compared raw it is 56 per cent similar to its own title. This
    normalises rather than loosens the comparison: the venue, year and author
    checks are unchanged.
    """
    import html
    import re

    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", text))).strip()


def _similarity(one: str, two: str) -> float:
    from difflib import SequenceMatcher

    return SequenceMatcher(None, _plain(one).lower(), _plain(two).lower()).ratio()


def by_doi(entry: dict, *, session) -> dict:
    """Verify a stated DOI, which is what section 12.3 actually asks for.

    Used where Crossref's bibliographic search will not match the title. The
    DOI is resolved and the returned title, venue and year are all checked
    against what the reference claims, so a mistyped or wrong DOI fails here
    rather than appearing on the page as a verified citation.
    """
    doi = entry["doi"]
    try:
        response = session.get(f"{CROSSREF}/{doi}",
                               headers={"User-Agent": USER_AGENT}, timeout=30)
        response.raise_for_status()
        item = response.json().get("message", {})
    except Exception as exc:                                    # noqa: BLE001
        return {**entry, "resolved": False, "note": f"DOI did not resolve: {exc}"[:120]}

    title = _plain((item.get("title") or [""])[0])
    venue = (item.get("container-title") or [""])[0]
    issued = (item.get("issued", {}).get("date-parts") or [[None]])[0]
    year = issued[0] if issued else None
    authors = [str(a.get("family", "")) for a in (item.get("author") or [])]
    score = _similarity(entry["title"], title)

    problems = []
    if score < 0.6:
        problems.append(f"title is {score:.0%} similar: {title[:60]}")
    if entry.get("expect_venue") and entry["expect_venue"].lower() not in venue.lower():
        problems.append(f"venue is {venue!r}")
    if entry.get("expect_year") and year and abs(int(year) - int(entry["expect_year"])) > 1:
        problems.append(f"year is {year}")
    if entry["author"] and not any(
            entry["author"].lower() == name.lower() for name in authors):
        problems.append(f"{entry['author']} is not an author")
    if problems:
        return {**entry, "resolved": False,
                "note": f"DOI {doi} resolves but " + "; ".join(problems)}

    return {**entry, "resolved": True, "doi": doi, "crossref_title": title,
            "title_similarity": round(score, 3), "venue": venue, "year": year,
            "first_author": authors[0] if authors else "",
            "n_authors": len(authors), "type": item.get("type"),
            "verified_by": "doi"}


def resolve(entry: dict, *, session) -> dict:
    """One reference, as Crossref describes it.

    Candidates are ranked by title similarity rather than taking Crossref's own
    first result, because its relevance order puts commentary records above the
    papers they comment on. Whether the expected first author appears is
    recorded rather than used to reject: GPCRdb's database papers have dozens of
    authors and the one the spec names is not always first, and rejecting on
    that produced five false failures out of fifteen.

    What goes in the file is what Crossref returned, including the venue and
    year, which is the point of doing this at all: the spec flags two citations
    as uncertain and this is how they get settled.
    """
    try:
        response = session.get(
            CROSSREF, params={"query.bibliographic": entry["title"], "rows": 10},
            headers={"User-Agent": USER_AGENT}, timeout=30)
        response.raise_for_status()
        items = response.json().get("message", {}).get("items", [])
    except Exception as exc:                                    # noqa: BLE001
        return {**entry, "resolved": False, "note": f"Crossref error: {exc}"[:120]}

    scored = []
    for item in items:
        title = _plain((item.get("title") or [""])[0])
        if not title:
            continue
        if any(title.lower().startswith(prefix) for prefix in COMMENTARY_PREFIXES):
            continue
        if item.get("type") not in ACCEPTED_TYPES:
            continue
        # A component's DOI is the paper's with a suffix, which is the other way
        # these appear when the type is missing or mislabelled.
        if str(item.get("DOI", "")).rsplit(".", 1)[-1].startswith("s00"):
            continue
        scored.append((_similarity(entry["title"], title), title, item))
    if not scored:
        return {**entry, "resolved": False, "note": "no Crossref result"}

    score, title, item = max(scored, key=lambda row: row[0])
    # Below this the match is a different paper rather than a different wording.
    if score < 0.6:
        return {**entry, "resolved": False,
                "note": f"best Crossref match was only {score:.0%} similar: "
                        f"{title[:70]}"}

    authors = [str(a.get("family", "")) for a in (item.get("author") or [])]
    issued = (item.get("issued", {}).get("date-parts") or [[None]])[0]
    year = issued[0] if issued else None
    venue = (item.get("container-title") or [""])[0]

    # A year two or more out, or a venue that is not the one the spec names, is
    # a different work however close the title. Meeting abstracts carry the
    # paper's title in a different journal a year later, which is exactly the
    # shape this catches.
    expect_year = entry.get("expect_year")
    if expect_year and year and abs(int(year) - int(expect_year)) > 1:
        return {**entry, "resolved": False,
                "note": f"title matched {score:.0%} but Crossref gives {year}, "
                        f"not {expect_year}: {venue}"}
    expect_venue = entry.get("expect_venue")
    if expect_venue and expect_venue.lower() not in venue.lower():
        return {**entry, "resolved": False,
                "note": f"title matched {score:.0%} but venue is {venue!r}, "
                        f"not {expect_venue!r}"}
    found = {
        **entry, "resolved": True,
        "doi": item.get("DOI"),
        "crossref_title": title,
        "title_similarity": round(score, 3),
        "venue": venue,
        "year": year,
        "first_author": authors[0] if authors else "",
        "n_authors": len(authors),
        "type": item.get("type"),
    }
    if entry["author"] and not any(
            entry["author"].lower() == name.lower() for name in authors):
        found["note"] = (f"Crossref does not list {entry['author']} as an "
                         f"author; first author is {found['first_author']}")
    return found


def main() -> int:
    import requests

    session = requests.Session()
    out = []
    for entry in REFERENCES:
        found = (by_doi(entry, session=session) if entry.get("doi")
                 else resolve(entry, session=session))
        mark = "ok " if found["resolved"] else "NOT"
        detail = (f"{found.get('venue', '')} {found.get('year', '')} "
                  f"doi:{found.get('doi')}" if found["resolved"]
                  else found.get("note", ""))
        log(f"  {mark} {entry['key']:12s} {detail}")
        out.append(found)
        time.sleep(0.3)                     # Crossref asks for a gentle rate

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    resolved = sum(1 for r in out if r["resolved"])
    log(f"{resolved} of {len(out)} verified; wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
