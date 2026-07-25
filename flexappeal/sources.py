"""Fetching structures from the public databases the Prepare tab offers.

Everything here runs server-side when a bundle is built, not when it is run --
so a user's simulation never depends on the RCSB being reachable eight hours
into a job. The fetched coordinates are embedded in the bundle itself.

Each fetcher enforces a timeout and a size cap and raises ``SourceError`` with a
message aimed at the person who typed the accession, not at a log file.
"""

from __future__ import annotations

import gzip
import urllib.error
import urllib.request
from dataclasses import dataclass

# A structure large enough to exceed this is not something the droplet should be
# introspecting, and almost certainly is not a single protein.
MAX_DOWNLOAD_BYTES = 25 * 1024 * 1024
TIMEOUT_SECONDS = 20

USER_AGENT = "FlexAppeal/0.1 (+https://flexappeal.mdeller.com)"


class SourceError(RuntimeError):
    """A structure could not be retrieved."""


@dataclass(frozen=True)
class FetchResult:
    data: bytes
    filename: str
    source: str
    url: str
    citation: str = ""


def _get(url: str, *, what: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT_SECONDS) as response:
            # Read one byte past the cap so an oversized body is detected rather
            # than silently truncated into a file that parses but is incomplete.
            payload = response.read(MAX_DOWNLOAD_BYTES + 1)
            encoding = response.headers.get("Content-Encoding", "")
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise SourceError(f"{what} was not found. Check the accession.") from None
        raise SourceError(
            f"the server returned HTTP {exc.code} for {what}. It may be temporarily "
            f"unavailable; try again, or upload the file directly."
        ) from None
    except (urllib.error.URLError, TimeoutError) as exc:
        raise SourceError(
            f"could not reach the server to fetch {what} ({exc}). Check your "
            f"connection, or upload the file directly."
        ) from None

    if len(payload) > MAX_DOWNLOAD_BYTES:
        raise SourceError(
            f"{what} is larger than {MAX_DOWNLOAD_BYTES // (1024 * 1024)} MB. "
            f"FlexAppeal is built for single proteins and small complexes."
        )
    if url.endswith(".gz") or encoding == "gzip":
        try:
            payload = gzip.decompress(payload)
        except OSError:
            pass  # Some servers set the header without actually compressing.
    if not payload.strip():
        raise SourceError(f"{what} came back empty.")
    return payload


def fetch_rcsb(pdb_id: str, assembly: str = "asymmetric") -> FetchResult:
    """Fetch a deposited entry from the RCSB PDB.

    The biological assembly is a different file rather than a flag, and is
    frequently what the user actually wants: a physiological dimer deposited as
    a single chain in the asymmetric unit is a real and easy mistake to make.
    """
    pdb_id = pdb_id.strip().upper()
    if assembly == "biological":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
        what = f"PDB entry {pdb_id} (biological assembly)"
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        what = f"PDB entry {pdb_id}"

    try:
        data = _get(url, what=what)
    except SourceError:
        if assembly == "biological":
            # Not every entry has a deposited assembly file; fall back rather
            # than dead-ending, and let the caller's warning explain.
            data = _get(
                f"https://files.rcsb.org/download/{pdb_id}.pdb",
                what=f"PDB entry {pdb_id}",
            )
            return FetchResult(data, f"{pdb_id}.pdb", "rcsb",
                               f"https://www.rcsb.org/structure/{pdb_id}",
                               citation=f"RCSB PDB {pdb_id} (no biological assembly "
                                        f"file deposited; asymmetric unit used)")
        raise

    return FetchResult(data, f"{pdb_id}.pdb", "rcsb",
                       f"https://www.rcsb.org/structure/{pdb_id}",
                       citation=f"RCSB PDB {pdb_id}")


def fetch_opm(pdb_id: str) -> FetchResult:
    """Fetch a membrane-oriented structure from OPM.

    OPM returns the entry rotated so the membrane normal is the z axis, which is
    precisely what ``Modeller.addMembrane`` assumes. Getting this from anywhere
    else means orienting the protein yourself, and building a membrane protein
    into a bilayer sideways is the most common way a membrane setup fails
    silently.
    """
    lower = pdb_id.strip().lower()
    data = _get(
        f"https://opm-assets.storage.googleapis.com/pdb/{lower}.pdb",
        what=f"OPM entry {lower.upper()}",
    )
    return FetchResult(data, f"{lower}_opm.pdb", "opm",
                       f"https://opm.phar.umich.edu/proteins/?search={lower}",
                       citation=f"OPM {lower.upper()} (membrane-oriented)")


def fetch_alphafold(uniprot: str, version: int = 4) -> FetchResult:
    """Fetch a predicted model from the AlphaFold database."""
    accession = uniprot.strip().upper()
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb"
    data = _get(url, what=f"AlphaFold model for {accession}")
    return FetchResult(data, f"AF-{accession}.pdb", "alphafold",
                       f"https://alphafold.ebi.ac.uk/entry/{accession}",
                       citation=f"AlphaFold DB {accession} (v{version})")


def fetch(cfg: dict) -> FetchResult:
    """Dispatch on the configured input source."""
    source = cfg.get("input_source")
    if source == "rcsb":
        return fetch_rcsb(str(cfg.get("pdb_id", "")), str(cfg.get("assembly", "asymmetric")))
    if source == "opm":
        return fetch_opm(str(cfg.get("pdb_id", "")))
    if source == "alphafold":
        return fetch_alphafold(str(cfg.get("uniprot_id", "")))
    raise SourceError(f"{source!r} is not a fetchable source; upload the file instead.")


# Warnings that depend on where a structure came from rather than what is in it.
SOURCE_WARNINGS = {
    "alphafold": (
        "AlphaFold models contain no waters, no ligands and no cofactors, and their "
        "loop conformations and relative domain orientations are the least reliable "
        "part of the prediction. Check the per-residue confidence before simulating."
    ),
    "opm": (
        "OPM coordinates are oriented for membrane insertion, with the bilayer normal "
        "along z. Do not re-align or re-centre them before building the membrane."
    ),
}
