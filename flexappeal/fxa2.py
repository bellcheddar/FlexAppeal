"""Writing and reading the ``.fxa`` version 2 package file.

Version 1 (``flexappeal/fxa.py``) is one system and one replicate, and DeepFlex
renders it as a flexibility run with no grade. Version 2 is a different kind of
file rather than a bigger version of the same one: it carries a whole package of
receptor-ligand complexes with the six sub-score inputs, so DeepFlex renders it
as a fully graded report.

The reason it has to be a package is that the central panels are comparative.
``S_E`` is a percentile within a receptor's own cohort, and the league table and
the two rank columns are places in that cohort, so a file with one complex in it
has no cohort and those quantities are undefined rather than zero. FlexAppeal runs
one system per job, which is why ``flexappeal.package`` collects several finished
runs into one file instead of a job writing one on its own.

Specified in ``DEEPFLEX_FXA2_SPEC.md`` in the repo root. This module owns the
format: the identifiers, the member layout, the archive limits, the chemistry
bans, and a reader that can read back everything the writer wrote. The reader
exists for the round trip rather than for the web app, which still reads version
1: a writer that cannot read its own output is a writer with a unit bug.

**Units.** Every length in this file is angstrom and every time is nanoseconds,
and the ``_a`` and ``_ns`` suffixes are the contract. Version 1 reported
nanometres. MDTraj works in nanometres, so every number crossing from
``flexappeal.measure`` into here has already been multiplied by ten; getting that
wrong puts every length in the report out by a factor of ten, which is invisible
until somebody checks a number by hand.
"""

from __future__ import annotations

import hashlib
import hmac
import io
import json
import os
import re
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

# The layout this module writes and reads. A reader seeing a version it does not
# know must refuse by number rather than attempt a partial parse.
FXA2_VERSION = 2

# --- archive limits, section 1 ---------------------------------------------
# These are DeepFlex's own numbers, enforced on read at the far end. Checking them
# here as well means a package that would be refused on upload fails on the
# machine that built it, where the trajectory it came from is still to hand.
MAX_MEMBERS = 20_000
MAX_TOTAL_UNCOMPRESSED = 2_500 * 1024 * 1024    # 2.5 GB
MAX_SINGLE_MEMBER = 512 * 1024 * 1024
MAX_COMPRESSION_RATIO = 200
# Below this, a member's ratio is meaningless: a 200-byte JSON object of repeated
# key names compresses 10:1 without being an attack.
MIN_SIZE_FOR_RATIO_CHECK = 64 * 1024

# --- what must never be in the file, section 8 ------------------------------
# Connectivity and atom typing are what must not travel. The coordinates of one
# binding site, under a hashed filename, are not: a crop says where atoms sat in
# one pocket, a structure string says what the molecule is, and that is the line.
#
# DeepFlex discards these on ingest and tells the uploader it did, so writing one
# costs effort and buys a warning. Refusing here instead means the mistake is
# caught by the test suite rather than by somebody's upload.
CHEMISTRY_SUFFIXES = frozenset({
    ".sdf", ".sd", ".mol", ".mol2", ".smi", ".smiles", ".inchi", ".pdbqt",
    ".frcmod", ".lib", ".prmtop", ".inpcrd", ".rst7", ".xtc", ".dcd", ".trr",
})
CHEMISTRY_ROOTS = ("smiles", "inchi", "molblock", "molfile", "tripos", "sdf")

COMPLEX_ID_RE = re.compile(r"^[0-9a-f]{16}_[A-Z0-9_-]+$")
RECEPTOR_ID_RE = re.compile(r"^[A-Z0-9_-]+$")

# Residues that must not appear in a pocket crop. Section 6.1 bans waters,
# counter-ions and lipids, and conformance check 7 looks for the first two.
SOLVENT_RESNAMES = frozenset({"HOH", "WAT", "TIP3", "SOL", "H2O"})
ION_RESNAMES = frozenset({
    "NA", "CL", "K", "MG", "CA", "ZN", "SOD", "CLA", "POT", "BR", "IOD", "LI",
    "CS", "RB", "F", "MN", "FE", "CU", "NI", "CO",
})

# The ligand residue name in every pocket crop. Fixed by section 6.1, and the 3D
# viewer selects on it.
LIGAND_RESNAME = "LIG"

# Caps on the payload arrays, section 5. Reducing before writing rather than after
# reading, because a 2,000-frame 2D RMSD is four million floats and the panel is a
# few hundred pixels wide.
MAX_MATRIX_SIDE = 500
MAX_BARCODE_ROWS = 28
MAX_BARCODE_COLS = 400


class Fxa2Error(ValueError):
    """A package could not be written or read. The message is shown to the user."""


# ===========================================================================
#  Identifiers, section 2
# ===========================================================================


def new_salt() -> bytes:
    """32 random bytes, one per package file.

    Never recorded in the file. These ids exist so that URLs, filenames and
    database rows route on something opaque and stable rather than on a compound
    name, and stability is only needed *within* one package, so a fresh salt per
    file costs nothing. Writing the salt down is the one thing that could turn an
    id back into a name.
    """
    return os.urandom(32)


def ligand_id(salt: bytes, compound_name: str) -> str:
    """16 lowercase hex characters, salted.

    An unsalted digest of a short compound code is reversed by brute force in
    seconds, which is the whole reason for the HMAC.
    """
    if not salt:
        raise Fxa2Error("a ligand id needs a salt; an unsalted digest is reversible")
    if not compound_name:
        raise Fxa2Error("a ligand id needs a compound name")
    digest = hmac.new(salt, compound_name.encode("utf-8"), hashlib.sha256)
    return digest.hexdigest()[:16]


def receptor_token(*candidates: str | None) -> str:
    """A short readable receptor id, from the first candidate that yields one.

    Callers pass their candidates in order of how much each can be trusted. It
    appears in URLs and becomes part of a filename, so it is uppercased and
    reduced to the character class section 2 allows.
    """
    for candidate in candidates:
        if not candidate:
            continue
        token = re.sub(r"[^A-Za-z0-9_-]+", "-", str(candidate)).strip("-_").upper()
        if token:
            return token[:24]
    raise Fxa2Error("no usable receptor identifier: pass a PDB id, a filename or a job name")


def complex_id(ligand: str, receptor: str) -> str:
    """``<ligand_id>_<receptor_id>``, validated, because it becomes a filename."""
    cid = f"{ligand}_{receptor}"
    if not COMPLEX_ID_RE.match(cid):
        raise Fxa2Error(
            f"{cid!r} is not a usable complex id: it must be sixteen hex characters, "
            f"an underscore, then an uppercase receptor token"
        )
    return cid


# ===========================================================================
#  Writing
# ===========================================================================


def _check_member_name(name: str) -> None:
    """Refuse a member name that a hostile-archive check at the far end would."""
    if not name or name.endswith("/"):
        raise Fxa2Error(f"{name!r} is not a usable member name")
    if "\x00" in name:
        raise Fxa2Error(f"{name!r} contains a null byte")
    if "\\" in name:
        raise Fxa2Error(f"{name!r} uses backslash separators")
    if name.startswith("/") or (len(name) > 1 and name[1] == ":"):
        raise Fxa2Error(f"{name!r} is an absolute path")
    if ".." in name.split("/"):
        raise Fxa2Error(f"{name!r} escapes the archive root")
    basename = name.rsplit("/", 1)[-1]
    suffix = "." + basename.rsplit(".", 1)[-1].lower() if "." in basename else ""
    if suffix in CHEMISTRY_SUFFIXES:
        raise Fxa2Error(
            f"{name!r} carries chemistry or a trajectory, which section 8 of the "
            f"spec bans: connectivity and atom typing must not travel, and DeepFlex "
            f"discards the member and says so"
        )


def _check_no_chemistry_fields(obj: Any, path: str = "manifest") -> None:
    """No field name at any depth may promise a structure string.

    A substring test on a lowercased name, so a prefixed or suffixed spelling such
    as a canonical or isomeric variant is caught without listing every one.
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            lowered = str(key).lower()
            for root in CHEMISTRY_ROOTS:
                if root in lowered:
                    raise Fxa2Error(
                        f"{path}.{key} names a structure string, which section 8 "
                        f"of the spec bans. A compound name says which compound; "
                        f"a structure says what it is, and that is the line."
                    )
            _check_no_chemistry_fields(value, f"{path}.{key}")
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            _check_no_chemistry_fields(value, f"{path}[{index}]")


def _json_bytes(obj: Any) -> bytes:
    """JSON with numpy scalars reduced to Python ones.

    ``json.dumps`` refuses a ``numpy.float32``, and every number here has come
    through numpy. Failing on the whole document because one value is the wrong
    flavour of float would be a poor trade.
    """
    def default(value: Any) -> Any:
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        raise TypeError(f"{type(value).__name__} is not JSON-serialisable")

    return json.dumps(obj, default=default).encode("utf-8")


@dataclass
class Package:
    """A version 2 package, accumulated in memory and then written as one zip.

    Held in memory rather than staged on disk because the members are a few
    hundred kilobytes each: the pocket crops and the decimated multi-model PDBs
    are the large ones, and a study of a few dozen complexes is tens of megabytes.
    A package large enough for that to matter would trip the 2.5 GB limit first.
    """

    package_name: str
    flexappeal_version: str
    engine: str = "unknown"
    run_settings: dict[str, Any] = field(default_factory=dict)
    salt: bytes = field(default_factory=new_salt)
    kind: str = "graded"
    created_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))

    ligand_names: dict[str, str] = field(default_factory=dict)
    receptors: dict[str, dict[str, Any]] = field(default_factory=dict)
    complexes: list[str] = field(default_factory=list)
    members: dict[str, bytes] = field(default_factory=dict)

    # --- identifiers ------------------------------------------------------

    def register_ligand(self, compound_name: str) -> str:
        """Hash a compound name and record it in clear, as section 2 allows."""
        lid = ligand_id(self.salt, compound_name)
        self.ligand_names[lid] = compound_name
        return lid

    def register_receptor(self, receptor: str, **facts: Any) -> str:
        if not RECEPTOR_ID_RE.match(receptor):
            raise Fxa2Error(f"{receptor!r} is not a usable receptor id")
        entry = self.receptors.setdefault(receptor, {"label": receptor})
        entry.update({k: v for k, v in facts.items() if v not in (None, "")})
        return receptor

    # --- members ----------------------------------------------------------

    def add_bytes(self, name: str, payload: bytes) -> None:
        _check_member_name(name)
        if name in self.members:
            raise Fxa2Error(f"{name} is already in this package")
        if len(payload) > MAX_SINGLE_MEMBER:
            raise Fxa2Error(f"{name} is {len(payload) / 1e6:.0f} MB, over the member limit")
        self.members[name] = payload

    def add_json(self, name: str, obj: Any) -> None:
        _check_no_chemistry_fields(obj, name)
        self.add_bytes(name, _json_bytes(obj))

    def add_text(self, name: str, text: str) -> None:
        self.add_bytes(name, text.encode("utf-8"))

    def add_arrays(self, name: str, **arrays: np.ndarray) -> None:
        """One ``savez_compressed`` member.

        Keys, dtypes and shapes are the contract: an endpoint reads a key by name,
        and a key it cannot find is a panel that does not draw.
        """
        buffer = io.BytesIO()
        np.savez_compressed(buffer, **arrays)
        self.add_bytes(name, buffer.getvalue())

    def add_complex(self, cid: str) -> None:
        if cid in self.complexes:
            raise Fxa2Error(f"{cid} is already in this package")
        if not COMPLEX_ID_RE.match(cid):
            raise Fxa2Error(f"{cid!r} is not a usable complex id")
        self.complexes.append(cid)

    # --- the manifest -----------------------------------------------------

    def manifest(self) -> dict[str, Any]:
        return {
            "fxa_version": FXA2_VERSION,
            "flexappeal_version": self.flexappeal_version,
            "kind": self.kind,
            "package_name": self.package_name,
            "created_at": self.created_at,
            "engine": self.engine,
            "run_settings": dict(self.run_settings),
            "ligand_names": dict(self.ligand_names),
            "receptors": {k: dict(v) for k, v in sorted(self.receptors.items())},
            "complexes": list(self.complexes),
            # True, and the salt itself is deliberately absent: see new_salt.
            "salt_is_per_bundle": True,
            # Must be null. Section 11 and section 8 between them: there is nowhere
            # in this file for chemistry to live.
            "chemistry": None,
        }

    # --- writing ----------------------------------------------------------

    def write(self, path: Path | str) -> Path:
        """Write the zip, after checking everything the far end checks on read."""
        path = Path(path)
        manifest = self.manifest()
        _check_no_chemistry_fields(manifest)

        for cid in self.complexes:
            missing = [
                member for member in (
                    f"metrics/{cid}.json", f"pocket/{cid}.json",
                    f"traj/{cid}_pocket.pdb",
                )
                if member not in self.members
            ]
            if missing:
                raise Fxa2Error(
                    f"{cid} is listed in the manifest without {', '.join(missing)}; "
                    f"those three are required for every complex"
                )
            lid = cid.split("_", 1)[0]
            if lid not in self.ligand_names:
                raise Fxa2Error(f"{cid}'s ligand half is not in manifest.ligand_names")

        payload = {"manifest.json": json.dumps(manifest, indent=2).encode("utf-8")}
        payload.update(self.members)

        if len(payload) > MAX_MEMBERS:
            raise Fxa2Error(f"{len(payload)} members, over the limit of {MAX_MEMBERS}")
        total = sum(len(blob) for blob in payload.values())
        if total > MAX_TOTAL_UNCOMPRESSED:
            raise Fxa2Error(
                f"this package expands to {total / 1e9:.1f} GB, over the 2.5 GB limit")

        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
            for name in sorted(payload):
                archive.writestr(name, payload[name])
        return path


# ===========================================================================
#  Reading
# ===========================================================================


@dataclass
class ReadPackage:
    """A package read back from disk, for the round trip and for inspection."""

    manifest: dict[str, Any]
    members: dict[str, bytes]
    warnings: list[str] = field(default_factory=list)

    @property
    def complexes(self) -> list[str]:
        return list(self.manifest.get("complexes") or [])

    @property
    def ligand_names(self) -> dict[str, str]:
        return dict(self.manifest.get("ligand_names") or {})

    @property
    def run_settings(self) -> dict[str, Any]:
        return dict(self.manifest.get("run_settings") or {})

    def json(self, name: str) -> Any | None:
        blob = self.members.get(name)
        if blob is None:
            return None
        try:
            return json.loads(blob)
        except json.JSONDecodeError as exc:
            raise Fxa2Error(f"{name} is not valid JSON: {exc}") from None

    def metrics(self, cid: str) -> dict[str, Any] | None:
        return self.json(f"metrics/{cid}.json")

    def interactions(self, cid: str) -> dict[str, Any] | None:
        return self.json(f"metrics/{cid}.interactions.json")

    def pocket(self, cid: str) -> dict[str, Any] | None:
        return self.json(f"pocket/{cid}.json")

    def summary(self, cid: str) -> list[dict[str, Any]] | None:
        return self.json(f"summary/{cid}.json")

    def arrays(self, name: str) -> dict[str, np.ndarray] | None:
        """One ``.npz`` member's arrays.

        ``allow_pickle`` stays off for every member except the barcode, whose
        ``labels`` key is an object array by specification (section 5.3) and
        cannot be read without it. An object array in an untrusted npz is
        arbitrary code execution on load, so it is enabled for the one member
        that needs it and the one key that is allowed to be an object.
        """
        blob = self.members.get(name)
        if blob is None:
            return None
        pickled = name.startswith("barcode/")
        try:
            with np.load(io.BytesIO(blob), allow_pickle=pickled) as data:
                return {key: data[key] for key in data.files}
        except (ValueError, OSError, zipfile.BadZipFile) as exc:
            raise Fxa2Error(f"{name} could not be read: {exc}") from None

    def pocket_pdb(self, cid: str) -> str | None:
        blob = self.members.get(f"traj/{cid}_pocket.pdb")
        return None if blob is None else blob.decode("utf-8", "replace")


def read(source: Path | str | bytes) -> ReadPackage:
    """Read a version 2 package, refusing anything else by number."""
    content = Path(source).read_bytes() if not isinstance(source, bytes) else source
    if not content:
        raise Fxa2Error("that file is empty.")

    try:
        archive = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile:
        raise Fxa2Error("this file is not a valid .fxa package.") from None

    warnings: list[str] = []
    with archive:
        entries = archive.infolist()
        if len(entries) > MAX_MEMBERS:
            raise Fxa2Error(f"this archive declares {len(entries)} members, over the limit.")
        total = 0
        for entry in entries:
            if entry.is_dir():
                continue
            _check_member_name(entry.filename)
            mode = entry.external_attr >> 16
            if mode and (mode & 0xF000) == 0xA000:
                raise Fxa2Error("this archive contains a symbolic link; refusing to read it.")
            if entry.file_size > MAX_SINGLE_MEMBER:
                raise Fxa2Error(f"{entry.filename} is larger than any member should be.")
            if (entry.file_size >= MIN_SIZE_FOR_RATIO_CHECK and entry.compress_size > 0
                    and entry.file_size / entry.compress_size > MAX_COMPRESSION_RATIO):
                raise Fxa2Error(
                    f"{entry.filename} claims a compression ratio real data does not "
                    f"achieve. Refusing to expand it."
                )
            total += entry.file_size
        if total > MAX_TOTAL_UNCOMPRESSED:
            raise Fxa2Error(f"this archive expands to {total / 1e9:.1f} GB, over the limit.")

        members: dict[str, bytes] = {}
        for entry in entries:
            if entry.is_dir():
                continue
            data = archive.read(entry)
            # The declared size is attacker-controlled; check what came out.
            if len(data) > MAX_SINGLE_MEMBER:
                raise Fxa2Error(f"{entry.filename} is larger than its header claimed.")
            members[entry.filename] = data

    if "manifest.json" not in members:
        raise Fxa2Error("this package has no manifest.json.")
    try:
        manifest = json.loads(members["manifest.json"])
    except json.JSONDecodeError as exc:
        raise Fxa2Error(f"the package manifest is corrupt: {exc}") from None
    if not isinstance(manifest, dict):
        raise Fxa2Error("the package manifest is not in the expected form.")

    version = manifest.get("fxa_version")
    if version != FXA2_VERSION:
        raise Fxa2Error(
            f"this file declares fxa_version {version!r} and the package reader "
            f"handles {FXA2_VERSION}. Version 1 files are read by flexappeal.fxa."
        )

    if manifest.get("chemistry") is not None:
        warnings.append("manifest.chemistry is not null; DeepFlex will discard it.")

    return ReadPackage(manifest=manifest, members=members, warnings=warnings)


# ===========================================================================
#  Conformance, section 10
# ===========================================================================

# The six sub-score inputs plus the quality flags, as section 4 requires them of a
# complex with status "ok" in a graded file. Asserted as a list of what is missing
# rather than as a boolean, so a failure says which.
REQUIRED_METRIC_KEYS = (
    "complex_id", "receptor", "ligand_id", "status", "n_frames",
    "equil_start_frame", "equil_method",
    "ligand_rmsd_median", "ligand_rmsd_p90", "ligand_rmsd_max",
    "ligand_rmsd_drift_slope", "inplace_rmsd_median",
    "n_native_contacts", "in_pocket_fraction", "n_exit_events",
    "com_drift_mean_a", "com_drift_max_a", "pocket_offset_mean_a",
    "conv_js_mean", "n_pose_clusters", "dominant_cluster_population",
    "changed_pose_cluster", "pose_clusters",
    "pbc_ok", "extra_ligand_fragments", "box_volume_a3",
    "receptor_rmsd_mean", "ligand_rmsf_mean", "ligand_rmsf_max",
)

# Keys that must lie in 0 to 1 inclusive wherever they appear.
FRACTION_KEYS = (
    "in_pocket_fraction", "dominant_cluster_population", "tanimoto_t0_mean",
    "pose_memory_mean", "key_interaction_occupancy_mean", "hbond_occupancy_mean",
)

_SERIES_REQUIRED = {
    "time_ns": "float64",
    "pose_rmsd_a": "float32",
    "inplace_rmsd_a": "float32",
    "receptor_rmsd_a": "float32",
    "com_drift_a": "float32",
    "pocket_offset_a": "float32",
    "contacts": "int16",
    "in_pocket": "int8",
    "ligand_rmsf_a": "float32",
    "equil_start_frame": "int64",
    "native_contact_threshold": "int64",
}

_SERIES_SHARE_LENGTH = (
    "pose_rmsd_a", "inplace_rmsd_a", "receptor_rmsd_a", "com_drift_a",
    "pocket_offset_a", "contacts", "in_pocket",
)


def conformance(package: ReadPackage) -> list[str]:
    """The fourteen checks in section 10. An empty list is a conformant file.

    Shipped as a test rather than only as a tool: it is a cheaper proxy for
    DeepFlex's own validator and it runs without DeepFlex being checked out.
    """
    problems: list[str] = []
    manifest = package.manifest
    members = package.members

    # 1. The archive opens, has a manifest, and declares version 2. `read` has
    #    already refused anything else, so reaching here settles it.
    if manifest.get("fxa_version") != FXA2_VERSION:
        problems.append(f"fxa_version is {manifest.get('fxa_version')!r}, not 2")

    cids = package.complexes
    if not cids:
        problems.append("manifest.complexes is empty")

    for cid in cids:
        # 2. Every id has its three required members.
        for member in (f"metrics/{cid}.json", f"pocket/{cid}.json",
                       f"traj/{cid}_pocket.pdb"):
            if member not in members:
                problems.append(f"{member} is missing")

        # 3. The id's shape, and its ligand half named in clear.
        if not COMPLEX_ID_RE.match(cid):
            problems.append(f"{cid} does not match ^[0-9a-f]{{16}}_[A-Z0-9_-]+$")
        elif cid.split("_", 1)[0] not in package.ligand_names:
            problems.append(f"{cid}'s ligand half is absent from manifest.ligand_names")

    # 4 and 5. Member names: traversal, and the chemistry suffixes.
    for name in members:
        try:
            _check_member_name(name)
        except Fxa2Error as exc:
            problems.append(str(exc))

    # 6. No manifest field name promises a structure string, at any depth.
    try:
        _check_no_chemistry_fields(manifest)
    except Fxa2Error as exc:
        problems.append(str(exc))

    for cid in cids:
        problems.extend(_check_pocket_pdb(package, cid))          # 7
        problems.extend(_check_arrays(package, cid))              # 8, 11, 12
        problems.extend(_check_metrics(package, cid))             # 9, 10, 13

    return problems


def _check_pocket_pdb(package: ReadPackage, cid: str) -> list[str]:
    """Check 7: the ligand is ``LIG`` in ``ATOM`` records, and no solvent.

    The ``TER`` records are checked too, and that is not pedantry. A ``TER``
    carries the residue name in the same columns as an ``ATOM`` record, so
    renaming only the atom lines leaves the real het code in the file. That exact
    bug shipped once on the DeepFlex side.
    """
    text = package.pocket_pdb(cid)
    if text is None:
        return []

    problems: list[str] = []
    atom_resnames: set[str] = set()
    ter_resnames: set[str] = set()
    for line in text.splitlines():
        record = line[:6].strip()
        if record not in ("ATOM", "HETATM", "TER"):
            continue
        resname = line[17:20].strip().upper()
        if not resname:
            continue
        if record == "TER":
            ter_resnames.add(resname)
        else:
            atom_resnames.add(resname)
            if resname == LIGAND_RESNAME and record == "HETATM":
                problems.append(
                    f"traj/{cid}_pocket.pdb names the ligand {LIGAND_RESNAME} on a "
                    f"HETATM record; section 6.1 asks for ATOM")
        if resname in SOLVENT_RESNAMES:
            problems.append(f"traj/{cid}_pocket.pdb contains {resname}; no waters")
        if resname in ION_RESNAMES:
            problems.append(f"traj/{cid}_pocket.pdb contains the ion {resname}")

    if LIGAND_RESNAME not in atom_resnames:
        problems.append(
            f"traj/{cid}_pocket.pdb has no {LIGAND_RESNAME} residue, so the viewer "
            f"has no ligand to show")

    # A TER naming a residue that no ATOM record names is the renaming bug: the
    # atoms were relabelled and the TER beneath them was not, so the real het code
    # is still in the file.
    for stray in sorted(ter_resnames - atom_resnames):
        problems.append(
            f"traj/{cid}_pocket.pdb has a TER record for residue {stray}, which no "
            f"ATOM record names: the ligand was renamed and its TER was not")
    return problems


def _check_arrays(package: ReadPackage, cid: str) -> list[str]:
    """Checks 8, 11 and 12: the npz keys, dtypes, shared lengths and caps."""
    problems: list[str] = []

    series = package.arrays(f"series/{cid}.npz")
    if series is not None:
        for key, dtype in _SERIES_REQUIRED.items():
            if key not in series:
                problems.append(f"series/{cid}.npz has no {key}")
            elif series[key].dtype.name != dtype:
                problems.append(
                    f"series/{cid}.npz {key} is {series[key].dtype.name}, not {dtype}")
        if "time_ns" in series:
            n = series["time_ns"].shape[0]
            for key in _SERIES_SHARE_LENGTH:
                if key in series and series[key].shape[0] != n:
                    problems.append(
                        f"series/{cid}.npz {key} has {series[key].shape[0]} values "
                        f"against time_ns's {n}")

    panels = package.arrays(f"panels/{cid}.npz")
    if panels is not None:
        if "rmsf_resid" in panels and "rmsf_a" in panels:
            if panels["rmsf_resid"].shape != panels["rmsf_a"].shape:
                problems.append(
                    f"panels/{cid}.npz rmsf_resid and rmsf_a are different lengths")
        for key in ("rmsd_2d", "cross_correlation"):
            matrix = panels.get(key)
            if matrix is None:
                continue
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                problems.append(f"panels/{cid}.npz {key} is not square")
            elif matrix.shape[0] > MAX_MATRIX_SIDE:
                problems.append(
                    f"panels/{cid}.npz {key} is {matrix.shape[0]} on a side, "
                    f"over {MAX_MATRIX_SIDE}")

    barcode = package.arrays(f"barcode/{cid}.npz")
    if barcode is not None:
        z = barcode.get("z")
        labels = barcode.get("labels")
        if z is None or labels is None:
            problems.append(f"barcode/{cid}.npz needs both z and labels")
        else:
            rows, cols = z.shape
            if rows > MAX_BARCODE_ROWS:
                problems.append(f"barcode/{cid}.npz has {rows} rows, over {MAX_BARCODE_ROWS}")
            if cols > MAX_BARCODE_COLS:
                problems.append(f"barcode/{cid}.npz has {cols} columns, over {MAX_BARCODE_COLS}")
            if labels.shape[0] != rows:
                problems.append(
                    f"barcode/{cid}.npz has {labels.shape[0]} labels for {rows} rows")

    return problems


def _check_metrics(package: ReadPackage, cid: str) -> list[str]:
    """Checks 9, 10 and 13: the fractions, the window, and the required keys."""
    problems: list[str] = []
    metrics = package.metrics(cid)
    if metrics is None:
        return problems

    status = str(metrics.get("status") or "")
    graded = package.manifest.get("kind") == "graded"

    if graded and status == "ok":
        missing = [key for key in REQUIRED_METRIC_KEYS if metrics.get(key) is None]
        if missing:
            problems.append(f"metrics/{cid}.json is missing {', '.join(missing)}")

    for key in FRACTION_KEYS:
        value = metrics.get(key)
        if isinstance(value, (int, float)) and not 0.0 <= float(value) <= 1.0:
            problems.append(f"metrics/{cid}.json {key} is {value}, outside 0 to 1")

    inter = package.interactions(cid)
    for key in FRACTION_KEYS:
        value = (inter or {}).get(key)
        if isinstance(value, (int, float)) and not 0.0 <= float(value) <= 1.0:
            problems.append(f"metrics/{cid}.interactions.json {key} is {value}, outside 0 to 1")

    for entry in (package.summary(cid) or []):
        occupancy = entry.get("occupancy")
        if isinstance(occupancy, (int, float)) and not 0.0 <= float(occupancy) <= 1.0:
            problems.append(f"summary/{cid}.json has an occupancy of {occupancy}")

    start = metrics.get("equil_start_frame")
    frames = metrics.get("n_frames")
    if isinstance(start, int) and isinstance(frames, int):
        if not 0 <= start < frames:
            problems.append(
                f"metrics/{cid}.json equil_start_frame is {start} against "
                f"n_frames {frames}")

    return problems
