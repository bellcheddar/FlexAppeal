"""Reading and validating ``.fxa`` results files.

The ``.fxa`` is the interchange contract between the two tabs: the bundle writes
one on the user's machine, the Analysis tab reads it here. Its shape is fixed by
``runtime/analyse.py.j2`` and mirrored by ``FXA_VERSION`` below; when the two
disagree, this module is the one that has to say so clearly.

Everything here treats the file as hostile input. It arrives over an upload form
from an unauthenticated user, and a zip is a classic vehicle for path traversal
and decompression bombs, so extraction is bounded on every axis: entry count,
declared size, actual size, and resolved path. The checks run as a full pass over
the member list *before* anything is read, so a hostile archive cannot get its
safe-looking half extracted before the dangerous entry is noticed.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# The .fxa layout this module understands. runtime/analyse.py.j2 writes it.
FXA_VERSION = 1

# --- extraction limits -----------------------------------------------------
# A results file is a manifest, some JSON, an npz and one decimated trajectory.
# Even the "full" tier is a handful of members, so a generous cap is still
# orders of magnitude tighter than anything an attacker needs.
MAX_ENTRIES = 64
MAX_TOTAL_UNCOMPRESSED = 1024 * 1024 * 1024   # 1 GB
MAX_SINGLE_MEMBER = 512 * 1024 * 1024
# A ratio this high is not achievable by real trajectory data; it is the
# signature of a zip bomb padded with zeros.
MAX_COMPRESSION_RATIO = 200

REQUIRED_MEMBERS = ("manifest.json", "metrics.json")


class FxaError(ValueError):
    """A results file could not be read. The message is shown to the user."""


@dataclass
class Results:
    """A validated results file, ready to plot."""

    manifest: dict[str, Any]
    metrics: dict[str, Any]
    arrays: dict[str, np.ndarray] = field(default_factory=dict)
    members: dict[str, bytes] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    # --- convenience accessors, so templates never dig through the manifest ---

    @property
    def job_name(self) -> str:
        return str(self.manifest.get("job_name") or "unnamed run")

    @property
    def config(self) -> dict[str, Any]:
        return self.manifest.get("config") or {}

    @property
    def run(self) -> dict[str, Any]:
        return self.manifest.get("run") or {}

    @property
    def trajectory(self) -> dict[str, Any]:
        return self.manifest.get("trajectory") or {}

    @property
    def tier(self) -> str:
        return str(self.manifest.get("payload_tier") or "unknown")

    @property
    def duration_ns(self) -> float:
        return float(self.trajectory.get("duration_ns") or 0.0)

    @property
    def has_membrane(self) -> bool:
        return bool(self.config.get("use_membrane"))

    @property
    def has_ligands(self) -> bool:
        return bool(self.config.get("has_ligands"))

    @property
    def computed(self) -> list[str]:
        return list(self.manifest.get("metrics_computed") or [])

    def has(self, *keys: str) -> bool:
        """Whether every named metric is present and non-empty."""
        return all(self.metrics.get(k) not in (None, [], {}) for k in keys)

    def array(self, name: str) -> np.ndarray | None:
        return self.arrays.get(name)

    @property
    def topology_pdb(self) -> bytes | None:
        return self.members.get("topology.pdb")

    @property
    def trajectory_xtc(self) -> bytes | None:
        return self.members.get("traj.xtc")


# ===========================================================================
#  Safe extraction
# ===========================================================================


def _inspect(archive: zipfile.ZipFile) -> list[zipfile.ZipInfo]:
    """Vet every member before reading any of them.

    Two passes on purpose. Checking each entry as it is extracted means a
    hostile archive gets its harmless prefix written to disk before the check
    trips -- so the whole member list is validated first, and only then is
    anything read.
    """
    entries = archive.infolist()
    if len(entries) > MAX_ENTRIES:
        raise FxaError(
            f"this archive contains {len(entries)} entries; a results file has a "
            f"handful. It is probably not a .fxa."
        )

    total = 0
    for entry in entries:
        name = entry.filename
        if entry.is_dir():
            continue
        # Absolute paths, parent traversal, and Windows drive letters are all
        # ways of naming somewhere outside the destination.
        if name.startswith("/") or name.startswith("\\") or ".." in name.replace("\\", "/").split("/"):
            raise FxaError(f"this archive contains an unsafe path ({name!r}); refusing to read it.")
        if ":" in name:
            raise FxaError(f"this archive contains an unsafe path ({name!r}); refusing to read it.")
        if entry.file_size > MAX_SINGLE_MEMBER:
            raise FxaError(
                f"{name} expands to {entry.file_size / 1e6:.0f} MB, which is larger "
                f"than any results file should contain."
            )
        if entry.compress_size > 0:
            ratio = entry.file_size / entry.compress_size
            if ratio > MAX_COMPRESSION_RATIO and entry.file_size > 10 * 1024 * 1024:
                raise FxaError(
                    f"{name} claims a {ratio:.0f}x compression ratio, which real "
                    f"trajectory data does not achieve. Refusing to expand it."
                )
        total += entry.file_size

    if total > MAX_TOTAL_UNCOMPRESSED:
        raise FxaError(
            f"this archive expands to {total / 1e9:.1f} GB. Results files are tens "
            f"of megabytes; rebuild with a smaller payload tier."
        )
    return [e for e in entries if not e.is_dir()]


def _read_members(content: bytes) -> dict[str, bytes]:
    try:
        archive = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile:
        raise FxaError(
            "this file is not a valid .fxa. It should be the file your run produced "
            "in its output directory, not the trajectory or the bundle."
        ) from None

    with archive:
        entries = _inspect(archive)
        members: dict[str, bytes] = {}
        for entry in entries:
            try:
                data = archive.read(entry)
            except (zipfile.BadZipFile, RuntimeError, OSError) as exc:
                raise FxaError(f"{entry.filename} could not be read: {exc}") from None
            # The declared size is attacker-controlled; check what actually came out.
            if len(data) > MAX_SINGLE_MEMBER:
                raise FxaError(f"{entry.filename} is larger than its header claimed.")
            members[entry.filename] = data
    return members


# ===========================================================================
#  Loading
# ===========================================================================


def load(content: bytes, *, verify_checksums: bool = True) -> Results:
    """Parse and validate a results file. Raises FxaError with a usable message."""
    if not content:
        raise FxaError("that file is empty.")

    members = _read_members(content)
    warnings: list[str] = []

    missing = [name for name in REQUIRED_MEMBERS if name not in members]
    if missing:
        raise FxaError(
            f"this .fxa is missing {', '.join(missing)}. It may be truncated, or "
            f"produced by a different tool."
        )

    try:
        manifest = json.loads(members["manifest.json"])
        metrics = json.loads(members["metrics.json"])
    except json.JSONDecodeError as exc:
        raise FxaError(f"the results file's metadata is corrupt: {exc}") from None

    # json.loads("3") is an int and json.loads("[]") is a list -- both valid
    # JSON, neither a mapping, and .get on either is an AttributeError.
    if not isinstance(manifest, dict) or not isinstance(metrics, dict):
        raise FxaError("the results file's metadata is not in the expected form.")

    version = manifest.get("fxa_version")
    if version is None:
        warnings.append("this file declares no format version; reading it as version 1.")
    elif not isinstance(version, int):
        raise FxaError(f"the results file declares an invalid format version ({version!r}).")
    elif version > FXA_VERSION:
        raise FxaError(
            f"this results file is format version {version}, but this server "
            f"understands up to {FXA_VERSION}. The run was made with a newer "
            f"FlexAppeal than this site; rebuild the bundle here and re-run."
        )

    if verify_checksums:
        warnings.extend(_verify(manifest, members))

    arrays: dict[str, np.ndarray] = {}
    if "arrays.npz" in members:
        try:
            # allow_pickle stays off: an .npz is a zip of .npy files, and a
            # pickled object array in one is arbitrary code execution on load.
            with np.load(io.BytesIO(members["arrays.npz"]), allow_pickle=False) as data:
                arrays = {name: data[name] for name in data.files}
        except (ValueError, OSError, zipfile.BadZipFile) as exc:
            warnings.append(f"the bulk arrays could not be read ({exc}); "
                            f"heat-map panels will be unavailable.")

    if not manifest.get("metrics_computed"):
        warnings.append("this file records no metric list; showing whatever it contains.")

    return Results(manifest=manifest, metrics=metrics, arrays=arrays,
                   members=members, warnings=warnings)


def _verify(manifest: dict[str, Any], members: dict[str, bytes]) -> list[str]:
    """Check the recorded checksums.

    A mismatch is a warning rather than an error: the data is still readable and
    a truncated upload will fail more informatively downstream. What matters is
    that the user is told the file is not what its manifest says it is.
    """
    warnings: list[str] = []
    declared = manifest.get("files")
    if not isinstance(declared, dict):
        return ["this file records no checksums, so its integrity cannot be confirmed."]

    for name, info in declared.items():
        if not isinstance(info, dict) or name not in members:
            continue
        expected = info.get("sha256")
        if not expected:
            continue
        actual = hashlib.sha256(members[name]).hexdigest()
        if actual != expected:
            warnings.append(
                f"{name} does not match the checksum recorded when it was written; "
                f"the file may have been modified or corrupted in transfer."
            )
    return warnings


def summarise(results: Results) -> dict[str, Any]:
    """The headline facts, for the top of the Analysis page."""
    run = results.run
    trajectory = results.trajectory
    config = results.config

    return {
        "job_name": results.job_name,
        "description": config.get("job_description") or "",
        "author": config.get("job_author") or "",
        "created_at": results.manifest.get("created_at", ""),
        "duration_ns": results.duration_ns,
        "frames_analysed": trajectory.get("frames_analysed"),
        "frames_packed": trajectory.get("frames_packed"),
        "atoms_packed": trajectory.get("atoms_packed"),
        "tier": results.tier,
        "force_field": config.get("protein_ff", ""),
        "water_model": config.get("water_model", ""),
        "temperature": config.get("temperature"),
        "timestep": config.get("timestep"),
        "platform": run.get("platform", ""),
        "ns_per_day": run.get("achieved_ns_per_day"),
        "wall_hours": (round(run["wall_seconds"] / 3600, 2)
                       if isinstance(run.get("wall_seconds"), (int, float)) else None),
        "openmm_version": run.get("openmm_version", ""),
        "flexappeal_version": results.manifest.get("flexappeal_version", ""),
        "structure": results.manifest.get("structure", {}),
        "machine": run.get("machine", {}),
    }
