"""Assemble a self-contained run bundle from a validated configuration.

The output is one file: a bash self-extractor with a base64 tar payload. It
carries the prepared structure, a pinned environment specification, a readable
OpenMM script, the analysis script, and the installer that ties them together.

Two decisions worth stating, because both could reasonably have gone the other way:

**The generated ``run.py`` inlines every value** rather than reading config.json
at run time. A user should be able to open it and see exactly what will happen --
``310.0 * unit.kelvin``, not ``cfg["temperature"]``. config.json travels alongside
as machine-readable provenance and is copied into the results manifest.

**The payload is base64 rather than appended binary.** It costs 33% in size, and
buys a file that survives being emailed, pasted, or served by something that
decides to be helpful about line endings. A truncated binary payload fails in
confusing ways; a truncated base64 one fails immediately and says so.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import io
import json
import re
import tarfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from . import options as opts
from . import schema
from .options import FLEXAPPEAL_VERSION

RUNTIME_DIR = Path(__file__).resolve().parent / "runtime"

# Pin the bundle's OpenMM to the minor series FlexAppeal's tests run against.
# A generated script is validated against exactly one API surface; letting the
# bundle solve a newer major version would silently invalidate that.
OPENMM_PIN = "8.5.*"

SITE_URL = "https://flexappeal.mdeller.com"

# Implicit-solvent XML file -> the openmm.app symbol createSystem wants.
_IMPLICIT_SYMBOL = {
    "implicit/gbn2.xml": "GBn2",
    "implicit/gbn.xml": "GBn",
    "implicit/obc1.xml": "OBC1",
    "implicit/obc2.xml": "OBC2",
    "implicit/hct.xml": "HCT",
}

_TRAJ_EXTENSION = {"xtc": "xtc", "dcd": "dcd", "hdf5": "h5"}

_MUTATION_RE = re.compile(r"^([A-Za-z0-9]+):([A-Z]{3})-(\d+)-([A-Z]{3})$")


class BundleError(ValueError):
    """A bundle could not be built from this configuration."""


@dataclass(frozen=True)
class Bundle:
    filename: str
    content: bytes
    manifest: dict[str, Any]

    @property
    def size_human(self) -> str:
        n = len(self.content)
        for unit_name in ("B", "KB", "MB"):
            if n < 1024:
                return f"{n:.0f} {unit_name}"
            n /= 1024
        return f"{n:.1f} GB"


# ===========================================================================
#  Template context
# ===========================================================================


def _forcefield_files(cfg: dict[str, Any]) -> list[str]:
    """The XML files the generated script loads, in load order."""
    files = [cfg["protein_ff"]]

    if cfg["solvent_mode"] == "explicit":
        water_xml, _ = opts.WATER_MODEL_XML[cfg["water_model"]]
        files.append(water_xml)
        if cfg.get("use_membrane"):
            # addMembrane needs lipid parameters. CHARMM36 carries its own; an
            # AMBER protein force field needs lipid17 loaded alongside.
            if cfg["protein_ff"] != "charmm36.xml":
                files.append("amber14/lipid17.xml")
    elif cfg["solvent_mode"] == "implicit":
        files.append(cfg["implicit_model"])

    return files


def _mutations(cfg: dict[str, Any]) -> dict[str, list[str]]:
    """Parse the mutation textarea into PDBFixer's per-chain form."""
    grouped: dict[str, list[str]] = {}
    for line in str(cfg.get("mutations") or "").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = _MUTATION_RE.match(line)
        if not match:
            raise BundleError(
                f"cannot parse the mutation {line!r}. Use CHAIN:WT-RESID-MUT with "
                f"three-letter codes, for example A:ALA-57-GLY."
            )
        chain, wt, resid, mut = match.groups()
        grouped.setdefault(chain, []).append(f"{wt}-{resid}-{mut}")
    return grouped


def _restraint_schedule(cfg: dict[str, Any]) -> list[float]:
    raw = str(cfg.get("restraint_schedule") or "")
    values: list[float] = []
    for part in raw.replace(",", " ").split():
        try:
            values.append(float(part))
        except ValueError:
            raise BundleError(
                f"the restraint release schedule contains {part!r}, which is not a "
                f"number. Give a comma-separated list of force constants."
            ) from None
    return values or [float(cfg.get("restraint_force") or 1000.0), 0.0]


def _vector3(raw: str, what: str) -> list[float]:
    parts = [p for p in str(raw).replace(",", " ").split() if p]
    if len(parts) != 3:
        raise BundleError(f"{what} needs exactly three numbers, got {len(parts)}.")
    try:
        return [float(p) for p in parts]
    except ValueError:
        raise BundleError(f"{what} must be three numbers.") from None


def _solvent_summary(cfg: dict[str, Any]) -> str:
    if cfg["solvent_mode"] == "vacuum":
        return "vacuum"
    if cfg["solvent_mode"] == "implicit":
        model = cfg["implicit_model"].split("/")[-1].replace(".xml", "").upper()
        return f"implicit ({model}, {cfg['implicit_salt']} M salt)"
    water = cfg["water_model"]
    if cfg.get("use_membrane"):
        return (f"{cfg['lipid_type']} bilayer in {water}, "
                f"{cfg['ionic_strength']} M {cfg['positive_ion']}{cfg['negative_ion']}")
    return (f"{water} in a {cfg['box_shape']}, {cfg['padding']} nm padding, "
            f"{cfg['ionic_strength']} M {cfg['positive_ion']}{cfg['negative_ion']}")


def build_context(cfg: dict[str, Any], structure_filename: str,
                  citation: str = "") -> dict[str, Any]:
    """Everything the templates need, derived once so they stay declarative."""
    derived = schema.derive(cfg)
    wall = schema.estimate_wall_time(cfg)

    needs_pdbfixer = bool(
        cfg.get("fix_missing_residues") or cfg.get("replace_nonstandard")
        or cfg.get("add_missing_atoms") or cfg.get("strip_hydrogens")
        or cfg.get("mutations") or cfg.get("chains")
    )

    ff_files = _forcefield_files(cfg)
    # ff19SB ships with openmmforcefields rather than OpenMM itself, so it pulls
    # that dependency into the bundle even when there is no ligand at all.
    ff19sb = cfg["protein_ff"].startswith("amber/")
    needs_ommff = bool(cfg.get("has_ligands")) or ff19sb
    if ff19sb and not cfg.get("has_ligands"):
        reason = "ff19SB ships with openmmforcefields, not OpenMM itself."
    elif cfg.get("has_ligands"):
        reason = f"ligand parameters ({cfg.get('ligand_ff')})."
    else:
        reason = ""

    traj_selection = cfg["traj_selection"]
    if traj_selection == "custom":
        # MDTraj selection strings cannot be honoured by the OpenMM-only run
        # script, so a custom selection saves everything and the analysis step
        # applies the expression. Stated rather than silently downgraded.
        traj_selection = "all"

    mutations = _mutations(cfg)

    traj_gb = derived["traj_bytes"] / 1e9

    return {
        "cfg": cfg,
        "derived": derived,
        "version": FLEXAPPEAL_VERSION,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "structure_filename": structure_filename,
        "structure_display": structure_filename,
        "citation": citation,
        "site_url": SITE_URL,
        "bundle_filename": bundle_filename(cfg),
        "ff_files": ff_files,
        "water_geometry": (
            opts.WATER_MODEL_XML[cfg["water_model"]][1]
            if cfg["solvent_mode"] == "explicit" else "tip3p"
        ),
        "implicit_symbol": _IMPLICIT_SYMBOL.get(cfg.get("implicit_model", ""), "GBn2"),
        "needs_pdbfixer": needs_pdbfixer,
        "needs_openmmforcefields": needs_ommff,
        "openmmforcefields_reason": reason,
        "openmm_pin": OPENMM_PIN,
        "allow_linux": False,
        "keep_chains": list(cfg.get("chains") or []),
        "drop_chains": bool(cfg.get("chains")) and "*" not in (cfg.get("chains") or []),
        "mutations": mutations,
        "mutation_count": sum(len(v) for v in mutations.values()),
        "restraint_schedule": _restraint_schedule(cfg),
        "box_vectors": (_vector3(cfg.get("box_vectors", ""), "the box vectors")
                        if cfg.get("box_sizing") == "vectors" else [0, 0, 0]),
        "anisotropic_pressure": (
            _vector3(cfg.get("anisotropic_pressure", ""), "the per-axis pressure")
            if cfg.get("barostat") == "MonteCarloAnisotropic" else [1.0, 1.0, 1.0]
        ),
        "traj_selection": traj_selection,
        "traj_extension": _TRAJ_EXTENSION.get(cfg["traj_format"], "xtc"),
        "solvent_summary": _solvent_summary(cfg),
        "summary_line": (
            f"{cfg['production_duration']} ns · {cfg['protein_ff'].split('/')[-1]} · "
            f"{_solvent_summary(cfg)}"
        ),
        "wall_estimate": wall["human"],
        "traj_gb": f"{traj_gb:.1f}",
        # Environment (~3 GB) plus the trajectory, rounded up, minimum 5 GB.
        "required_gb": max(5, int(traj_gb + 4)),
    }


# ===========================================================================
#  Rendering
# ===========================================================================


def _environment() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(RUNTIME_DIR)),
        # StrictUndefined so a template referencing a key the context does not
        # provide fails here, loudly, rather than rendering an empty string into
        # someone's simulation script.
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )


def render(template_name: str, context: dict[str, Any]) -> str:
    return _environment().get_template(template_name).render(**context)


def bundle_filename(cfg: dict[str, Any]) -> str:
    return f"flexappeal_{cfg['job_name']}.command"


# ===========================================================================
#  Assembly
# ===========================================================================


def build(cfg: dict[str, Any], structure_bytes: bytes,
          structure_filename: str = "input.pdb", citation: str = "",
          ligands: dict[str, bytes] | None = None) -> Bundle:
    """Build the complete self-extracting bundle.

    Raises BundleError for anything wrong with the configuration; the caller is
    expected to have run schema.validate first, so anything reaching here is a
    parse failure in a free-text field rather than a policy violation.
    """
    result = schema.validate(cfg)
    if not result.ok:
        raise BundleError(
            "this configuration does not validate: "
            + "; ".join(f"{i.option_id}: {i.message}" for i in result.errors)
        )
    # validate() rebuilds the config from the registry's defaults, which drops
    # the internal underscore keys the structure analysis put there. Without
    # carrying them across, every size and wall-time estimate in the generated
    # README reads "unknown" -- they are inputs to schema.derive, not options.
    validated = result.config
    for key, value in cfg.items():
        if key.startswith("_"):
            validated[key] = value
    cfg = validated

    # The structure is embedded under a fixed name so run.py can reference it
    # literally, regardless of what the user called their upload.
    embedded_name = "input" + (".cif" if structure_filename.lower().endswith(
        (".cif", ".mmcif")) else ".pdb")

    context = build_context(cfg, embedded_name, citation)
    # Code refers to the embedded name; prose should name the file the user knows.
    context["structure_display"] = structure_filename

    files: dict[str, bytes] = {
        "run.py": render("run.py.j2", context).encode(),
        "analyse.py": render("analyse.py.j2", context).encode(),
        "pixi.toml": render("pixi.toml.j2", context).encode(),
        "README.md": render("README.md.j2", context).encode(),
        embedded_name: structure_bytes,
        "config.json": json.dumps({
            "flexappeal_version": FLEXAPPEAL_VERSION,
            "generated_at": context["generated_at"],
            "structure": {"filename": structure_filename, "citation": citation,
                          "sha256": hashlib.sha256(structure_bytes).hexdigest()},
            "config": {k: v for k, v in cfg.items() if not k.startswith("_")},
            "derived": context["derived"],
        }, indent=2).encode(),
    }
    for name, payload in (ligands or {}).items():
        files[f"ligands/{name}"] = payload

    payload_b64 = _pack(files)
    header = render("bootstrap.sh.j2", context)
    content = header.encode() + payload_b64

    manifest = {
        "filename": context["bundle_filename"],
        "job_name": cfg["job_name"],
        "files": sorted(files),
        "payload_sha256": hashlib.sha256(payload_b64).hexdigest(),
        "bundle_bytes": len(content),
        "generated_at": context["generated_at"],
        "flexappeal_version": FLEXAPPEAL_VERSION,
    }
    return Bundle(context["bundle_filename"], content, manifest)


def _pack(files: dict[str, bytes]) -> bytes:
    """Tar, gzip and base64 the payload, deterministically.

    Fixed mtimes and uid/gid so the same configuration always produces the same
    bytes -- which is what makes the golden-file tests meaningful and lets a user
    verify two bundles are identical.
    """
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w") as tar:
        for name in sorted(files):
            data = files[name]
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            entry.mtime = 0
            entry.uid = entry.gid = 0
            entry.uname = entry.gname = ""
            entry.mode = 0o755 if name.endswith(".py") else 0o644
            tar.addfile(entry, io.BytesIO(data))

    compressed = gzip.compress(raw.getvalue(), compresslevel=9, mtime=0)
    encoded = base64.b64encode(compressed)

    # Wrapped at 76 columns: an unwrapped multi-megabyte line breaks editors,
    # some mail transports, and `awk` implementations with a line-length limit.
    return b"\n".join(
        encoded[i:i + 76] for i in range(0, len(encoded), 76)
    ) + b"\n"


def unpack(content: bytes) -> dict[str, bytes]:
    """Extract a bundle's payload. Used by the tests and by `FlexAppeal.py unpack`."""
    marker = b"\n__FLEXAPPEAL_PAYLOAD__\n"
    index = content.find(marker)
    if index < 0:
        raise BundleError("this does not look like a FlexAppeal bundle: no payload marker.")
    payload = content[index + len(marker):]
    try:
        decoded = base64.b64decode(b"".join(payload.split()))
        raw = gzip.decompress(decoded)
    except (ValueError, OSError) as exc:
        raise BundleError(f"the payload is corrupt or truncated: {exc}") from None

    files: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r") as tar:
        for entry in tar.getmembers():
            if not entry.isfile():
                continue
            handle = tar.extractfile(entry)
            if handle is not None:
                files[entry.name] = handle.read()
    return files
