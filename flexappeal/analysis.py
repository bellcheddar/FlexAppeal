"""Bounded server-side re-analysis of an already-packed results file.

This is the one place FlexAppeal does science on the droplet, and it exists for
the case the bundle's standard panel did not anticipate: a metric that was not
requested at build time, or the same metric restricted to part of the molecule
("RMSD of just the active-site loop", "contacts within this domain").

Everything here is bounded, because the alternative is a shared 3.8 GB box
running someone's contact map over two thousand frames:

* only ``full``-tier results carry enough trajectory to be worth re-analysing
* an atom x frame budget rejects the request before any array is allocated
* one job at a time, enforced by an exclusive lock file
* a wall-clock timeout kills anything that gets past the estimate

The work runs in a detached subprocess rather than a request thread, following
AlphaFraud's pattern: MDTraj holds the GIL in C extensions and would block a
gunicorn worker for the duration, and a detached process survives a restart.
"""

from __future__ import annotations

import json
import os
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

# An atom x frame product above this is refused. 10 million is 2000 frames of
# 5,000 atoms -- comfortably more than any `full` tier packs for a typical
# protein.
#
# MUST stay consistent with MemoryMax in deploy/flexappeal-web.service. The
# coordinate array alone is 12 bytes per atom-frame (120 MB here), and MDTraj's
# atom_slice and superpose each take a copy, so the real ceiling is several
# times that. The number is chosen so the largest permitted job fits inside the
# 1200M cgroup limit rather than being OOM-killed halfway through.
BUDGET_ATOM_FRAMES = 10_000_000

# Pairwise work scales with the square of the residue count, so it gets its own
# ceiling rather than riding on the atom budget.
MAX_RESIDUE_PAIRS = 250_000

MAX_SELECTION_LENGTH = 200
TIMEOUT_SECONDS = 240

# MDTraj's selection grammar is a pyparsing DSL and rejects Python syntax
# outright -- verified against injection attempts, not assumed. These patterns
# are defence in depth on top of that, because the parser does ultimately
# compile what it parsed, and a cheap gate costs nothing.
#
# The allowed set covers everything the grammar uses: keywords, residue and atom
# names, comparison operators, numbers and quoted strings.
_SELECTION_ALLOWED = re.compile(r"^[A-Za-z0-9_ ()<>=!.,\-'\"]+$")

# And the two shapes that would matter if the parser ever loosened: a dunder
# (the front door to __import__ and friends) and attribute access. A decimal
# point between digits is fine -- `mass > 12.5` is a legitimate selection.
_SELECTION_FORBIDDEN = re.compile(r"__|\.[A-Za-z_]|[A-Za-z_]\.")

METRICS = {
    "rmsd": "RMSD against the reference",
    "rmsf": "Per-residue fluctuation",
    "rgyr": "Radius of gyration",
    "sasa": "Solvent-accessible surface area",
    "dssp": "Secondary structure",
    "hbonds": "Hydrogen bonds",
    "contacts": "Residue contact map",
    "pca": "Principal component analysis",
    "distance": "Distance between two selections",
}


class ReanalysisError(ValueError):
    """A re-analysis request cannot be honoured. The message reaches the user."""


def _reimage(traj) -> bool:
    """Undo periodic wrapping, anchoring on the largest molecule.

    MDTraj's own anchor heuristic cannot be used here. Left to choose for
    itself it looks for a molecule with *more* atoms than the largest one,
    which is unsatisfiable whenever the system is a protein plus one ligand:

        ValueError: Could not find any anchor molecules. Based on our
        heuristic, those should be molecules with more than 3220 atoms.

    Naming the anchor explicitly sidesteps that. Without this step a bound
    ligand sits a whole box length from its protein in wrapped coordinates --
    a 4 A contact reads as 65 A -- and because compute_contacts applies the
    minimum image convention internally, the contact map looks perfectly
    correct beside the wrong distances.
    """
    if traj.unitcell_vectors is None:
        return False
    try:
        molecules = traj.topology.find_molecules()
        if not molecules:
            return False
        traj.image_molecules(inplace=True, anchor_molecules=[max(molecules, key=len)])
        return True
    except (ValueError, RuntimeError, IndexError):
        return False


@dataclass
class Request:
    metrics: list[str] = field(default_factory=list)
    selection: str = "protein"
    selection_b: str = ""
    reference: str = "first"
    label: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _check_selection(expression: str, what: str = "selection") -> str:
    expression = (expression or "").strip()
    if not expression:
        raise ReanalysisError(f"the {what} is empty.")
    if len(expression) > MAX_SELECTION_LENGTH:
        raise ReanalysisError(
            f"the {what} is longer than {MAX_SELECTION_LENGTH} characters."
        )
    if not _SELECTION_ALLOWED.match(expression) or _SELECTION_FORBIDDEN.search(expression):
        raise ReanalysisError(
            f"the {what} contains characters MDTraj's selection language does "
            f"not use. Try something like 'protein and resid 40 to 60'."
        )
    return expression


def parse_request(raw: dict[str, Any]) -> Request:
    """Validate a submitted re-analysis request."""
    if not isinstance(raw, dict):
        raise ReanalysisError("the request was not in the expected form.")

    metrics = raw.get("metrics") or []
    if isinstance(metrics, str):
        metrics = [metrics]
    metrics = [str(m) for m in metrics]
    unknown = [m for m in metrics if m not in METRICS]
    if unknown:
        raise ReanalysisError(
            f"unknown metric(s): {', '.join(unknown)}. Available: "
            f"{', '.join(sorted(METRICS))}."
        )
    if not metrics:
        raise ReanalysisError("choose at least one metric to compute.")

    request = Request(
        metrics=metrics,
        selection=_check_selection(raw.get("selection", "protein")),
        reference=str(raw.get("reference", "first")),
        label=str(raw.get("label", ""))[:80],
    )
    if "distance" in metrics:
        request.selection_b = _check_selection(
            raw.get("selection_b", ""), "second selection")
    return request


def check_budget(n_atoms: int, n_frames: int, n_residues: int = 0,
                 metrics: list[str] | None = None) -> None:
    """Refuse work that is too large, before allocating anything."""
    product = n_atoms * n_frames
    if product > BUDGET_ATOM_FRAMES:
        raise ReanalysisError(
            f"this would process {product:,} atom-frames, above the "
            f"{BUDGET_ATOM_FRAMES:,} limit for server-side analysis. Narrow the "
            f"selection, or re-run the bundle locally with the metrics you want."
        )
    if metrics and n_residues:
        pairwise = {"contacts"} & set(metrics)
        if pairwise and (n_residues * n_residues) > MAX_RESIDUE_PAIRS:
            raise ReanalysisError(
                f"a contact map over {n_residues} residues is {n_residues ** 2:,} "
                f"pairs, above the {MAX_RESIDUE_PAIRS:,} limit. Restrict the "
                f"selection to the region you care about."
            )


# ===========================================================================
#  Concurrency
# ===========================================================================


class Busy(RuntimeError):
    """Another re-analysis is already running."""


def acquire_lock(root: Path, ttl: int = TIMEOUT_SECONDS * 2) -> Path:
    """Take the single re-analysis slot, or raise Busy.

    An exclusive-create lock file rather than a database or a semaphore: this
    is one droplet with one worker slot, and a stale lock from a killed process
    has to expire on its own, which a mtime check gives for free.
    """
    lock = root / ".reanalyse.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            age = time.time() - lock.stat().st_mtime
        except OSError:
            age = 0
        if age > ttl:
            lock.unlink(missing_ok=True)
            return acquire_lock(root, ttl)
        raise Busy(
            "another analysis is running. This server runs one at a time; "
            "try again in a minute."
        ) from None
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    return lock


def release_lock(root: Path) -> None:
    (root / ".reanalyse.lock").unlink(missing_ok=True)


# ===========================================================================
#  The work
# ===========================================================================


def run(fxa_path: Path, request: Request) -> dict[str, Any]:
    """Compute the requested metrics. Imports MDTraj only when actually used."""
    import mdtraj as md
    import numpy as np

    from . import fxa as fxa_module

    results = fxa_module.load(fxa_path.read_bytes(), verify_checksums=False)
    if not results.trajectory_xtc or not results.topology_pdb:
        raise ReanalysisError(
            "this results file carries no trajectory, so there is nothing to "
            "re-analyse. Rebuild the bundle with the 'full' payload tier."
        )

    # MDTraj needs real files; the scratch directory is already per-session.
    workdir = fxa_path.parent
    topology_path = workdir / "_reanalyse_top.pdb"
    trajectory_path = workdir / "_reanalyse_traj.xtc"
    topology_path.write_bytes(results.topology_pdb)
    trajectory_path.write_bytes(results.trajectory_xtc)

    try:
        traj = md.load(str(trajectory_path), top=str(topology_path))
    except Exception as exc:  # noqa: BLE001 -- MDTraj raises several types
        raise ReanalysisError(f"the trajectory could not be read: {exc}") from None
    finally:
        topology_path.unlink(missing_ok=True)
        trajectory_path.unlink(missing_ok=True)

    # Re-image before measuring anything. Trajectories are written with
    # enforcePeriodicBox, which wraps each molecule into the primary cell
    # independently -- so a bound ligand can sit a full box length from its
    # protein in the stored coordinates. Distances then come back as ~65 A for
    # a contact that is really 4 A. compute_contacts applies the minimum image
    # convention itself and is unaffected, which is precisely why this is easy
    # to miss: the contact map looks right while the distances are nonsense.
    _reimage(traj)

    try:
        indices = traj.topology.select(request.selection)
    except Exception:  # noqa: BLE001 -- the parser raises bare ValueError
        raise ReanalysisError(
            f"{request.selection!r} is not a valid MDTraj selection. Try "
            f"'protein', 'name CA', or 'resid 40 to 60'."
        ) from None

    if len(indices) == 0:
        raise ReanalysisError(
            f"{request.selection!r} matched no atoms in this trajectory."
        )

    residues = {traj.topology.atom(i).residue.index for i in indices}
    check_budget(len(indices), traj.n_frames, len(residues), request.metrics)

    subset = traj.atom_slice(indices)
    reference = subset[0]
    subset.superpose(reference)

    metrics: dict[str, Any] = {
        "time_ns": (subset.time / 1000.0).tolist(),
        "selection": request.selection,
        "n_atoms": int(len(indices)),
        "n_frames": int(subset.n_frames),
        "n_residues": len(residues),
        "label": request.label,
    }

    if "rmsd" in request.metrics:
        metrics["rmsd_nm"] = md.rmsd(subset, reference).tolist()

    if "rgyr" in request.metrics:
        metrics["rgyr_nm"] = md.compute_rg(subset).tolist()

    if "rmsf" in request.metrics:
        ca = subset.topology.select("name CA")
        target = ca if len(ca) else np.arange(subset.n_atoms)
        metrics["rmsf_nm"] = md.rmsf(subset, subset, 0, atom_indices=target).tolist()
        metrics["rmsf_resids"] = [subset.topology.atom(i).residue.resSeq for i in target]
        metrics["rmsf_resnames"] = [subset.topology.atom(i).residue.name for i in target]

    if "sasa" in request.metrics:
        sasa = md.shrake_rupley(subset, mode="residue")
        metrics["sasa_total_nm2"] = sasa.sum(axis=1).tolist()

    if "dssp" in request.metrics:
        try:
            dssp = md.compute_dssp(subset, simplified=True)
            metrics["dssp_helix_fraction"] = (dssp == "H").mean(axis=1).tolist()
            metrics["dssp_sheet_fraction"] = (dssp == "E").mean(axis=1).tolist()
        except Exception as exc:  # noqa: BLE001
            metrics["dssp_error"] = str(exc)

    if "hbonds" in request.metrics:
        try:
            bonds = md.baker_hubbard(subset, freq=0.1, periodic=False)
            metrics["hbond_count"] = int(len(bonds))
            metrics["hbonds"] = [
                {"donor": str(subset.topology.atom(d)),
                 "hydrogen": str(subset.topology.atom(h)),
                 "acceptor": str(subset.topology.atom(a))}
                for d, h, a in bonds
            ][:200]
        except Exception as exc:  # noqa: BLE001
            metrics["hbond_error"] = str(exc)

    if "contacts" in request.metrics:
        n = subset.topology.n_residues
        pairs = [(i, j) for i in range(n) for j in range(i + 3, n)]
        distances, _ = md.compute_contacts(subset, contacts=pairs,
                                           scheme="closest-heavy")
        occupancy = (distances < 0.45).mean(axis=0)
        matrix = np.zeros((n, n), dtype=float)
        for (i, j), value in zip(pairs, occupancy):
            matrix[i, j] = matrix[j, i] = value
        metrics["contact_map"] = matrix.round(4).tolist()
        metrics["contact_residues"] = [r.resSeq for r in subset.topology.residues]

    if "pca" in request.metrics:
        ca = subset.topology.select("name CA")
        target = ca if len(ca) >= 3 else np.arange(subset.n_atoms)
        coords = subset.atom_slice(target).xyz.reshape(subset.n_frames, -1)
        coords = coords - coords.mean(axis=0)
        _, singular, components = np.linalg.svd(coords, full_matrices=False)
        variance = singular ** 2 / max(subset.n_frames - 1, 1)
        metrics["pca_variance_ratio"] = (variance[:10] / variance.sum()).tolist()
        metrics["pca_projection"] = (coords @ components[:3].T).round(4).tolist()

    if "distance" in request.metrics and request.selection_b:
        try:
            other = traj.topology.select(request.selection_b)
        except Exception:  # noqa: BLE001
            raise ReanalysisError(
                f"{request.selection_b!r} is not a valid MDTraj selection."
            ) from None
        if len(other) == 0:
            raise ReanalysisError(
                f"{request.selection_b!r} matched no atoms in this trajectory."
            )
        # Centre-of-geometry separation, which is what "distance between these
        # two things" means when either side is more than one atom.
        a = traj.xyz[:, indices, :].mean(axis=1)
        b = traj.xyz[:, other, :].mean(axis=1)
        metrics["distance_nm"] = np.linalg.norm(a - b, axis=1).tolist()
        metrics["distance_selection_b"] = request.selection_b

    return metrics


def run_to_file(fxa_path: Path, request_path: Path, output_path: Path) -> int:
    """Entry point for the detached worker. Never raises; writes status out."""
    try:
        raw = json.loads(request_path.read_text())
        request = parse_request(raw)
        metrics = run(fxa_path, request)
        output_path.write_text(json.dumps({"status": "ready", "metrics": metrics}))
        return 0
    except ReanalysisError as exc:
        output_path.write_text(json.dumps({"status": "error", "message": str(exc)}))
        return 1
    except Exception as exc:  # noqa: BLE001 -- the worker must never die silently
        output_path.write_text(json.dumps({
            "status": "error",
            "message": f"the analysis failed unexpectedly: {type(exc).__name__}: {exc}",
        }))
        return 1
    finally:
        release_lock(fxa_path.parent.parent)
