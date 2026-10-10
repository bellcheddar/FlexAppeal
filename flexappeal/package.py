"""`flexappeal package`: collect finished runs into one `.fxa` version 2 file.

FlexAppeal runs **one system per job**, and a version 2 file needs **many
complexes**, because the league table, both rank columns and the energy sub-score
are all comparative: `S_E` is a percentile within a receptor's own cohort, so one
complex has no cohort and the sub-score is undefined rather than zero.

Rather than change how a job runs, packaging is a separate, deliberate act:

    flexappeal package --out study.fxa run-dir-1 run-dir-2 ...

That keeps the job model untouched, nothing about an existing run changes, no
version 1 file is altered, and a study can be assembled from runs done weeks
apart.

**Two source kinds**, detected per argument:

* a **finished run directory**, which carries the complete trajectory and the
  job's configuration. This is the better source and the only one from which the
  ligand strain can be measured, because that needs the chemical definition the
  bundle fetched.
* a **version 1 `.fxa`**, which carries whatever its payload tier packed. A
  `light` tier is 100 frames, which is enough for the pose and occupancy numbers
  and thin for the convergence ones, and the command says so.

A directory holding several of either is expanded, so a whole output tree can be
handed over at once.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import fxa2, measure
from .options import FLEXAPPEAL_VERSION

# Frames used for the measurements, at most. Every scalar here is a statistic
# over the scored window, and a statistic over 2,000 frames and one over 10,000
# agree to more digits than anybody reads; what the extra frames cost is real,
# because three superpositions each hold a copy of the coordinates.
#
# FlexAppeal's default run is 100 ns at 10 ps, which is 10,000 frames, so this
# cap is reached by an ordinary job rather than by an unusual one. At 3,000 atoms
# that was a gigabyte for 2,000 frames and would have been five for 10,000; on a
# membrane system, where `protein_ligand` keeps the lipids too, it was tens of
# gigabytes and would simply have failed.
#
# 2,000 is also the length of DeepFlex's own production runs, so a package built
# at the cap carries the frame count its report was designed around.
MAX_MEASURE_FRAMES = 2000


# Trajectory selections that keep no ligand. A run written with one of these
# cannot become a receptor-ligand complex however well it ran, and the failure
# arrives as "no ligand in the trajectory", which is true and says nothing about
# the cause. `protein_ligand` is the default and keeps the ligand.
LIGANDLESS_SELECTIONS = frozenset({"protein", "protein_heavy", "backbone", "ca"})


# How a version 1 payload tier maps onto a warning. `full` keeps every atom and
# 2,000 frames, which is as good as a run directory for everything but the
# chemistry; the lighter tiers are decimated and the numbers say so.
TIER_NOTE = {
    "light": "100 frames: enough for the pose and occupancy numbers, thin for "
             "the convergence block split",
    "standard": "500 frames",
    "full": "2,000 frames",
}


class PackageError(ValueError):
    """A package could not be assembled. The message is shown to the user."""


@dataclass
class Source:
    """One finished run, loaded and ready to measure."""

    label: str
    origin: str
    trajectory: Any                       # an mdtraj.Trajectory
    config: dict[str, Any] = field(default_factory=dict)
    structure: dict[str, Any] = field(default_factory=dict)
    run_manifest: dict[str, Any] = field(default_factory=dict)
    frame_spacing_ps: float | None = None
    ligand_dir: Path | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def job_name(self) -> str:
        return str(self.config.get("job_name") or self.label)

    @property
    def ligand_resnames(self) -> list[str]:
        """The heteroatom residues this job chose to keep.

        This is the only compound identity a FlexAppeal job has: a PDB chemical
        component code, which is a name rather than a structure, so it travels in
        clear in `manifest.ligand_names` as section 2 allows.
        """
        return [str(n).strip().upper()
                for n in (self.config.get("keep_heteroatoms") or []) if str(n).strip()]

    def receptor(self) -> str:
        """A short readable receptor token, in the order each can be trusted.

        `config.pdb_id` is what the user asked for by accession and is therefore
        the most reliable. The prepared structure's filename is next: it is what
        was actually simulated, but it may be a local name rather than an
        accession. The job name is last, because it is free text.
        """
        filename = str(self.structure.get("filename") or "")
        stem = Path(filename).stem if filename else ""
        return fxa2.receptor_token(self.config.get("pdb_id"), stem, self.job_name)


# ===========================================================================
#  Loading sources
# ===========================================================================


def _load_config(start: Path) -> tuple[dict, dict]:
    """A job's `config.json`, found by walking up from the run directory.

    A bundle unpacks to `<bundle>/config.json` with the run directories under
    `<bundle>/flexappeal_output/<job>/`, so the config is two levels above the
    trajectory. Three levels of search covers that and the case of somebody
    pointing at a run directory they moved.
    """
    for candidate in (start, *start.parents[:3]):
        path = candidate / "config.json"
        if path.is_file():
            try:
                document = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if isinstance(document, dict):
                return (document.get("config") or {}), (document.get("structure") or {})
    return {}, {}


def _frame_spacing_ps(config: dict, run_manifest: dict, n_frames: int) -> float | None:
    """Picoseconds between the frames that were actually written.

    `traj_interval` is the configured spacing, and it is the right answer unless
    the payload was decimated on the way into a version 1 file, in which case the
    caller multiplies by the stride. The report prints this and uses it to turn
    frame indices into nanoseconds, so a wrong value mislabels every x axis.
    """
    interval = config.get("traj_interval")
    if isinstance(interval, (int, float)) and interval > 0:
        return float(interval)
    # Fall back to the production length divided by the frames, which is right
    # whenever the run completed.
    duration = config.get("production_duration")
    if isinstance(duration, (int, float)) and duration > 0 and n_frames > 1:
        return float(duration) * 1000.0 / n_frames
    steps = run_manifest.get("production_steps")
    timestep = run_manifest.get("timestep_ps")
    if all(isinstance(v, (int, float)) for v in (steps, timestep)) and n_frames > 1:
        return float(steps) * float(timestep) / n_frames
    return None


def _reimage(traj, label: str, notes: list[str]) -> None:
    """Undo periodic wrapping, in place, where there is a box to do it with.

    A trajectory written with `enforcePeriodicBox` wraps every molecule into the
    primary cell independently, so a bound ligand can sit a whole box length from
    its receptor in the stored coordinates. Every distance that spans the two is
    then wrong by a box length, which is exactly the quantity this module
    measures. The anchor has to be named: MDTraj's own heuristic looks for a
    molecule with *more* atoms than the largest one, which can never be satisfied
    for a receptor plus a ligand.
    """
    if traj.unitcell_vectors is None:
        return
    try:
        molecules = traj.topology.find_molecules()
        if molecules:
            traj.image_molecules(inplace=True,
                                 anchor_molecules=[max(molecules, key=len)])
    except (ValueError, RuntimeError, IndexError) as exc:
        notes.append(f"{label}: could not re-image molecules ({exc}); distances "
                     f"that span two molecules may be wrong by a box length")


def measurement_selection(topology, resnames: list[str]) -> Any:
    """The atoms the measurements need: the protein and the kept heteroatoms.

    Everything else is dropped before the coordinates are read. That is not only
    a memory decision, although it is mostly that: section 6.1 bans waters,
    counter-ions and lipids from the pocket crop, and a 10 A shell drawn around a
    ligand in a solvated bilayer is mostly bilayer. Nothing measured here reads a
    solvent atom, so there is no reason to carry one.

    Returns None when the selection would match nothing, which means the
    topology is not shaped the way this expects and the caller should keep
    everything rather than measure an empty trajectory.
    """
    import numpy as np  # noqa: PLC0415

    wanted = ["protein"] + [f"resname {name}" for name in resnames]
    indices = topology.select(" or ".join(wanted))
    if len(indices) == 0 or len(indices) == topology.n_atoms:
        return None
    return np.asarray(indices)


def measurement_stride(n_frames: int) -> int:
    """Frames to skip so the measurement stays inside `MAX_MEASURE_FRAMES`."""
    if n_frames <= MAX_MEASURE_FRAMES:
        return 1
    return int(-(-n_frames // MAX_MEASURE_FRAMES))


def _frame_count(path: Path) -> int | None:
    """How many frames a trajectory holds, without reading its coordinates."""
    import mdtraj as md  # noqa: PLC0415

    try:
        with md.open(str(path)) as handle:
            return len(handle)
    except (OSError, ValueError, TypeError, NotImplementedError):
        return None


def load_run_directory(path: Path) -> Source:
    """A finished run directory: `topology.pdb` plus a `trajectory.*`.

    Reduced on the way in rather than afterwards. A default FlexAppeal job is
    100 ns at 10 ps, which is 10,000 frames of whatever `traj_selection` kept, so
    loading it whole and then throwing most of it away is the difference between
    a gigabyte and a job that does not finish.
    """
    import mdtraj as md  # noqa: PLC0415

    topology = path / "topology.pdb"
    candidates = sorted(path.glob("trajectory.*"))
    if not topology.is_file() or not candidates:
        raise PackageError(f"{path} is not a finished run directory")

    config, structure = _load_config(path)
    run_manifest: dict[str, Any] = {}
    manifest_path = path / "run_manifest.json"
    if manifest_path.is_file():
        try:
            run_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            run_manifest = {}

    notes: list[str] = []
    resnames = [str(n).strip().upper()
                for n in (config.get("keep_heteroatoms") or []) if str(n).strip()]
    atom_indices = measurement_selection(md.load_topology(str(topology)), resnames)

    stride = measurement_stride(_frame_count(candidates[0]) or 0)
    if stride > 1:
        notes.append(
            f"{path.name}: read every {stride} frames, which keeps the "
            f"measurement inside {MAX_MEASURE_FRAMES} frames. Every scalar here "
            f"is a statistic over the scored window and is unchanged by it; the "
            f"frame spacing in the manifest describes the frames actually written")

    traj = md.load(str(candidates[0]), top=str(topology),
                   atom_indices=atom_indices, stride=stride)
    _reimage(traj, path.name, notes)

    # The bundle fetched a chemical definition for each heteroatom it kept, and it
    # is still beside the run. Needed for the MMFF94 strain and for nothing else,
    # and it never enters the package: section 8 is absolute about that.
    ligand_dir = None
    for candidate in (path, *path.parents[:3]):
        if (candidate / "ligands").is_dir():
            ligand_dir = candidate / "ligands"
            break

    # The spacing of the frames this package carries, which is the configured
    # spacing times whatever was skipped. The report prints it and uses it to
    # turn a frame index into nanoseconds, so a stride left out of it mislabels
    # every x axis by that factor.
    spacing = _frame_spacing_ps(config, run_manifest, traj.n_frames * stride)
    if spacing is not None:
        spacing *= stride

    return Source(
        label=path.name, origin="run directory", trajectory=traj,
        config=config, structure=structure, run_manifest=run_manifest,
        frame_spacing_ps=spacing, ligand_dir=ligand_dir, notes=notes,
    )


def load_fxa1(path: Path) -> Source:
    """A version 1 `.fxa`, read with FlexAppeal's own reader.

    Its topology and trajectory are written out to a temporary directory because
    MDTraj loads from paths rather than from bytes. They are deleted as soon as
    the trajectory is in memory.
    """
    import tempfile  # noqa: PLC0415

    import mdtraj as md  # noqa: PLC0415

    from . import fxa as fxa1  # noqa: PLC0415

    try:
        results = fxa1.load(path.read_bytes())
    except fxa1.FxaError as exc:
        raise PackageError(f"{path.name}: {exc}") from None

    if results.topology_pdb is None or results.trajectory_xtc is None:
        raise PackageError(
            f"{path.name} carries no coordinates, so nothing can be measured from "
            f"it. Point at the run directory instead."
        )

    with tempfile.TemporaryDirectory(prefix="flexappeal-fxa1-") as work:
        work = Path(work)
        (work / "topology.pdb").write_bytes(results.topology_pdb)
        (work / "traj.xtc").write_bytes(results.trajectory_xtc)
        traj = md.load(str(work / "traj.xtc"), top=str(work / "topology.pdb"))

    notes = [f"{path.name}: read from a version 1 .fxa at the "
             f"{results.tier} payload tier ({TIER_NOTE.get(results.tier, 'decimated')})"]
    _reimage(traj, path.name, notes)

    stride = results.trajectory.get("stride") or 1
    spacing = _frame_spacing_ps(results.config, results.run, traj.n_frames)
    if spacing is not None and isinstance(stride, int) and stride > 1:
        # The configured interval describes the frames the run wrote; the file
        # carries every nth of them.
        spacing *= stride

    # The same reduction as a run directory, after the fact rather than before:
    # the payload is already in memory and a decimated one is small. It still
    # matters at the `full` tier, which keeps every atom including the solvent
    # and the bilayer.
    resnames = [str(n).strip().upper()
                for n in (results.config.get("keep_heteroatoms") or [])
                if str(n).strip()]
    atom_indices = measurement_selection(traj.topology, resnames)
    if atom_indices is not None:
        traj = traj.atom_slice(atom_indices)
    measure_stride = measurement_stride(traj.n_frames)
    if measure_stride > 1:
        traj = traj[::measure_stride]
        if spacing is not None:
            spacing *= measure_stride
        notes.append(
            f"{path.name}: read every {measure_stride} frames, which keeps the "
            f"measurement inside {MAX_MEASURE_FRAMES} frames")

    return Source(
        label=path.stem, origin="version 1 .fxa", trajectory=traj,
        config=results.config, structure=results.manifest.get("structure") or {},
        run_manifest=results.run, frame_spacing_ps=spacing, notes=notes,
    )


def discover(paths: list[Path]) -> list[Path]:
    """Expand each argument into the sources it names.

    A run directory or a `.fxa` is taken as itself. Any other directory is
    searched for both, one level into `flexappeal_output` as well, so a whole
    output tree can be handed over at once.
    """
    found: list[Path] = []
    for path in paths:
        path = Path(path)
        if not path.exists():
            raise PackageError(f"{path} does not exist")
        if path.is_file():
            found.append(path)
            continue
        if (path / "topology.pdb").is_file() and any(path.glob("trajectory.*")):
            found.append(path)
            continue
        here: list[Path] = sorted(path.glob("*.fxa"))
        for parent in (path, path / "flexappeal_output"):
            if not parent.is_dir():
                continue
            here.extend(sorted(
                child for child in parent.iterdir()
                if child.is_dir() and (child / "topology.pdb").is_file()
                and any(child.glob("trajectory.*"))
            ))
        if not here:
            raise PackageError(
                f"{path} holds no finished runs and no .fxa files. A source is "
                f"either a run directory (topology.pdb plus a trajectory) or a "
                f"version 1 .fxa."
            )
        found.extend(here)
    return found


def load(path: Path) -> Source:
    path = Path(path)
    if path.is_dir():
        return load_run_directory(path)
    if path.suffix.lower() == ".fxa":
        return load_fxa1(path)
    raise PackageError(f"{path.name} is not a run directory or a .fxa")


# ===========================================================================
#  Assembly
# ===========================================================================


def _run_settings(sources: list[Source]) -> dict[str, Any]:
    """The run settings for the manifest, from the sources that agree on them.

    One manifest describes a whole package, so a study assembled from runs with
    different settings has to report something. The first source's settings are
    reported and any disagreement becomes a warning: silently reporting one run's
    temperature for twenty is how a reader ends up trusting a number nothing in
    the file supports.
    """
    first = sources[0]
    config = first.config
    settings = {
        "protein_ff": config.get("protein_ff"),
        "ligand_ff": config.get("ligand_ff"),
        "water_model": config.get("water_model"),
        "temperature": config.get("temperature"),
        "pressure": config.get("pressure"),
        "timestep_fs": config.get("timestep"),
        "n_strides": 1,
        "stride_time_ns": config.get("production_duration"),
        "prod_traj_interval_ps": first.frame_spacing_ps,
        "contact_cutoff": measure.CONTACT_CUTOFF_A,
    }
    return {key: value for key, value in settings.items() if value is not None}


def _settings_disagreements(sources: list[Source]) -> list[str]:
    watched = ("protein_ff", "water_model", "temperature", "timestep",
               "production_duration")
    out: list[str] = []
    first = sources[0]
    for key in watched:
        values = {json.dumps(source.config.get(key)) for source in sources}
        if len(values) > 1:
            out.append(
                f"the sources disagree on {key}; the manifest reports "
                f"{first.config.get(key)!r} from {first.label} and the report will "
                f"print that for the whole package"
            )
    return out


def annotate(package: fxa2.Package, *, progress=None) -> list[str]:
    """Add `annotations/<receptor>.json` for every receptor, over the network.

    The only part of this format that needs the internet, and the only reason to
    make it opt-in: a packaging run should work on a train. A failure costs the
    receptor page its accession and nothing else, so each one is caught
    separately and reported rather than aborting the package.

    No consensus numbering. `resid_to_generic` is a GPCR convention built from an
    alignment to a reference receptor, and FlexAppeal has neither the reference
    nor the alignment; `n_generic_mapped: 0` says so rather than leaving the
    field out and implying the mapping was simply not attempted.
    """
    from . import sources as sources_mod  # noqa: PLC0415

    warnings: list[str] = []
    for receptor in sorted(package.receptors):
        if progress is not None:
            progress(f"resolving {receptor} against UniProt")
        try:
            uniprot = sources_mod.fetch_uniprot_for_pdb(receptor)
        except sources_mod.SourceError as exc:
            warnings.append(
                f"{receptor}: no UniProt annotation ({exc}). The receptor page "
                f"will show its identifier and no accession.")
            continue
        package.add_json(f"annotations/{receptor}.json", {
            "label": receptor,
            "uniprot": uniprot,
            "n_generic_mapped": 0,
            "notes": ["resolved from the RCSB entry's first polymer entity; no "
                      "consensus numbering was computed"],
        })
        package.register_receptor(
            receptor, protein_name=uniprot.get("protein_name"))
    return warnings


def assemble(sources: list[Source], *, package_name: str,
             energies: dict[str, dict[str, float]] | None = None,
             progress=None) -> tuple[fxa2.Package, list[str]]:
    """Measure every source and build the package. Returns it and any warnings."""
    if not sources:
        raise PackageError("no sources to package")

    warnings: list[str] = []
    if len(sources) < 2:
        warnings.append(
            "this package carries one complex. The league table, both rank columns "
            "and the energy sub-score are percentiles within a receptor's own "
            "cohort, so a single complex has no cohort and those panels stay empty. "
            "A version 1 .fxa is the complete, supported path for one system."
        )

    engine = next((str(source.run_manifest.get("openmm_version") or "")
                   for source in sources if source.run_manifest.get("openmm_version")), "")
    package = fxa2.Package(
        package_name=package_name,
        flexappeal_version=FLEXAPPEAL_VERSION,
        engine=f"openmm {engine}" if engine else "openmm",
        run_settings=_run_settings(sources),
    )
    warnings.extend(_settings_disagreements(sources))

    seen: dict[str, str] = {}
    with_energies: list[str] = []
    panel_lines: list[dict[str, Any]] = []
    conformer_lines: list[dict[str, Any]] = []
    painted: list[str] = []
    for source in sources:
        warnings.extend(source.notes)
        receptor = source.receptor()
        resnames = source.ligand_resnames
        if not resnames:
            warnings.append(
                f"{source.label}: the job kept no heteroatoms, so it has no ligand "
                f"and no complex. An apo run belongs in a version 1 .fxa.")
            continue

        # A compound per kept heteroatom would be wrong: a job simulates one
        # complex, and several kept residues are one ligand plus its cofactors.
        # The first is taken as the ligand and the rest are noted.
        compound = resnames[0]
        if len(resnames) > 1:
            warnings.append(
                f"{source.label}: kept {', '.join(resnames)}; {compound} is taken "
                f"as the ligand and the rest are treated as part of the receptor")

        ligand = package.register_ligand(compound)
        cid = fxa2.complex_id(ligand, receptor)
        if cid in seen:
            raise PackageError(
                f"{source.label} and {seen[cid]} both describe {compound} against "
                f"{receptor}. A complex is one receptor and one ligand, and two "
                f"runs of the same pair cannot share a file: rename one job, or "
                f"package them separately."
            )
        seen[cid] = source.label

        selection = str(source.config.get("traj_selection") or "")
        if selection in LIGANDLESS_SELECTIONS:
            warnings.append(
                f"{source.label}: the run wrote its trajectory with "
                f"traj_selection={selection!r}, which keeps no ligand, so there is "
                f"nothing here to measure {compound} from. Re-run with "
                f"'protein_ligand' (the default) or 'all'."
            )

        if progress is not None:
            progress(f"measuring {compound} against {receptor}")

        measured = measure.measure_complex(
            source.trajectory, complex_id=cid, receptor=receptor,
            ligand_id=ligand, ligand_resnames=resnames,
            frame_spacing_ps=source.frame_spacing_ps,
        )
        supplied = _energies_for(energies, compound, receptor)
        if supplied:
            measured.metrics.update(supplied)
            with_energies.append(cid)
        terms = _energy_terms(energies, compound, receptor)
        if terms:
            measured.panels_summary = {**(measured.panels_summary or {}),
                                       "mmpbsa_terms": terms}
        _add(package, measured, panel_lines, conformer_lines)
        if measured.depiction_painted:
            painted.append(cid)
        package.register_receptor(
            receptor,
            protein_name=str(source.structure.get("filename") or "") or None,
        )
        if measured.status != "ok":
            warnings.append(
                f"{source.label}: recorded as {measured.status}. A missing result "
                f"and a poor result are different findings, so it carries no grade.")

    if not package.complexes:
        raise PackageError(
            "nothing to package: none of the sources carried a ligand. A version 2 "
            "file is a package of receptor-ligand complexes."
        )

    # Neither layer manifest is in the spec's member layout, and DeepFlex reads
    # the receptor flexibility rows, the energy decomposition and the four strain
    # columns out of them rather than out of `metrics/`. The same numbers are
    # also in `metrics/` where the spec puts them, so a future DeepFlex that
    # reads them there finds them; FXA2_PROGRESS.md lists this as a gap for the
    # other side rather than a decision taken here.
    # Only the diagrams that actually carry mobility shading. The legend is a
    # claim about colours, so listing a plain diagram here would draw a key for
    # something that is not there.
    if painted:
        package.add_json("depictions/painted.json", painted)

    for name, lines in (("panels", panel_lines), ("conformer", conformer_lines)):
        if lines:
            package.add_text(
                f"manifests/{name}.jsonl",
                "".join(json.dumps(line, default=_plain) + "\n" for line in lines))

    if not with_energies:
        warnings.append(
            "no binding energies: FlexAppeal cannot compute MM/GBSA or an "
            "interaction energy, because both need every frame rescored in "
            "implicit solvent from a rebuilt parameterised system, and nothing in "
            "the pipeline produces one. The keys are omitted rather than written "
            "as zero, so the energy sub-score is recorded as pending and DeepFlex "
            "says so on the page. The cost is the energy panel, the MM/GBSA column "
            "of the league table and the 0.15 energy weight, which DeepFlex "
            "renormalises away. Supply them with --energies if the study has them."
        )
    elif len(with_energies) < len(package.complexes):
        warnings.append(
            f"binding energies were supplied for {len(with_energies)} of "
            f"{len(package.complexes)} complexes. The energy sub-score is a "
            f"percentile within a receptor's own cohort, so the complexes without "
            f"one take no percentile rather than a last place they did not earn."
        )

    warnings.extend(_set_kind(package))
    return package, warnings


def _plain(value: Any) -> Any:
    """numpy scalars reduced to Python ones, for the JSON-lines manifest."""
    import numpy as np  # noqa: PLC0415

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"{type(value).__name__} is not JSON-serialisable")


def _set_kind(package: fxa2.Package) -> list[str]:
    """`graded` only when every analysable complex carries every required input.

    Section 9: a `graded` file missing required inputs is refused rather than
    part-graded, because a grade with holes in it is worse than no grade, and
    `flexibility` always works and is never wrong, only less. Deciding it from
    what was actually measured rather than from an intention means a run that
    lost a measurement downgrades itself instead of claiming a grade it cannot
    support.
    """
    short: dict[str, list[str]] = {}
    for cid in package.complexes:
        metrics = json.loads(package.members[f"metrics/{cid}.json"])
        if not str(metrics.get("status") or "").startswith("ok"):
            continue
        missing = [key for key in fxa2.REQUIRED_METRIC_KEYS if metrics.get(key) is None]
        if missing:
            short[cid] = missing
    if not short:
        package.kind = "graded"
        return []

    package.kind = "flexibility"
    example = next(iter(short.values()))
    return [
        f"{len(short)} of {len(package.complexes)} complexes are missing a required "
        f"sub-score input ({', '.join(example[:4])}...), so this file declares "
        f"kind: flexibility rather than claiming a grade it cannot support. See "
        f"DEEPFLEX_FXA2_SPEC.md section 9."
    ]


def _add(package: fxa2.Package, measured: measure.Measured,
         panel_lines: list[dict[str, Any]],
         conformer_lines: list[dict[str, Any]]) -> None:
    """Write one complex's members into the package."""
    cid = measured.complex_id
    package.add_complex(cid)
    package.add_json(f"metrics/{cid}.json", measured.metrics)
    # Also as `harvest/<cid>.json`. DeepFlex's ingester refuses a package with no
    # `harvest/` and reads the binding energies, the box volume and
    # `is_reference` from there rather than from `metrics/`, while the spec's
    # member layout does not mention it at all. Written in both places, which is
    # what DeepFlex's own conformance fixture does; FXA2_PROGRESS.md lists this as
    # a spec gap for the other side rather than a decision taken here.
    package.add_json(f"harvest/{cid}.json", measured.metrics)
    package.add_json(f"pocket/{cid}.json", measured.pocket)
    for name, text in measured.structures.items():
        package.add_text(f"traj/{name}", text)
    if measured.interactions is not None:
        package.add_json(f"metrics/{cid}.interactions.json", measured.interactions)
    if measured.summary is not None:
        package.add_json(f"summary/{cid}.json", measured.summary)
    if measured.series is not None:
        package.add_arrays(f"series/{cid}.npz", **measured.series)
    if measured.panels is not None:
        package.add_arrays(f"panels/{cid}.npz", **measured.panels)
    if measured.barcode is not None:
        package.add_arrays(f"barcode/{cid}.npz", **measured.barcode)
    if measured.conformer is not None:
        package.add_arrays(f"conformer/{cid}.npz", **measured.conformer)
    if measured.panels_summary:
        panel_lines.append({"complex_id": cid, "receptor": measured.receptor,
                            **measured.panels_summary})
    if measured.conformer_summary:
        conformer_lines.append({"complex_id": cid, "receptor": measured.receptor,
                                **measured.conformer_summary})
    if measured.depiction_svg:
        package.add_text(f"depictions/{cid}.svg", measured.depiction_svg)


def build(paths: list[Path], *, out: Path, package_name: str | None = None,
          energies: Path | None = None, annotations: bool = False,
          progress=None) -> tuple[Path, list[str]]:
    """The whole command: discover, load, measure, write."""
    sources_paths = discover([Path(p) for p in paths])
    sources = [load(path) for path in sources_paths]
    supplied = _read_energies(energies) if energies else None
    package, warnings = assemble(
        sources, package_name=package_name or Path(out).stem,
        energies=supplied, progress=progress)
    if annotations:
        warnings.extend(annotate(package, progress=progress))
    return package.write(out), warnings


def _read_energies(path: Path) -> dict[str, dict[str, float]]:
    """Binding energies supplied from outside, keyed by compound name.

    FlexAppeal does not compute these: `mmpbsa_gb_dg` and the interaction-energy
    terms need an MM/GBSA or LIE rescoring of every frame, which means rebuilding
    the parameterised system in implicit solvent, and nothing in the pipeline
    produces it. Omitting the keys is the honest default and leaves `S_E`
    pending, which DeepFlex says on the page. Where a study has the numbers from
    elsewhere, this carries them in.

    The file is `{"<compound name>": {"mmpbsa_gb_dg": -40.2, ...}}`, or the same
    keyed by `<compound>/<receptor>` where one compound was run against several.
    """
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PackageError(f"{path} could not be read: {exc}") from None
    if not isinstance(document, dict):
        raise PackageError(f"{path} should be an object keyed by compound name")
    return document


# The seven keys section 4 lists under `S_E`. Only these are taken from a
# supplied file: anything else in it is ignored rather than written through, so a
# sidecar cannot be used to smuggle a field into the package.
ENERGY_KEYS = (
    "mmpbsa_gb_dg", "mmpbsa_gb_sd", "mmpbsa_pb_dg", "mmpbsa_pb_sd",
    "lie_total_mean", "lie_total_sd", "lie_drift_slope",
)


def _energies_for(energies: dict | None, compound: str,
                  receptor: str) -> dict[str, float]:
    """The supplied energies for one complex, by compound and receptor.

    Keyed on the compound name rather than the complex id, because the id is a
    salted hash that is fresh for every package file and nobody can write it down
    in advance. `<compound>/<receptor>` wins over `<compound>` where both are
    given, so one compound run against two receptors can carry two numbers.

    Report the **difference** block: the three absolute blocks are each around
    -100,000 kcal/mol and only their difference is a binding quantity. A value
    that size is refused rather than ranked.
    """
    if not energies:
        return {}
    entry = energies.get(f"{compound}/{receptor}")
    if entry is None:
        entry = energies.get(compound)
    if not isinstance(entry, dict):
        return {}

    out: dict[str, float] = {}
    for key in ENERGY_KEYS:
        value = entry.get(key)
        if value is None:
            continue
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise PackageError(f"{compound}: {key} is {value!r}, which is not a number")
        if key.endswith("_dg") and abs(float(value)) > 1000.0:
            raise PackageError(
                f"{compound}: {key} is {value}, which is an absolute MM/PBSA block "
                f"rather than the difference. Report complex minus receptor minus "
                f"ligand: only the difference is a binding quantity."
            )
        out[key] = float(value)
    return out


def _energy_terms(energies: dict | None, compound: str,
                  receptor: str) -> dict[str, dict[str, list[float]]]:
    """An optional per-term decomposition, for the energy breakdown panel.

    `{"gb": {"VDWAALS": [mean, sd, sem], "EEL": [...]}, "pb": {...}}`, the
    difference block only, exactly as `MMPBSA.py` reports it. Carried through
    rather than interpreted: these are somebody else's numbers and this writer
    has no basis for arguing with them beyond the sanity check above.
    """
    entry = energies.get(f"{compound}/{receptor}") if energies else None
    if entry is None and energies:
        entry = energies.get(compound)
    terms = (entry or {}).get("terms") if isinstance(entry, dict) else None
    if not isinstance(terms, dict):
        return {}

    out: dict[str, dict[str, list[float]]] = {}
    for method, block in terms.items():
        if not isinstance(block, dict):
            continue
        cleaned: dict[str, list[float]] = {}
        for term, values in block.items():
            if isinstance(values, (int, float)):
                values = [values]
            if not isinstance(values, list) or not values:
                continue
            if not all(isinstance(v, (int, float)) for v in values):
                raise PackageError(
                    f"{compound}: the {method} {term} term is not a list of numbers")
            cleaned[str(term)] = [float(v) for v in values[:3]]
        if cleaned:
            out[str(method).lower()] = cleaned
    return out
