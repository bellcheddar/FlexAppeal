"""The FlexAppeal web application.

Flask app factory plus blueprints, following the house pattern established by
AlphaFraud's ``webapp.py`` and reused by BoltzMaker-web: server-rendered Jinja,
form POSTs rather than a JSON API for anything that carries a file, vanilla JS
only, and a per-request scratch directory cleaned in a ``finally``.

The one place this departs from BoltzMaker-web is session state. Preparing a run
is genuinely multi-step -- load a structure, then configure against what was
found in it -- so a structure has to outlive the request that uploaded it. It is
kept in a scratch directory keyed by an opaque token that the form carries in a
hidden field. No cookies, no server-side session, no database: the token is the
only handle, and a systemd timer sweeps anything left behind.
"""

from __future__ import annotations

import io
import json
import os
import re
import secrets
import shutil
import string
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from flask import (
    Blueprint, Flask, abort, current_app, jsonify, redirect,
    render_template, request, send_file, send_from_directory, url_for,
)

from . import options as opts
from . import analysis, bundle, fxa, measure, plots, schema, sources, structure
from .options import FLEXAPPEAL_VERSION
from .refstates import panels as refstates_panels

# MUST stay in sync with `client_max_body_size` in deploy/nginx-flexappeal.conf.
# The two limits protect different things -- nginx rejects the body before it
# reaches Python, Flask catches anything that slips past -- and a mismatch means
# one of them produces an unhelpful error instead of the friendly 413 page.
MAX_CONTENT_LENGTH = 250 * 1024 * 1024

# Uploaded structures are capped far lower than the results files the Analysis
# tab accepts. A coordinate file bigger than this is not a single protein.
MAX_STRUCTURE_BYTES = 25 * 1024 * 1024

# How long a prepared structure stays in scratch before the sweeper may remove
# it. Long enough to configure a run without being rushed; short enough that an
# abandoned upload does not sit on disk overnight.
SESSION_TTL_SECONDS = 4 * 60 * 60

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent
DEFAULT_SCRATCH_ROOT = REPO_ROOT / "web_scratch"


# ===========================================================================
#  Scratch sessions
# ===========================================================================


def _scratch_root() -> Path:
    return Path(current_app.config["SCRATCH_ROOT"])


def new_session() -> tuple[str, Path]:
    """Create a scratch directory and the opaque token that addresses it."""
    token = secrets.token_urlsafe(16)
    path = _scratch_root() / token
    path.mkdir(parents=True, exist_ok=False)
    return token, path


def session_path(token: str) -> Path:
    """Resolve a token to its directory, refusing anything that escapes scratch.

    The token comes back from the browser in a hidden field, so it is untrusted
    input that names a path -- exactly the shape of a directory-traversal bug.
    Resolving and then checking containment is the check that actually holds,
    rather than filtering for '..' and hoping.
    """
    if not token or len(token) > 64:
        abort(400, "invalid session token")
    root = _scratch_root().resolve()
    path = (root / token).resolve()
    if not path.is_relative_to(root) or not path.is_dir():
        abort(400, "this preparation session has expired -- please load your structure again")
    return path


def sweep_scratch(root: Path, ttl: int = SESSION_TTL_SECONDS) -> int:
    """Remove expired sessions. Called by the systemd timer and on startup."""
    removed = 0
    cutoff = time.time() - ttl
    if not root.is_dir():
        return 0
    for entry in root.iterdir():
        try:
            if entry.is_dir() and entry.stat().st_mtime < cutoff:
                shutil.rmtree(entry, ignore_errors=True)
                removed += 1
        except OSError:
            continue
    return removed


def _save_session(path: Path, report: structure.StructureReport,
                  data: bytes, filename: str, citation: str = "") -> None:
    (path / "input").with_suffix(Path(filename).suffix or ".pdb").write_bytes(data)
    (path / "report.json").write_text(json.dumps({
        "report": report.to_dict(),
        "filename": filename,
        "citation": citation,
    }, indent=2))


def _load_session(path: Path) -> dict[str, Any]:
    try:
        return json.loads((path / "report.json").read_text())
    except (OSError, json.JSONDecodeError):
        abort(400, "this preparation session is incomplete -- please load your structure again")


# ===========================================================================
#  Prepare
# ===========================================================================

prepare_bp = Blueprint("prepare", __name__)


def _render_prepare(cfg: dict[str, Any], *, report: dict | None = None,
                    token: str = "", issues: list | None = None,
                    dynamic: dict | None = None, citation: str = "",
                    source_warning: str = "") -> str:
    estimate = None
    if report:
        estimate = report.get("_estimate")
    return render_template(
        "prepare.html",
        active="prepare",
        groups=opts.GROUPS,
        by_group=opts.BY_GROUP,
        cfg=cfg,
        report=report,
        token=token,
        issues=issues or [],
        dynamic=dynamic or {},
        citation=citation,
        source_warning=source_warning,
        estimate=estimate,
        derived=schema.derive(cfg),
        wall=schema.estimate_wall_time(cfg),
        is_active=schema.is_active,
        active_choices=schema.active_choices,
    )


@prepare_bp.route("/", methods=["GET"])
def index():
    return _render_prepare(opts.defaults())


@prepare_bp.route("/prepare/structure", methods=["POST"])
def load_structure():
    """Step one: acquire a structure and report what is in it."""
    raw = request.form.to_dict(flat=True)
    source = raw.get("input_source", "upload")
    cfg = opts.defaults() | {k: v for k, v in raw.items() if k in opts.BY_ID}

    try:
        if source == "upload":
            uploaded = request.files.get("input_file")
            if not uploaded or not uploaded.filename:
                return _render_prepare(cfg, issues=[schema.Issue(
                    "error", "choose a PDB or mmCIF file to upload.", "input_file")])
            data = uploaded.read(MAX_STRUCTURE_BYTES + 1)
            if len(data) > MAX_STRUCTURE_BYTES:
                return _render_prepare(cfg, issues=[schema.Issue(
                    "error",
                    f"that file is larger than {MAX_STRUCTURE_BYTES // (1024 * 1024)} MB. "
                    f"FlexAppeal is built for single proteins and small complexes.",
                    "input_file")])
            filename = uploaded.filename
            citation = ""
            source_warning = ""
        else:
            validation = schema.validate(cfg, form_post=True)
            accession_errors = [
                i for i in validation.errors if i.option_id in ("pdb_id", "uniprot_id")
            ]
            if accession_errors:
                return _render_prepare(cfg, issues=accession_errors)
            result = sources.fetch(cfg)
            data, filename, citation = result.data, result.filename, result.citation
            source_warning = sources.SOURCE_WARNINGS.get(result.source, "")

        report = structure.analyse(
            data, filename, disulfide_cutoff=float(cfg.get("disulfide_cutoff") or 2.5)
        )
    except (sources.SourceError, structure.StructureError) as exc:
        return _render_prepare(cfg, issues=[schema.Issue("error", str(exc), "input_file")])

    token, path = new_session()
    _save_session(path, report, data, filename, citation)

    # Default to keeping every chain and every real ligand -- dropping things is
    # a deliberate act, and a user who does nothing should get the whole
    # molecule. `default_heteroatoms` owns which of them count: it used to be a
    # prefix match on the help text here, which could not see that a haem has an
    # iron in it and so defaulted metalloproteins into a configuration that
    # cannot be parameterised.
    dynamic = structure.dynamic_choices(report)
    cfg["chains"] = [c["value"] for c in dynamic["chains"]]
    cfg["keep_heteroatoms"], hetero_notes = structure.default_heteroatoms(report)
    cfg = schema.normalise(cfg)
    # A kept ligand brings its own two measurements with it.
    #
    # Both are offered by the registry and gated on `has_ligands`, and neither
    # was ticked, so a user who uploaded a structure with a ligand, changed
    # nothing, and ran it got every binding-site panel and no ligand RMSD or
    # contact occupancy. Keeping the ligand is already a default rather than a
    # question, which is the same reasoning: a measurement of the thing the run
    # is about should not need finding.
    #
    # Seeded here rather than in the registry's own default, because the
    # choices are only valid for a run that has a ligand and an apo run would
    # fail validation on them.
    if cfg.get("has_ligands"):
        chosen = list(cfg.get("analysis_metrics") or [])
        for metric in ("ligand_rmsd", "ligand_contacts"):
            if metric not in chosen:
                chosen.append(metric)
        cfg["analysis_metrics"] = chosen
    if hetero_notes:
        source_warning = " ".join(
            [source_warning] + hetero_notes).strip() if source_warning \
            else " ".join(hetero_notes)

    report_dict = report.to_dict()
    report_dict["_estimate"] = structure.estimate_system_size(report, cfg)
    cfg["_estimated_atoms"] = report_dict["_estimate"]["total_atoms"]
    cfg["_solute_atoms"] = report.solute_atoms
    if not cfg.get("job_name") or cfg["job_name"] == "flexappeal_run":
        cfg["job_name"] = Path(filename).stem.replace(".", "_")[:64] or "flexappeal_run"

    return _render_prepare(cfg, report=report_dict, token=token, dynamic=dynamic,
                           citation=citation, source_warning=source_warning)


def _dynamic_from_report(report_dict: dict[str, Any]) -> dict[str, list[dict[str, str]]]:
    """Rebuild the dynamic choice lists from a stored report.

    The full StructureReport is not rehydrated because only these two lists are
    needed to re-render the form, and reconstructing the dataclass from JSON
    would couple the session format to the dataclass's field order.
    """
    return {
        "chains": [
            {"value": c["id"], "label": f"Chain {c['id']}",
             "help": f"{c['kind']}, {c['observed_residues']} residues"}
            for c in report_dict["chains"] if c["kind"] != "water"
        ],
        "keep_heteroatoms": [
            {"value": h["name"],
             "label": f"{h['name']}{' × ' + str(h['count']) if h['count'] > 1 else ''}",
             "help": h["description"] or h["category"]}
            for h in report_dict["heteroatoms"] if h["category"] != "water"
        ],
    }


@prepare_bp.route("/prepare/build", methods=["POST"])
def build():
    """Step two: validate the configuration and emit the bundle."""
    token = request.form.get("_token", "")
    path = session_path(token)
    saved = _load_session(path)

    # A multiselect posts one key repeatedly; everything else posts it once.
    # to_dict() would silently keep only the first value of a multiselect, so
    # the widget type decides whether we take the list or the scalar.
    raw: dict[str, Any] = {}
    for key in request.form:
        values = request.form.getlist(key)
        opt = opts.BY_ID.get(key)
        raw[key] = values if (opt and opt.widget == "multiselect") else values[0]

    result = schema.validate(raw, form_post=True)
    report_dict = saved["report"]
    cfg = result.config
    cfg["_solute_atoms"] = report_dict["solute_atoms"]
    cfg["_estimated_atoms"] = (report_dict.get("_estimate") or {}).get("total_atoms", 0)

    if not result.ok:
        return _render_prepare(cfg, report=report_dict, token=token,
                               issues=result.issues,
                               dynamic=_dynamic_from_report(report_dict))

    structure_files = [p for p in path.glob("input.*")]
    if not structure_files:
        abort(400, "the uploaded structure is no longer in this session -- please load it again")

    try:
        # The chemical definitions, fetched here rather than on the user's
        # machine: a PDB HETATM record has no bond orders, and the run script
        # transfers them from these onto the deposited coordinates. Without them
        # it warns that it has no definition and then fails at system creation
        # with a missing template, several minutes into somebody's run.
        ligands, ligand_warnings = bundle.collect_ligands(cfg, report_dict)
        built = bundle.build(
            cfg,
            structure_files[0].read_bytes(),
            saved["filename"],
            citation=saved.get("citation", ""),
            ligands=ligands,
        )
    except bundle.BundleError as exc:
        return _render_prepare(
            cfg, report=report_dict, token=token,
            issues=[schema.Issue("error", str(exc))],
            dynamic=_dynamic_from_report(report_dict),
        )

    return send_file(
        io.BytesIO(built.content),
        # Not application/x-sh: some browsers and mail gateways will try to be
        # helpful about shell scripts. An opaque stream is downloaded verbatim,
        # which is what a base64 payload needs.
        mimetype="application/octet-stream",
        as_attachment=True,
        download_name=built.filename,
    )


# ===========================================================================
#  API -- live readouts
# ===========================================================================

api_bp = Blueprint("api", __name__, url_prefix="/api")


@api_bp.route("/estimate", methods=["POST"])
def estimate():
    """Recompute the derived readouts as the user changes the form.

    The only JSON endpoint in the app, and it exists because these three numbers
    -- atom count, trajectory size, wall time -- are what stop someone
    accidentally configuring a three-week run. Computing them server-side means
    the arithmetic lives in schema.derive() and nowhere else.
    """
    raw = request.get_json(silent=True)
    if not isinstance(raw, dict):
        # json.loads("3") returns an int, and a bare list is valid JSON too --
        # guard the type before calling .get on it.
        return jsonify({"error": "expected a JSON object"}), 400

    cfg = schema.validate(raw).config
    cfg["_estimated_atoms"] = int(raw.get("_estimated_atoms") or 0)
    cfg["_solute_atoms"] = int(raw.get("_solute_atoms") or 0)

    result = schema.validate(cfg)
    derived = schema.derive(cfg)
    wall = schema.estimate_wall_time(cfg)

    return jsonify({
        "derived": derived,
        "wall": wall,
        "errors": [{"option": i.option_id, "message": i.message} for i in result.errors],
        "warnings": [{"option": i.option_id, "message": i.message} for i in result.warnings],
        "active": [o.id for o in schema.active_options(cfg)],
    })


# ===========================================================================
#  Analysis (Phase 4)
# ===========================================================================

analysis_bp = Blueprint("analysis", __name__)

MAX_FXA_BYTES = MAX_CONTENT_LENGTH


@analysis_bp.route("/analysis", methods=["GET", "POST"])
def analysis_page():
    # Named `analysis_page`, not `analysis`: a view function called `analysis`
    # sits at module scope and shadows the imported `analysis` module, so
    # `analysis.METRICS` resolves against the function and raises
    # AttributeError. The collision only shows up when a module attribute is
    # actually touched, which is why it survived until the re-analysis route.
    if request.method == "GET":
        return render_template("analysis.html", active="analysis", results=None)

    uploaded = request.files.get("fxa_file")
    if not uploaded or not uploaded.filename:
        return render_template("analysis.html", active="analysis", results=None,
                               error="Choose a .fxa results file to upload.")

    content = uploaded.read(MAX_FXA_BYTES + 1)
    if len(content) > MAX_FXA_BYTES:
        return render_template(
            "analysis.html", active="analysis", results=None,
            error=f"That file is larger than {MAX_FXA_BYTES // (1024 * 1024)} MB. "
                  f"Rebuild the bundle with a smaller payload tier.")

    try:
        results = fxa.load(content)
    except fxa.FxaError as exc:
        return render_template("analysis.html", active="analysis", results=None,
                               error=str(exc))

    # The token lets the browser fetch the structure and trajectory for the
    # viewer without re-uploading them, and gives /reanalyse something to work
    # from in phase 7.
    token, path = new_session()
    (path / "results.fxa").write_bytes(content)

    # Kept in the library only when the uploader asks. A results file can carry
    # a receptor and a compound somebody has not published, and this site is
    # public, so nothing an upload brings is listed for other visitors without
    # that being an explicit choice.
    kept = None
    if request.form.get("keep_run"):
        kept = keep_run(content, label=request.form.get("run_label", ""))

    return render_template(
        "analysis.html",
        active="analysis",
        kept_run=kept,
        **_results_context(results, content, token))


# The binding panels are measured from the trajectory rather than read out of
# the file, which costs seconds, and the example page re-renders on every
# request. Keyed on a digest of the file so the example is measured once per
# process and a re-uploaded file is not measured twice; the built panels are
# cached rather than the measurement, because they are the part the template
# needs and they are two orders of magnitude smaller.
_BINDING_CACHE: dict[str, dict[str, Any]] = {}
_BINDING_CACHE_MAX = 4


# The same record on disk, so the measurement is paid once rather than once per
# worker per restart.
#
# Measured on the live site: a run with a ligand takes 10 to 13 seconds to
# measure its binding site, and 0.5 to render once the answer is in hand. The
# memory cache above is per worker, so with two workers the first two visitors
# to each bound run paid it after every deploy, and the committed example was
# as slow as the two reference-state rows. A page that takes twelve seconds is
# a page a reader assumes is broken.
#
# Keyed on the SHA-256 of the results file, which is what the memory cache is
# keyed on: the measurement is a pure function of the file, so a file that has
# not changed cannot have a different answer, and a file that has is a
# different key and a different entry. Nothing is invalidated, because nothing
# can go stale.
def _binding_cache_dir() -> Path:
    return Path(current_app.config["RUNS_ROOT"]) / ".binding-cache"


def _binding_from_disk(key: str) -> dict[str, Any] | None:
    """The stored record, or None. Never raises: this is a cache."""
    try:
        return json.loads((_binding_cache_dir() / f"{key}.json").read_text())
    except (OSError, ValueError, KeyError, RuntimeError):
        return None


def _binding_to_disk(key: str, record: dict[str, Any]) -> None:
    """Store the record beside the library, atomically and quietly.

    Written to a temporary name and moved into place, so a worker reading the
    file while another writes it sees either the old one or the new one and
    never half of either. A failure here costs nothing but the next request's
    time, so it is swallowed: a read-only disk must not take down a page that
    would otherwise render.
    """
    try:
        directory = _binding_cache_dir()
        directory.mkdir(parents=True, exist_ok=True)
        temporary = directory / f".{key}.{os.getpid()}.tmp"
        temporary.write_text(json.dumps(record))
        temporary.replace(directory / f"{key}.json")
    except (OSError, ValueError, TypeError, RuntimeError):
        current_app.logger.warning("could not cache the binding panels")

# The mobility legend's swatches, read from the measurement module rather than
# written out here. The diagram is painted with those exact values, so a second
# list would be one edit away from a legend that describes colours the diagram
# does not have.
MOBILITY_BANDS = measure.MOBILITY_BANDS
MOBILITY_CEILING = measure.MOBILITY_CEILING


def _binding_panels(results, content: bytes, label: str) -> dict[str, Any]:
    """The binding-site panels for one results file, with the ligand diagram.

    Returns an empty record for a run with no ligand: that is not a failure, it
    is a run with no binding site. A measurement that raises costs the panels
    and not the page -- a traceback here would take down a page that works.
    """
    import base64  # noqa: PLC0415
    import hashlib  # noqa: PLC0415

    key = hashlib.sha256(content).hexdigest()
    if key in _BINDING_CACHE:
        return _BINDING_CACHE[key]
    stored = _binding_from_disk(key)
    if stored is not None:
        _BINDING_CACHE[key] = stored
        return stored

    record: dict[str, Any] = {"panels": [], "notes": [], "depiction": None,
                              "painted": False, "grade": None}
    try:
        from flexappeal import grade as grade_mod  # noqa: PLC0415
        from flexappeal import package as package_mod  # noqa: PLC0415

        measured, reasons = package_mod.measure_results(results, label)
        if measured is None:
            record["notes"] = [r for r in reasons if "no heteroatoms" not in r]
        else:
            record["panels"] = plots.build_binding(measured)
            record["notes"] = list(measured.notes)
            record["painted"] = bool(measured.depiction_painted)
            # Graded here and not in the package: spec section 7 reserves the
            # grade to DeepFlex on ingest so there is only ever one copy of the
            # judgement in a file. This one is rendered and discarded.
            record["grade"] = grade_mod.meter(grade_mod.score_one(measured))
            if measured.depiction_svg:
                # A data URI rather than a route: the diagram is a few kilobytes
                # and belongs to the page that measured it, so serving it
                # separately would mean a second endpoint and a second lookup
                # for something already in hand.
                encoded = base64.b64encode(
                    measured.depiction_svg.encode("utf-8")).decode("ascii")
                record["depiction"] = f"data:image/svg+xml;base64,{encoded}"
    except Exception as exc:                                  # noqa: BLE001
        current_app.logger.exception("binding panels failed")
        record["notes"] = [f"the binding-site panels could not be built ({exc})"]

    if len(_BINDING_CACHE) >= _BINDING_CACHE_MAX:
        _BINDING_CACHE.pop(next(iter(_BINDING_CACHE)))
    _BINDING_CACHE[key] = record
    _binding_to_disk(key, record)
    return record


# ---------------------------------------------------------------------------
#  The worked example
# ---------------------------------------------------------------------------
#
# A real 100 ns run of T4 lysozyme L99A with benzene, committed to the repository
# and rendered through the same code path as an upload. Deliberately not a
# screenshot or a written-up transcript: the panels below are built by
# plots.build_all() from the .fxa the run actually produced, so a change that
# breaks the Analysis tab breaks this page too and cannot ship looking fine.
#
# A complex rather than an apo protein because an apo run draws none of the
# binding-site panels, so an apo example would leave most of the tab untested.

# Which committed run the Example tab shows. One name, because the directory,
# the results file and the bundle are all named after it: changing the example
# used to mean finding three strings that had to agree.
EXAMPLE_NAME = "lysozyme_benzene_100ns"
EXAMPLE_DIR = REPO_ROOT / "examples" / EXAMPLE_NAME
EXAMPLE_FXA = EXAMPLE_DIR / "output" / f"{EXAMPLE_NAME}.fxa"
EXAMPLE_BUNDLE = f"examples/{EXAMPLE_NAME}/flexappeal_{EXAMPLE_NAME}.command"

# The Mol* viewer and /reanalyse address the results through a session token,
# so the example needs one. Reused across requests rather than minted per hit,
# which would fill the scratch directory with identical copies; recreated if the
# sweeper has since expired it.
_EXAMPLE_TOKEN: list[str] = []


def _example_session() -> str:
    if _EXAMPLE_TOKEN:
        # Deliberately NOT session_path(): that helper aborts with 400 when the
        # directory is missing, which is right for a token off a form and wrong
        # here, where a missing directory is the expected state after the
        # sweeper has run and the answer is to make a new one. Calling it turned
        # "the example session expired" into a 400 page four hours after the
        # first visit, rather than the transparent re-copy intended.
        candidate = _scratch_root() / _EXAMPLE_TOKEN[0]
        if (candidate / "results.fxa").is_file():
            return _EXAMPLE_TOKEN[0]
        _EXAMPLE_TOKEN.clear()
    token, path = new_session()
    shutil.copyfile(EXAMPLE_FXA, path / "results.fxa")
    _EXAMPLE_TOKEN.append(token)
    return token


def _example_option_tables(config: dict) -> list[dict]:
    """The run's settings, grouped exactly as the Prepare form groups them.

    Read out of the registry rather than hand-written, so a new option appears
    here the day it is added and cannot be forgotten. Options the run left at
    their default are marked rather than hidden -- for a reference example the
    interesting fact is often that a value was NOT changed.

    An option the file does not carry is listed too, and says so. The committed
    example was run before the reference-state settings existed, and skipping
    what it has no value for quietly undid the promise in the paragraph above:
    five options were added and the page that is supposed to list every one of
    them simply did not. What it must not do is print the registry's default in
    that row, because the run did not choose it; the value is left empty and
    the row is marked as not recorded.
    """
    defaults = opts.defaults()
    tables = []
    for group in opts.GROUPS:
        rows = []
        for opt in opts.OPTIONS:
            if opt.group != group.id:
                continue
            recorded = opt.id in config
            value = config.get(opt.id)
            if not recorded:
                shown = "—"
            elif isinstance(value, bool):
                shown = "yes" if value else "no"
            elif isinstance(value, list):
                shown = ", ".join(str(v) for v in value) or "—"
            else:
                shown = str(value) if str(value) != "" else "—"
            rows.append({
                "label": opt.label,
                "value": shown,
                "units": opt.units or "" if recorded else "",
                "openmm": (opt.openmm or "").split(".")[-1],
                "changed": recorded and value != defaults.get(opt.id),
                "recorded": recorded,
                "help": opt.help,
            })
        if rows:
            tables.append({"title": group.title, "icon": group.icon,
                           "blurb": group.blurb, "rows": rows,
                           "changed": sum(1 for r in rows if r["changed"])})
    return tables


def _example_captures() -> dict[str, list[str]]:
    """The recorded terminal output, grouped by which script produced it.

    Filenames carry their phase (run-01-startup, analyse-02-packing), so the
    grouping is read off the name rather than inferred from the numbering --
    which would silently mis-sort the moment a stage is added to either half.
    """
    from markupsafe import Markup

    directory = EXAMPLE_DIR / "screenshots"
    grouped: dict[str, list[str]] = {"run": [], "analyse": []}
    if not directory.is_dir():
        return grouped
    for path in sorted(directory.glob("*.html")):
        phase = path.name.split("-", 1)[0]
        if phase in grouped:
            grouped[phase].append(Markup(path.read_text()))
    return grouped


def _results_context(results, content: bytes, token: str) -> dict[str, Any]:
    """Everything a results file renders into, built once for every page.

    The Analysis tab, the Example tab and a run in the library all show the same
    measurements of the same kind of file, and they drifted the moment there
    were two copies of this dict: the binding panels reached the Analysis tab a
    release before the Example tab, because only one of the two routes had been
    given them. One builder, three callers.
    """
    summary = fxa.summarise(results)
    # Parsed once and used twice: the convergence panels shade the head and the
    # tiles print where it ends, and the two reading the file separately is how
    # a page ends up shading one window and naming another.
    state = (plots.parse_state_data(results.members["state_data.csv"])
             if "state_data.csv" in results.members else {})
    settled = plots.equilibration(state)
    return {
        "results": results,
        "summary": summary,
        "token": token,
        "tiles": plots.tiles(summary, results.metrics, settled),
        "binding": _binding_panels(results, content, results.job_name),
        "mobility_bands": MOBILITY_BANDS,
        "mobility_ceiling": MOBILITY_CEILING,
        "panels": plots.build_all(results)["panels"],
        "convergence": plots.convergence(state),
        "settled": settled,
        "has_viewer": bool(results.topology_pdb),
        "warnings": results.warnings,
        "reanalysis_metrics": analysis.METRICS,
        # The reference-state panels, for a results file that carries the
        # analysis they are built from. Empty for every other run, so a page
        # that has nothing extra to show simply does not show a section. In
        # here rather than in the two routes for the reason the docstring
        # above gives: three callers, one builder.
        "refstates": refstates_panels.from_results(results),
        # What the viewer's own extra controls are built from: the motif
        # colouring, one ghost per reference state, and the bilayer. Empty for
        # every run that carries no reference-state analysis, which is how the
        # viewer card stays one piece of markup for all of them.
        "refstates_viewer": refstates_panels.viewer_controls(results),
    }


@analysis_bp.route("/runs", methods=["GET"])
def runs_page():
    """The library: every kept run, newest first, the example always present."""
    return render_template("runs.html", active="runs", runs=list_runs())


@analysis_bp.route("/runs/<run_id>", methods=["GET"])
def run_page(run_id: str):
    """One run, rendered with the same panels an upload gets.

    Its own session token, so the viewer can fetch the structure and the
    trajectory the same way it does for an upload. A run in the library is
    durable; the scratch copy the viewer reads from is not, and is remade on
    demand when the sweeper has been through.
    """
    content = run_content(run_id)
    try:
        results = fxa.load(content, verify_checksums=False)
    except fxa.FxaError as exc:
        abort(404, str(exc))

    token, path = new_session()
    (path / "results.fxa").write_bytes(content)

    label = run_id
    if run_id != EXAMPLE_RUN_ID:
        try:
            label = json.loads(
                (run_directory(run_id) / "run.json").read_text()
            ).get("label") or results.job_name
        except (OSError, json.JSONDecodeError):
            label = results.job_name
    else:
        label = results.job_name

    return render_template(
        "run_detail.html", active="runs", run_id=run_id, run_label=label,
        **_results_context(results, content, token))


@analysis_bp.route("/example", methods=["GET"])
def example_page():
    if not EXAMPLE_FXA.is_file():
        return render_template(
            "example_missing.html", active="example",
            path=EXAMPLE_FXA.relative_to(REPO_ROOT),
            rebuild=f"examples/{EXAMPLE_NAME}/rebuild.sh"), 200

    content = EXAMPLE_FXA.read_bytes()
    results = fxa.load(content, verify_checksums=False)
    config = results.manifest.get("config", {})

    return render_template(
        "analysis.html",
        active="example",
        example={
            "tables": _example_option_tables(config),
            "captures": _example_captures(),
            "manifest": results.manifest,
            "bundle": EXAMPLE_BUNDLE,
        },
        **_results_context(results, content, _example_session()))


# ---------------------------------------------------------------------------
#  The runs library
# ---------------------------------------------------------------------------
#
# A directory of kept results files, one per run, listed on the Runs tab so a
# reader can open any of them and get the same panels an upload gets. Separate
# from `web_scratch`, which is the opposite thing: scratch is per-visit, opaque
# and swept, and this is durable and addressed by a name a person can read.
#
# Deliberately a directory rather than a database. A run is one `.fxa` plus a
# line of metadata, the listing is a handful of entries, and a schema migration
# for that would cost more than it saves.
DEFAULT_RUNS_ROOT = REPO_ROOT / "web_runs"

# The committed example is always in the library and does not live in that
# directory: it is served straight out of `examples/`, so the two cannot drift
# and a fresh checkout has a populated Runs tab with no setup at all.
EXAMPLE_RUN_ID = "example"

_RUN_ID = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}$")
# The characters `_RUN_ID` will accept in first position. A generated id has to
# start with one of these or it is refused by the two functions that resolve
# it, so the first character is chosen from here rather than hoped for.
_RUN_ID_FIRST = string.ascii_letters + string.digits


def _runs_root() -> Path:
    return Path(current_app.config["RUNS_ROOT"])


def run_directory(run_id: str) -> Path:
    """Resolve a run id to its directory, refusing anything that escapes.

    The id arrives in a URL, so it is untrusted: the pattern admits no dot and
    no separator, and the resolved path is checked against the root as well,
    which is belt and braces on purpose. A traversal here would serve any file
    on the box through the viewer's structure route.
    """
    if not _RUN_ID.match(run_id or ""):
        abort(404)
    root = _runs_root().resolve()
    path = (root / run_id).resolve()
    if path != root and root not in path.parents:
        abort(404)
    return path


def _run_record(run_id: str, content: bytes, meta: dict) -> dict:
    """One row of the listing, read from the results file itself.

    From the file rather than from the metadata beside it, so a row cannot
    describe something the file does not contain.
    """
    try:
        results = fxa.load(content, verify_checksums=False)
    except fxa.FxaError as exc:
        return {"id": run_id, "label": meta.get("label") or run_id,
                "broken": str(exc)}
    summary = fxa.summarise(results)
    return {
        "id": run_id,
        "label": meta.get("label") or results.job_name,
        "job_name": results.job_name,
        "duration_ns": summary.get("duration_ns"),
        "frames": summary.get("frames_packed"),
        "atoms": summary.get("atoms_packed"),
        "created_at": meta.get("added_at") or summary.get("created_at") or "",
        "tier": summary.get("tier"),
        "has_ligands": bool(results.has_ligands),
        "source": meta.get("source") or "upload",
        "broken": None,
    }


def list_runs() -> list[dict]:
    """Every run in the library, the committed example first.

    A directory that cannot be read is skipped rather than fatal: the listing is
    a convenience and one bad entry must not cost the page.
    """
    rows: list[dict] = []
    if EXAMPLE_FXA.is_file():
        rows.append(_run_record(
            EXAMPLE_RUN_ID, EXAMPLE_FXA.read_bytes(),
            {"label": EXAMPLE_NAME, "source": "committed example"}))

    root = _runs_root()
    if root.is_dir():
        for path in sorted(root.iterdir(), reverse=True):
            if not path.is_dir() or not _RUN_ID.match(path.name):
                continue
            stored = path / "results.fxa"
            if not stored.is_file():
                continue
            meta = {}
            try:
                meta = json.loads((path / "run.json").read_text())
            except (OSError, json.JSONDecodeError):
                pass
            try:
                rows.append(_run_record(path.name, stored.read_bytes(), meta))
            except OSError:
                continue
    return rows


def run_content(run_id: str) -> bytes:
    """The results file for one run, or a 404."""
    if run_id == EXAMPLE_RUN_ID:
        if not EXAMPLE_FXA.is_file():
            abort(404)
        return EXAMPLE_FXA.read_bytes()
    stored = run_directory(run_id) / "results.fxa"
    if not stored.is_file():
        abort(404)
    return stored.read_bytes()


def keep_run(content: bytes, label: str = "") -> str:
    """Add an uploaded results file to the library, and return its id.

    The id is a random token rather than the job name. Two runs can share a
    name, and a name chosen by an uploader would otherwise decide where a file
    lands on disk.

    The first character is drawn separately, which is not fussiness.
    `token_urlsafe` draws from a 64-character alphabet that includes `-` and
    `_`, and this used to be

        secrets.token_urlsafe(9).replace("_", "-").replace("=", "")[:16]

    so roughly one id in thirty-two began with a hyphen. `_RUN_ID` requires the
    first character to be alphanumeric, and it guards both `list_runs` and
    `run_directory`. A run with such an id was therefore written to disk, the
    page said "Kept in the library", and it appeared in no listing and resolved
    to no URL: kept, invisible and unreachable, about three per cent of the
    time. Measured at 3.03% over 20,000 draws against the 3.12% the alphabet
    predicts.

    It showed up as a test that failed about once in thirty full runs, which
    reads as flakiness and was a product fault the whole time.
    """
    run_id = (secrets.choice(_RUN_ID_FIRST)
              + secrets.token_urlsafe(9).replace("_", "-").replace("=", "")[:15])
    path = _runs_root() / run_id
    path.mkdir(parents=True, exist_ok=True)
    (path / "results.fxa").write_bytes(content)
    (path / "run.json").write_text(json.dumps({
        "label": (label or "").strip()[:120],
        "added_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "upload",
    }, indent=2))
    return run_id


def _load_results(token: str) -> fxa.Results:
    path = session_path(token)
    stored = path / "results.fxa"
    if not stored.is_file():
        abort(400, "this analysis session has expired -- please upload your results file again")
    return fxa.load(stored.read_bytes(), verify_checksums=False)


@analysis_bp.route("/analysis/<token>/structure", methods=["GET"])
def structure_file(token: str):
    """The topology, for the Mol* viewer."""
    results = _load_results(token)
    if not results.topology_pdb:
        abort(404)
    return send_file(io.BytesIO(results.topology_pdb), mimetype="chemical/x-pdb",
                     download_name="topology.pdb")


@analysis_bp.route("/analysis/<token>/trajectory", methods=["GET"])
def trajectory_file(token: str):
    """The decimated trajectory, for playback in the viewer."""
    results = _load_results(token)
    if not results.trajectory_xtc:
        abort(404)
    return send_file(io.BytesIO(results.trajectory_xtc),
                     mimetype="application/octet-stream",
                     download_name="traj.xtc")


# The only other members the viewer may fetch, and the type each is served as.
# A fixed allowlist rather than a pattern, which is the discipline the page this
# replaced used and for the same reason: the route takes a name from a URL, and
# a pattern is how the protection gets lost later. Nothing here can name
# `refstates.json`, the manifest, or the metrics.
VIEWER_MEMBERS = {
    "refstates_active.pdb": "chemical/x-pdb",
    "refstates_inactive.pdb": "chemical/x-pdb",
    "lipids.pdb": "chemical/x-pdb",
    "lipids.xtc": "application/octet-stream",
}


@analysis_bp.route("/analysis/<token>/member/<member>", methods=["GET"])
def viewer_member(token: str, member: str):
    """A reference overlay or the bilayer, for the viewer's own toggles.

    These ride inside the results file rather than beside it, so a run in the
    library is one file and a deploy is that file. The overlay was checked
    against the analysis when it was built, to 0.005 angstrom, and withheld
    where the two disagreed, so a ghost that is here is one the measurements
    endorse.
    """
    if member not in VIEWER_MEMBERS:
        abort(404)
    results = _load_results(token)
    payload = (results.members or {}).get(member)
    if not payload:
        abort(404)
    return send_file(io.BytesIO(payload), mimetype=VIEWER_MEMBERS[member],
                     download_name=member)


# ===========================================================================
#  Factory
# ===========================================================================


def create_app() -> Flask:
    app = Flask(
        __name__,
        template_folder=str(PACKAGE_ROOT / "templates"),
        static_folder=str(PACKAGE_ROOT / "static"),
    )
    app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH

    scratch_root = Path(os.environ.get("FLEXAPPEAL_SCRATCH_ROOT", str(DEFAULT_SCRATCH_ROOT)))
    scratch_root.mkdir(parents=True, exist_ok=True)
    app.config["SCRATCH_ROOT"] = scratch_root
    runs_root = Path(os.environ.get("FLEXAPPEAL_RUNS_ROOT", str(DEFAULT_RUNS_ROOT)))
    runs_root.mkdir(parents=True, exist_ok=True)
    app.config["RUNS_ROOT"] = runs_root
    sweep_scratch(scratch_root)

    app.jinja_env.globals["FLEXAPPEAL_VERSION"] = FLEXAPPEAL_VERSION
    app.jinja_env.globals["OPENMM_VERSION"] = opts.OPENMM_VERSION

    @app.context_processor
    def _asset_helper():
        """Cache-bust static files by mtime so nginx can serve them immutably.

        Lifted from AlphaFraud's webapp.py -- the vendored Mol* and Plotly
        bundles are several megabytes each, and far-future caching is only safe
        if the URL changes when the file does.
        """
        def asset(filename: str) -> str:
            try:
                version = int(os.path.getmtime(os.path.join(app.static_folder, filename)))
            except OSError:
                version = 0
            return url_for("static", filename=filename, v=version)
        return {"asset": asset}

    @app.get("/favicon.ico")
    def favicon():
        """Served at the root as well as from /static.

        Browsers request /favicon.ico whatever the link tags say, and a 404 per
        page view is noise in the access log that looks like a fault.
        """
        return send_from_directory(
            Path(app.root_path) / "static" / "icons", "favicon.ico",
            mimetype="image/vnd.microsoft.icon")

    app.register_blueprint(prepare_bp)
    app.register_blueprint(api_bp)
    app.register_blueprint(analysis_bp)
    # There was a fifth blueprint here, `refstates`, serving a tab of its own:
    # an Examples card, a trajectory explorer, a reference browser and a method
    # page, all reading one SQLite file. The four runs it existed for are rows
    # in the Runs library now, with their reference-state measurements as panels
    # at the end of their own pages, so the tab, the blueprint, its four
    # templates and its four JSON endpoints are gone and the analysis library
    # they all read stays. D035.

    @app.after_request
    def _no_html_caching(response):
        """Make HTML always revalidate.

        This is the other half of asset() and without it that function does
        nothing. Static files are served `immutable, max-age=31536000`, which
        is only safe because the ?v=mtime in their URL changes when the file
        does -- but the URL lives in the HTML. Flask sends no Cache-Control on
        a rendered template, so browsers fall back to heuristic freshness
        (commonly a tenth of the age since Last-Modified) and keep serving the
        cached page, which keeps asking for the *old* asset URL, which really
        is cached for a year. A CSS fix then stays invisible for days and looks
        like a deploy that did not take -- which is exactly how it presented.

        no-cache is not no-store: the browser may still hold the page, it just
        has to revalidate it, so an unchanged page is a 304 rather than a
        refetch.
        """
        if response.mimetype == "text/html":
            response.headers.setdefault("Cache-Control", "no-cache")
        return response

    @app.route("/healthz")
    def healthz():
        return {"status": "ok", "version": FLEXAPPEAL_VERSION}, 200

    @app.errorhandler(400)
    def bad_request(exc):
        return render_template(
            "error.html", code=400, title="That request did not make sense",
            message=getattr(exc, "description", str(exc)), active="",
        ), 400

    @app.errorhandler(404)
    def not_found(exc):
        return render_template(
            "error.html", code=404, title="Nothing here",
            message="That page does not exist.", active="",
        ), 404

    @app.errorhandler(413)
    def too_large(exc):
        return render_template(
            "error.html", code=413, title="That file is too large",
            message=f"Uploads are limited to "
                    f"{MAX_CONTENT_LENGTH // (1024 * 1024)} MB.",
            active="",
        ), 413

    @app.errorhandler(500)
    def server_error(exc):
        app.logger.exception("unhandled error")
        return render_template(
            "error.html", code=500, title="Something broke",
            message="That is a bug in FlexAppeal rather than anything you did.",
            active="",
        ), 500

    return app


# ---------------------------------------------------------------------------
#  Re-analysis (the hybrid path)
# ---------------------------------------------------------------------------


@analysis_bp.route("/analysis/<token>/reanalyse", methods=["POST"])
def reanalyse(token: str):
    """Start a bounded re-analysis in a detached subprocess.

    Never does the work inline. MDTraj holds the GIL inside C extensions, so a
    contact map computed in the request thread would block a gunicorn worker for
    the duration; and a detached process survives a service restart. This is
    AlphaFraud's `POST /calculate/run` pattern, with a lock file instead of a
    database row because there is no database here.
    """
    path = session_path(token)
    stored = path / "results.fxa"
    if not stored.is_file():
        abort(400, "this analysis session has expired -- please upload your results file again")

    raw = request.get_json(silent=True)
    if not isinstance(raw, dict):
        raw = request.form.to_dict(flat=False)
        raw = {k: (v if k == "metrics" else v[0]) for k, v in raw.items()}

    try:
        parsed = analysis.parse_request(raw)
    except analysis.ReanalysisError as exc:
        return jsonify({"status": "error", "message": str(exc)}), 400

    try:
        analysis.acquire_lock(_scratch_root())
    except analysis.Busy as exc:
        return jsonify({"status": "busy", "message": str(exc)}), 429

    request_path = path / "reanalyse_request.json"
    output_path = path / "reanalyse_result.json"
    request_path.write_text(json.dumps(parsed.to_dict()))
    output_path.unlink(missing_ok=True)

    script = REPO_ROOT / "FlexAppeal.py"
    try:
        subprocess.Popen(
            [sys.executable, str(script), "reanalyse",
             str(stored), str(request_path), str(output_path)],
            cwd=str(REPO_ROOT),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            # A new session so the worker outlives a gunicorn restart rather
            # than being killed with its parent's process group.
            start_new_session=True,
        )
    except OSError as exc:
        analysis.release_lock(_scratch_root())
        return jsonify({"status": "error",
                        "message": f"could not start the analysis: {exc}"}), 500

    return jsonify({"status": "running"}), 202


@analysis_bp.route("/analysis/<token>/reanalyse/status", methods=["GET"])
def reanalyse_status(token: str):
    """Poll for the detached worker's result."""
    path = session_path(token)
    output_path = path / "reanalyse_result.json"

    if not output_path.is_file():
        request_path = path / "reanalyse_request.json"
        if not request_path.is_file():
            return jsonify({"status": "idle"})
        age = time.time() - request_path.stat().st_mtime
        if age > analysis.TIMEOUT_SECONDS:
            analysis.release_lock(_scratch_root())
            return jsonify({
                "status": "error",
                "message": f"the analysis did not finish within "
                           f"{analysis.TIMEOUT_SECONDS} seconds and was abandoned.",
            })
        return jsonify({"status": "running", "seconds": round(age)})

    try:
        payload = json.loads(output_path.read_text())
    except (OSError, json.JSONDecodeError):
        # The worker may be mid-write; treat it as still running rather than
        # reporting a failure that will resolve itself on the next poll.
        return jsonify({"status": "running"})

    if payload.get("status") == "ready":
        payload["figures"] = plots.build_reanalysis(payload.get("metrics", {}))
    return jsonify(payload)
