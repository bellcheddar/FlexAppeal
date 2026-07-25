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
import secrets
import shutil
import time
from pathlib import Path
from typing import Any

from flask import (
    Blueprint, Flask, abort, current_app, jsonify, redirect,
    render_template, request, send_file, url_for,
)

from . import options as opts
from . import bundle, schema, sources, structure
from .options import FLEXAPPEAL_VERSION

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

    # Default to keeping every chain and every cofactor -- dropping things is a
    # deliberate act, and a user who does nothing should get the whole molecule.
    dynamic = structure.dynamic_choices(report)
    cfg["chains"] = [c["value"] for c in dynamic["chains"]]
    cfg["keep_heteroatoms"] = [
        c["value"] for c in dynamic["keep_heteroatoms"]
        if not c["help"].startswith(("crystallisation", "ion"))
    ]
    cfg = schema.normalise(cfg)

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
        built = bundle.build(
            cfg,
            structure_files[0].read_bytes(),
            saved["filename"],
            citation=saved.get("citation", ""),
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


@analysis_bp.route("/analysis", methods=["GET"])
def analysis():
    return render_template("analysis.html", active="analysis")


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

    app.register_blueprint(prepare_bp)
    app.register_blueprint(api_bp)
    app.register_blueprint(analysis_bp)

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
