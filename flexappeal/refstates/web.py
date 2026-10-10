"""The reference-states blueprint: read the explorer database, render the pages.

Spec section 6.2. Read-only at serve time, which is what section 7 asks for: the
database is written by `tools/ingest.py` during the build and this never writes
to it, so there is no migration to run on the droplet and no write lock to
contend with.

Every number on these pages comes out of the database. Section 6.2 is explicit
that the one-line finding on each card is "generated from the data, never written
by hand", and that rule is worth extending to the whole page: a hand-written
figure in a template is a figure that stops being true the next time a run is
re-analysed, and nothing would catch it.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from flask import (Blueprint, abort, current_app, jsonify, render_template,
                   send_from_directory)

from ..webapp import _RUN_ID

refstates_bp = Blueprint("refstates", __name__, url_prefix="/refstates")

# Section 6.1's state tokens. Shape as well as colour, because the palette has
# to work for a colour-blind reader and a state badge that is only a colour is
# not a label.
# What the ribbon can usefully draw. The example page learned this the hard
# way: an undecimated per-frame heatmap made one page 24 MB.
MAX_RIBBON_COLUMNS = 400

STATE_SHAPES = {"active": "circle", "inactive": "square",
                "intermediate": "diamond", "conflict": "triangle"}

SYSTEM_TITLES = {
    "abl1": ("ABL1 kinase in water", "Tyrosine kinase, DFG and activation loop"),
    "b2ar": ("β2AR in a POPC bilayer", "Class A GPCR, TM6 and NPxxY"),
}


def analysis_file(run_id: str) -> Path:
    return (Path(current_app.root_path).parent / "data" / "analysis"
            / f"{run_id}.json")


def references_file() -> Path:
    return Path(current_app.root_path).parent / "data" / "references.json"


def _database() -> Path | None:
    path = Path(current_app.config.get(
        "REFSTATES_DB",
        os.environ.get("FLEXAPPEAL_REFSTATES_DB",
                       Path(current_app.root_path).parent / "data" / "refstates.db")))
    return path if path.is_file() else None


def _rows(sql: str, params: tuple = ()) -> list[dict[str, Any]]:
    """One read-only query. An absent database gives no rows, not an error.

    The reference-states feature is additive: a deployment that has not had the
    database copied to it should show an empty section rather than a 500, which
    is also what makes the page safe to ship before the runs have finished.
    """
    from . import store

    path = _database()
    if path is None:
        return []
    with store.connect(path, read_only=True) as connection:
        return [dict(row) for row in connection.execute(sql, params)]


def finding(run_id: str) -> str:
    """The one-line finding for a card, generated from the per-frame data.

    Section 6.2's example is "TM3-TM6 distance fell from x to y". The honest
    general form is the motif that moved furthest in `s`, named with the
    direction it moved, because that is the statement the data supports whatever
    the run turned out to do.
    """
    frames = _rows(
        "SELECT motif_id, time_ns, s_score FROM motif_frames WHERE run_id = ? "
        "AND s_score IS NOT NULL ORDER BY frame", (run_id,))
    if not frames:
        return "Not yet analysed."

    by_motif: dict[str, list[tuple[float, float]]] = {}
    for row in frames:
        by_motif.setdefault(row["motif_id"], []).append(
            (row["time_ns"], row["s_score"]))

    moved = []
    for motif_id, series in by_motif.items():
        series.sort()
        if len(series) < 2:
            continue
        # The first and last tenth rather than the first and last frame: a
        # single frame at either end is one thermal fluctuation, and a finding
        # drawn from it would change every time the run was extended.
        tenth = max(1, len(series) // 10)
        start = sum(s for _, s in series[:tenth]) / tenth
        end = sum(s for _, s in series[-tenth:]) / tenth
        moved.append((abs(end - start), motif_id, start, end, series[-1][0]))
    if not moved:
        return "Not yet analysed."

    delta, motif_id, start, end, last_ns = max(moved)
    name = motif_id.split(":", 1)[-1]
    if delta < 0.05:
        return (f"Over {last_ns:.0f} ns no motif moved by more than 0.05 on the "
                f"state score; {name} changed least of all, holding at "
                f"{end:+.2f}.")
    direction = "towards active" if end > start else "towards inactive"
    return (f"Over {last_ns:.0f} ns the {name} state score moved from "
            f"{start:+.2f} to {end:+.2f}, {direction}.")


@refstates_bp.get("/")
def index():
    """Section 6.2's Examples tab: one card per system."""
    from . import store

    systems = _rows("SELECT * FROM systems ORDER BY system_id")
    cards = []
    for system in systems:
        runs = _rows("SELECT * FROM runs WHERE system_id = ? ORDER BY run_id",
                     (system["system_id"],))
        title, subtitle = SYSTEM_TITLES.get(
            system["system_id"], (system.get("name") or system["system_id"], ""))
        # The pair, where both halves exist. The point of each system is two
        # runs from opposite start states, so the result that matters is the
        # separation between them rather than either one on its own: two runs
        # that stay apart have held their states, and two that meet have not.
        scores = {run["run_id"]: store.whole_score(_database(), run["run_id"])
                  for run in runs} if _database() else {}
        by_state = {}
        for run in runs:
            score = scores.get(run["run_id"]) or {}
            if score.get("final_s") is not None and run.get("start_state"):
                by_state[run["start_state"]] = (run["run_id"], score)
        pair = None
        active = next((v for k, v in by_state.items()
                       if k in ("active", "DFG-in")), None)
        inactive = next((v for k, v in by_state.items()
                         if k in ("inactive", "DFG-out")), None)
        if active and inactive:
            gap = active[1]["final_s"] - inactive[1]["final_s"]
            pair = {
                "active_run": active[0], "inactive_run": inactive[0],
                "active_s": active[1]["final_s"],
                "inactive_s": inactive[1]["final_s"],
                "separation": gap,
                # Both on the correct side of the intermediate band is the
                # result the build is looking for; anything else is reported
                # as what it is rather than rounded towards it.
                "held": (active[1]["final_s"] > 0.2
                         and inactive[1]["final_s"] < -0.2),
            }

        cards.append({
            "system": system, "runs": runs, "title": title, "subtitle": subtitle,
            "scores": scores, "pair": pair,
            "findings": {run["run_id"]: finding(run["run_id"]) for run in runs},
            "n_references": (_rows(
                "SELECT COUNT(*) AS n FROM reference_structures WHERE "
                "system_id = ?", (system["system_id"],)) or [{"n": 0}])[0]["n"],
        })
    return render_template("refstates_index.html", cards=cards, active="refstates",
                           have_database=_database() is not None)


@refstates_bp.get("/run/<run_id>")
def run(run_id: str):
    """Section 6.2's trajectory explorer, for one run."""
    rows = _rows("SELECT * FROM runs WHERE run_id = ?", (run_id,))
    if not rows:
        abort(404)
    detail = rows[0]
    system = (_rows("SELECT * FROM systems WHERE system_id = ?",
                    (detail["system_id"],)) or [{}])[0]
    motifs = _rows("SELECT * FROM motifs WHERE system_id = ? ORDER BY name",
                   (detail["system_id"],))
    for motif in motifs:
        # Counted here rather than in the template: Jinja has no JSON filter in
        # this app, and adding one for a length is more surface than it is worth.
        try:
            residues = json.loads(motif.get("residues_json") or "[]")
        except json.JSONDecodeError:
            residues = []
        motif["n_residues"] = len(residues)
        # The span rather than the whole list: a 35-residue helix listed in
        # full is a wall of numbers, and the endpoints are what locate it.
        motif["span"] = (f"{min(residues)} to {max(residues)}"
                         if residues else "")
    membrane = _rows(
        "SELECT * FROM membrane WHERE run_id = ? ORDER BY time_ns", (run_id,))

    from . import store

    summary = store.summarise(_database(), run_id) if _database() else {}

    # Is there a reference map to draw? Asked here rather than letting the
    # panel render and the fetch 404: the browser logs a failed request as a
    # console error, so an empty panel costs a real error on a real page.
    has_map = any(
        json.loads(row.get("metrics_json") or "{}").get(
            "dfg_d_alphac" if detail["system_id"] == "abl1" else "tm3_tm6")
        is not None
        for row in _rows("SELECT metrics_json FROM reference_structures "
                         "WHERE system_id = ?", (detail["system_id"],)))

    viewer: dict[str, Any] = {}
    analysis = analysis_file(run_id)
    if analysis.is_file():
        try:
            viewer = json.loads(analysis.read_text()).get("viewer") or {}
        except json.JSONDecodeError:
            viewer = {}
    # Only offer the panel when both files are actually there: a viewer that
    # fails to load is worse than one that is not shown.
    directory = Path(current_app.root_path).parent / "data" / "analysis"
    if not all((directory / f"{run_id}.{kind}").is_file()
               for kind in ("pdb", "xtc")):
        viewer = {}
    return render_template(
        "refstates_run.html", run=detail, system=system, motifs=motifs,
        viewer=viewer, has_map=has_map,
        membrane=membrane, summary=summary, finding=finding(run_id),
        shapes=STATE_SHAPES, active="refstates",
        core=json.loads(system.get("notes") or "{}").get("core", {}))


@refstates_bp.get("/run/<run_id>/frames.json")
def frames(run_id: str):
    """Per-frame series for the timelines, one request for the whole run.

    One request rather than one per motif: section 6.3 wants first paint under
    two seconds, and a 50 ns run is a few thousand rows, which is far cheaper to
    send once than to fetch six times.
    """
    rows = _rows(
        "SELECT frame, time_ns, motif_id, metric, value, s_score, label "
        "FROM motif_frames WHERE run_id = ? ORDER BY frame", (run_id,))
    if not rows:
        abort(404)
    # Keyed on the frame, not on one metric's name. Keying on `d_act_a` meant a
    # motif that resolved only its inactive distance came back with eight
    # metric rows and no points at all: a 200, a trace with empty x and y, and
    # an invisible line on the timeline. The rule that an unmeasured key is
    # omitted rather than zeroed makes that a legal row, so the reader has to
    # cope with it.
    by_motif: dict[str, dict[int, dict]] = {}
    for row in rows:
        frame = by_motif.setdefault(row["motif_id"], {}).setdefault(
            row["frame"], {"time_ns": row["time_ns"], "s": row["s_score"],
                           "label": row["label"], "metrics": {}})
        frame["metrics"][row["metric"]] = row["value"]

    series: dict[str, dict[str, list]] = {}
    for motif_id, frames_by_index in by_motif.items():
        order = sorted(frames_by_index)
        metrics_seen = sorted({m for f in frames_by_index.values()
                               for m in f["metrics"]})
        series[motif_id] = {
            "time_ns": [frames_by_index[i]["time_ns"] for i in order],
            "s": [frames_by_index[i]["s"] for i in order],
            "label": [frames_by_index[i]["label"] for i in order],
            # None where a frame did not resolve that metric, so the arrays
            # stay the same length as the time axis and a gap stays a gap.
            "metrics": {name: [frames_by_index[i]["metrics"].get(name)
                               for i in order]
                        for name in metrics_seen},
        }
    return jsonify({"run_id": run_id, "motifs": series})


@refstates_bp.get("/references/<system_id>")
def references(system_id: str):
    """Section 6.2's reference browser, where the positive controls are visible."""
    rows = _rows("SELECT * FROM reference_structures WHERE system_id = ? "
                 "ORDER BY our_state, pdb", (system_id,))
    if not rows:
        abort(404)
    for row in rows:
        try:
            row["metrics"] = json.loads(row.get("metrics_json") or "{}")
        except json.JSONDecodeError:
            row["metrics"] = {}
    # Comparable means both sides state an opinion. Counting a structure this
    # build could not label as a disagreement would read as 17 of 20 where the
    # control's real result is 17 of 17 with three unlabelled, which are
    # different claims.
    comparable = [r for r in rows if r.get("source_state") and r.get("our_state")]
    agree = [r for r in comparable
             if r["source_state"].lower() == r["our_state"].lower()]
    unlabelled = [r for r in rows if r.get("source_state") and not r.get("our_state")]
    return render_template(
        "refstates_references.html", references=rows, system_id=system_id,
        n_agree=len(agree), n_comparable=len(comparable),
        n_unlabelled=len(unlabelled), active="refstates")

@refstates_bp.get("/map/<system_id>.json")
def state_map(system_id: str):
    """Section 5.6's 2D state map: the reference cloud, by the catalogue's label.

    For the kinase this is the Dunbrack plane, the DFG-Phe ring measured against
    its two anchors, which is the plot the published criterion is drawn on and
    the one where D010's threshold failure is visible: the two clusters sit 8.5 A
    apart on the difference and the published 11.0 A line cuts through the
    corner of one of them.

    For the receptor the x axis is the TM3-TM6 cleft distance, which is the axis
    the Dror 2011 state maps use.

    The reference coordinates come from the phase 1 manifest, where they were
    measured once, rather than being recomputed per request: they are properties
    of the structures and do not change between page loads.
    """
    rows = _rows("SELECT pdb, source_state, our_state, metrics_json FROM "
                 "reference_structures WHERE system_id = ? ORDER BY pdb",
                 (system_id,))
    if not rows:
        abort(404)

    points = []
    for row in rows:
        try:
            metrics = json.loads(row.get("metrics_json") or "{}")
        except json.JSONDecodeError:
            continue
        if system_id == "abl1":
            x, y = metrics.get("dfg_d_alphac"), metrics.get("dfg_d_lys")
            axes = ("DFG-Phe to alphaC-Glu+4 Calpha (A)",
                    "DFG-Phe to catalytic Lys Calpha (A)")
        else:
            x, y = metrics.get("tm3_tm6"), metrics.get("npxxy_to_inactive")
            axes = ("TM3-TM6 Calpha distance 3.50 to 6.34 (A)",
                    "NPxxY RMSD to the nearest inactive reference (A)")
        if x is None or (system_id != "abl1" and y is None):
            continue
        points.append({"pdb": row["pdb"], "x": x, "y": y,
                       "state": row.get("source_state") or row.get("our_state")})

    if not points:
        abort(404)
    return jsonify({"system_id": system_id, "axes": axes, "points": points,
                    # D010's line, drawn so a reader can see where it falls
                    # rather than being told about it.
                    "published_threshold": 11.0 if system_id == "abl1" else None})


@refstates_bp.get("/about")
def about():
    """Section 6.2's About page: methods, software and verified references.

    The reference list is read from `data/references.json`, which
    `tools/verify_references.py` wrote by resolving every citation through
    Crossref during the build. Serving this page touches no network, and a
    citation that could not be verified is shown as unverified rather than
    omitted, because dropping it silently would leave the method uncited and
    look deliberate.
    """
    rows: list[dict[str, Any]] = []
    path = references_file()
    if path.is_file():
        try:
            rows = json.loads(path.read_text())
        except json.JSONDecodeError:
            rows = []
    verified = [r for r in rows if r.get("resolved")]

    vendored = ""
    sha = Path(current_app.root_path).parent / "vendor" / "BOLTZMAKER_SHA.txt"
    if sha.is_file():
        vendored = sha.read_text().strip().split()[0][:12]

    return render_template("refstates_about.html", references=rows,
                           n_verified=len(verified), boltzmaker_sha=vendored,
                           active="refstates")


@refstates_bp.get("/run/<run_id>/dssp.json")
def dssp(run_id: str):
    """Per-frame secondary structure, for the ribbon in section 6.2.

    Served from the run's analysis file rather than from SQLite. Section 7 keeps
    per-residue secondary structure out of the database because it is a
    residues-by-frames matrix rather than a set of rows, and at one row per
    residue per frame a 50 ns run would be millions of them.

    Decimated to at most 400 columns, which is what the heatmap can draw. The
    full resolution stays in the analysis file for anything that needs it.
    """
    if not _RUN_ID.match(run_id or ""):
        abort(404)
    path = analysis_file(run_id)
    if not path.is_file():
        abort(404)
    try:
        payload = json.loads(path.read_text()).get("dssp") or {}
    except json.JSONDecodeError:
        abort(404)
    if not payload.get("frames"):
        abort(404)

    frames = payload["frames"]
    stride = max(1, len(frames) // MAX_RIBBON_COLUMNS)
    return jsonify({
        "run_id": run_id,
        "residues": payload.get("residues", []),
        "frames": frames[::stride],
        "stride": stride,
        "n_frames_total": payload.get("n_frames"),
        "motifs": payload.get("motifs", {}),
    })


# The only two things the viewer may fetch out of the analysis directory.
VIEWER_SUFFIXES = {"pdb": "chemical/x-pdb", "xtc": "application/octet-stream"}


@refstates_bp.get("/run/<run_id>/viewer.<kind>")
def viewer_file(run_id: str, kind: str):
    """The topology or the strided coordinates for the 3D viewer.

    Served by suffix from a fixed allowlist rather than by filename, so the
    route cannot be talked into returning the run's JSON, the database, or
    anything else that happens to sit in that directory. The run id is checked
    against the same pattern the rest of the app uses.
    """
    if not _RUN_ID.match(run_id or "") or kind not in VIEWER_SUFFIXES:
        abort(404)
    directory = (Path(current_app.root_path).parent / "data" / "analysis")
    name = f"{run_id}.{kind}"
    if not (directory / name).is_file():
        abort(404)
    return send_from_directory(directory, name,
                               mimetype=VIEWER_SUFFIXES[kind])
