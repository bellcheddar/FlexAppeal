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
import sqlite3
import os
from pathlib import Path
from typing import Any

from flask import (Blueprint, abort, current_app, jsonify, render_template,
                   request, send_from_directory, url_for)

from . import membrane as membrane_mod
from . import kinase as kinase_mod
from . import pocket, sources
from .state import DECISIVE_S
from ..webapp import _RUN_ID

refstates_bp = Blueprint("refstates", __name__, url_prefix="/refstates")

# Section 6.1's state tokens. Shape as well as colour, because the palette has
# to work for a colour-blind reader and a state badge that is only a colour is
# not a label.
# What the ribbon can usefully draw. The example page learned this the hard
# way: an undecimated per-frame heatmap made one page 24 MB.
MAX_RIBBON_COLUMNS = 400

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
    figures = Path(current_app.root_path).parent / "data" / "figures"
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
                "held": (active[1]["final_s"] > DECISIVE_S
                         and inactive[1]["final_s"] < -DECISIVE_S),
            }

        cards.append({
            "system": system, "runs": runs, "title": title, "subtitle": subtitle,
            "scores": scores, "pair": pair,
            "findings": {run["run_id"]: finding(run["run_id"]) for run in runs},
            "n_references": (_rows(
                "SELECT COUNT(*) AS n FROM reference_structures WHERE "
                "system_id = ?", (system["system_id"],)) or [{"n": 0}])[0]["n"],
            # Offered only when the file is actually there, the same rule the
            # viewer and the bilayer toggle follow: an <img> whose source 404s
            # costs a console error and a broken icon on a real page.
            "thumbnail": (
                url_for("refstates.thumbnail", system_id=system["system_id"])
                if (figures / f"{system['system_id']}_thumb.png").is_file()
                else None),
        })
    return render_template("refstates_index.html", cards=cards, active="refstates",
                           decisive_s=DECISIVE_S,
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
        # The endpoints as numbers too, for the ribbon's motif bands. The
        # string above is for a reader; a chart needs the ends to draw
        # between, and parsing them back out of the sentence would be absurd.
        motif["first"] = min(residues) if residues else None
        motif["last"] = max(residues) if residues else None
    membrane = _rows(
        "SELECT * FROM membrane WHERE run_id = ? ORDER BY time_ns", (run_id,))
    # Section 5.8 asks for a panel a user can read the bilayer's health off,
    # and five numbers against two windows held in the reader's head is not
    # that. The verdict is not a stored column, so it is recomputed from the
    # numbers by the module that owns the windows.
    from .membrane import row_verdicts

    for chunk in membrane:
        chunk["verdicts"] = row_verdicts(chunk.get("apl"),
                                         chunk.get("thickness"))

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

    # Section 6.2's motif table asks for "nearest reference PDB at end". The
    # store does not have it: section 7 gives `motif_frames` no column for it,
    # and the nearest reference is a property of a frame rather than of the
    # run. It is in the analysis payload, which this view already opens for the
    # viewer bundle, so it is read from there and merged into the summary rows
    # rather than the schema growing a column the spec does not have.
    #
    # Both are given, not one. The score is a comparison of two distances, and
    # a table that named only the nearer of them would hide exactly the case
    # the comparison exists for: a motif sitting between the two states.
    nearest: dict[str, dict[str, Any]] = {}
    analysis_path = analysis_file(run_id)
    if analysis_path.is_file():
        try:
            frames = json.loads(analysis_path.read_text()).get("frames") or []
        except json.JSONDecodeError:
            frames = []
        for motif in (frames[-1].get("motifs") if frames else []) or []:
            name = motif.get("motif")
            if not name:
                continue
            found = {key: motif.get(key) for key in
                     ("nearest_active", "nearest_inactive",
                      "d_act_a", "d_inact_a")
                     if motif.get(key) is not None}
            if found:
                nearest[f"{detail['system_id']}:{name}"] = found
    for key, found in nearest.items():
        if summary.get("motifs", {}).get(key) is not None:
            summary["motifs"][key].update(found)

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
    # Same rule for the bilayer: offered only when both of its files are
    # actually there, so the toggle cannot produce a failed request. A water
    # run has none, and the template simply does not draw the control.
    elif not all((directory / f"{run_id}.lipids.{kind}").is_file()
                 for kind in ("pdb", "xtc")):
        viewer.pop("lipids", None)
    return render_template(
        "refstates_run.html", run=detail, system=system, motifs=motifs,
        viewer=viewer, has_map=has_map,
        measurements=frame_measurements(run_id),
        sse_rows=sse_agreement(run_id),
        membrane=membrane, summary=summary, finding=finding(run_id),
        active="refstates",
        # The band boundary the score plot shades, passed rather than written
        # into the template. It is the same number the label rule uses, and a
        # plot whose shading disagreed with the labels beside it would be worse
        # than one with no shading at all.
        decisive_s=DECISIVE_S,
        dfg_published_a=kinase_mod.DFG_PUBLISHED_A,
        # The bilayer windows the badges are decided against, passed rather
        # than written into the prose. The note said "60 to 68" and the badge
        # read `APL_EXPECTED_A2`, so the two could part without either looking
        # wrong, and a reader checking a badge against the note would have been
        # checking it against a number nothing uses.
        apl_expected=membrane_mod.APL_EXPECTED_A2,
        apl_acceptable=membrane_mod.APL_ACCEPTABLE_A2,
        thickness_expected=membrane_mod.THICKNESS_EXPECTED_A,
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
    # Section 5.6's whole-protein score, residue-weighted, with the smoothed
    # version beside it rather than instead of it: the spec asks for smoothing
    # "for display only" and for the raw values to be what is stored.
    whole: dict[str, Any] = {}
    path = analysis_file(run_id)
    if path.is_file():
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError:
            payload = {}
        raw = payload.get("whole_s_raw")
        if raw:
            times = [f.get("time_ns") for f in payload.get("frames") or []]
            whole = {"time_ns": times[:len(raw)], "raw": raw,
                     "smoothed": payload.get("whole_s_smoothed"),
                     "window_ns": payload.get("smoothing_window_ns"),
                     "skipped": payload.get("smoothing_skipped")}

    return jsonify({"run_id": run_id, "motifs": series, "whole": whole})


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
        # The cuts the two catalogues were filtered at, from the filters that
        # applied them. Written out as prose the page could describe a set it
        # had not been given.
        kinase_max_resolution_a=sources.KINASE_MAX_RESOLUTION_A,
        kinase_min_quality=sources.KINASE_MIN_QUALITY,
        kinase_max_missing=sources.KINASE_MAX_MISSING_POCKET,
        gpcr_max_resolution_a=sources.GPCR_MAX_RESOLUTION_A,
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

    # Section 6.2 asks for the trajectory's path across the same plane, and for
    # clicking a point to jump the viewer to that frame, so each point carries
    # its frame index as well as its coordinates.
    path = []
    run_id = request.args.get("run")
    if run_id and _RUN_ID.match(run_id):
        analysis = analysis_file(run_id)
        if analysis.is_file():
            try:
                frames = json.loads(analysis.read_text()).get("frames") or []
            except json.JSONDecodeError:
                frames = []
            for index, frame in enumerate(frames):
                found = frame.get("measurements") or {}
                if system_id == "abl1":
                    plane = found.get("dfg_plane") or {}
                    x, y = plane.get("d_alphac_a"), plane.get("d_lys_a")
                else:
                    x = (found.get("micro_switches") or {}).get("tm3_tm6_a")
                    # The NPxxY distance to the nearest inactive reference is
                    # already computed as that motif's score, so it is read
                    # from there rather than measured twice.
                    y = next((m.get("d_inact_a") for m in frame.get("motifs", [])
                              if m.get("motif") == "NPxxY"), None)
                if x is None or y is None:
                    continue
                path.append({"frame": index, "time_ns": frame.get("time_ns"),
                             "x": x, "y": y})

    return jsonify({"system_id": system_id, "axes": axes, "points": points,
                    "path": path,
                    # D010's line, drawn so a reader can see where it falls
                    # rather than being told about it.
                    "published_threshold": (kinase_mod.DFG_PUBLISHED_A
                                            if system_id == "abl1" else None)})


def ensemble_counts() -> dict[str, dict[str, int]]:
    """How many references each ensemble actually holds, per system.

    The About page stated these in prose and one pair had gone stale: it said
    the receptor was served by 72 active and 28 inactive structures where the
    ensembles hold 84 and 32. Nothing could have caught that, which is the
    whole trouble with a count written into a sentence. A page comparing the
    two evidence bases should count them.

    The rule here is the one `sources` builds the ensembles with, not the
    catalogue's own column: ABL1's active ensemble is DFG-in **and** alphaC-in,
    which is 5 of the 6 entries KLIFS calls DFG-in, while its inactive ensemble
    is DFG-out. The receptor's is GPCRdb's state as given.

    Returns an empty dict when there is no database, which the page renders by
    leaving the sentence out rather than printing zeroes.
    """
    from . import store

    database = _database()
    if database is None:
        return {}
    counts: dict[str, dict[str, int]] = {}
    try:
        with store.connect(database, read_only=True) as connection:
            for row in connection.execute(
                    "SELECT system_id, source_state, metrics_json "
                    "FROM reference_structures"):
                try:
                    metrics = json.loads(row["metrics_json"] or "{}")
                except json.JSONDecodeError:
                    metrics = {}
                if row["system_id"] == "abl1":
                    if (metrics.get("their_dfg") == "in"
                            and metrics.get("their_alphac") == "in"):
                        state = "active"
                    elif metrics.get("their_dfg") == "out":
                        state = "inactive"
                    else:
                        continue
                else:
                    state = (row["source_state"] or "").lower()
                    if state not in ("active", "inactive"):
                        continue
                counts.setdefault(row["system_id"], {})
                counts[row["system_id"]][state] = (
                    counts[row["system_id"]].get(state, 0) + 1)
    except sqlite3.Error:
        return {}
    return counts


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
                           r_spine_published_a=pocket.R_SPINE_CONTACT_CRYSTAL_A,
                           r_spine_a=pocket.R_SPINE_CONTACT_A,
                           dfg_published_a=kinase_mod.DFG_PUBLISHED_A,
                           counts=ensemble_counts(),
                           active="refstates")


# How the per-frame measurements are presented: a label, where to find the
# value, and its unit. A table rather than code per measurement, so adding one
# to the analysis does not mean editing the page.
MEASUREMENT_ROWS = [
    ("alphaC salt bridge", ("alphac_salt_bridge", "distance_a"), "\u00c5",
     ("alphac_salt_bridge", "state")),
    ("alphaC axis rotation", ("alphac_displacement", "axis_rotation_deg"), "\u00b0", None),
    ("alphaC centroid shift", ("alphac_displacement", "centroid_shift_a"), "\u00c5", None),
    ("P-loop tip to hinge", ("p_loop", "tip_to_hinge_a"), "\u00c5", None),
    ("R-spine", ("r_spine", "state"), "", None),
    ("TM3\u2013TM6 cleft", ("micro_switches", "tm3_tm6_a"), "\u00c5", None),
    ("Ionic lock", ("micro_switches", "ionic_lock", "distance_a"), "\u00c5",
     ("micro_switches", "ionic_lock", "state")),
    ("Y5.58\u2013Y7.53", ("micro_switches", "y5x58_y7x53_oh_a"), "\u00c5", None),
    ("TM6 axis rotation", ("tm6_rotation", "axis_rotation_deg"), "\u00b0", None),
    ("TM6 centroid shift", ("tm6_rotation", "centroid_shift_a"), "\u00c5", None),
    ("PIF connector RMSD", ("connector", "rmsd_a"), "\u00c5", None),
    ("Helix 8 to the membrane plane", ("helix8", "angle_to_membrane_plane_deg"),
     "\u00b0", None),
]


def _dig(payload: dict, path: tuple):
    """Follow a key path, returning None rather than raising on a gap."""
    found: Any = payload
    for key in path:
        if not isinstance(found, dict):
            return None
        found = found.get(key)
    return found


def sse_agreement(run_id: str) -> list[dict[str, Any]]:
    """Per motif, how far its secondary structure matches each reference set.

    Section 5.6's comparison. Reported with the margin, because secondary
    structure is far less discriminating than geometry and a lean on a small
    margin is noise: most residues of a helix are helical in both states.
    """
    path = analysis_file(run_id)
    if not path.is_file():
        return []
    try:
        motifs = (json.loads(path.read_text()).get("dssp") or {}).get("motifs") or {}
    except json.JSONDecodeError:
        return []

    rows = []
    for name, found in sorted(motifs.items()):
        agreement = found.get("reference_agreement") or {}
        if not agreement.get("n_residues"):
            continue
        rows.append({
            "motif": name,
            "n_residues": agreement["n_residues"],
            "active": agreement["mean_active_agreement"],
            "inactive": agreement["mean_inactive_agreement"],
            "margin": agreement.get("margin"),
            "leans": agreement.get("leans"),
            "decisive": agreement.get("decisive"),
            "helix_fraction": found.get("helix_fraction_last"),
        })
    return rows


def frame_measurements(run_id: str) -> list[dict[str, Any]]:
    """The per-frame readings at the first and last frame, for the page.

    First and last rather than every frame: the table says what moved over the
    run, and the per-frame series belongs in a plot. A reading absent from both
    ends is left out rather than shown blank, so the table is what was measured
    rather than a list of what might have been.
    """
    path = analysis_file(run_id)
    if not path.is_file():
        return []
    try:
        frames = json.loads(path.read_text()).get("frames") or []
    except json.JSONDecodeError:
        return []
    usable = [f for f in frames if f.get("measurements")]
    if not usable:
        return []

    first, last = usable[0]["measurements"], usable[-1]["measurements"]
    rows = []
    for label, value_path, unit, state_path in MEASUREMENT_ROWS:
        start, end = _dig(first, value_path), _dig(last, value_path)
        if start is None and end is None:
            continue
        row = {"label": label, "unit": unit, "start": start, "end": end,
               "numeric": isinstance(end if end is not None else start,
                                     (int, float))}
        if state_path:
            row["state"] = _dig(last, state_path)
        if row["numeric"] and start is not None and end is not None:
            row["change"] = end - start
        rows.append(row)
    return rows


def sse_consensus(payload: dict, residues: list[int]) -> dict[str, list]:
    """The reference ensembles' secondary structure, on the ribbon's own axis.

    Section 6.2 asks the DSSP ribbon for "active/inactive reference consensus
    rows above it", so that a reader can see whether the trajectory's secondary
    structure looks like the active references or the inactive ones rather than
    only what it is.

    `report.sse_profiles` already holds the data, as a count of DSSP codes per
    state, per motif, per position within that motif. Two things have to happen
    to put it on the ribbon: the modal code has to be picked, and the motif
    position has to be resolved to a residue number through `motif_residues`,
    which is the same index space.

    The agreement fraction comes back with the code. A residue where 70 of 72
    references agree and one where 38 do are not the same claim, and the
    difference is the reader's to judge.
    """
    profiles = (payload.get("report") or {}).get("sse_profiles") or {}
    motif_residues = payload.get("motif_residues") or {}
    out: dict[str, list] = {}
    for state, by_motif in profiles.items():
        found: dict[int, dict] = {}
        for motif, positions in (by_motif or {}).items():
            numbers = motif_residues.get(motif) or []
            for index, counts in (positions or {}).items():
                try:
                    residue = numbers[int(index)]
                except (ValueError, IndexError):
                    continue
                if not counts:
                    continue
                code, best = max(counts.items(), key=lambda kv: kv[1])
                total = sum(counts.values())
                found[residue] = {"code": code.strip() or "C",
                                  "agreement": round(best / total, 3)
                                  if total else None}
        # Over the ribbon's residues, so the two rows line up with the heatmap
        # rather than needing the browser to join them.
        out[state] = [found.get(residue) for residue in residues]
    return out


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
        whole = json.loads(path.read_text())
    except json.JSONDecodeError:
        abort(404)
    payload = whole.get("dssp") or {}
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
        "consensus": sse_consensus(whole, payload.get("residues", [])),
    })


# The only things the viewer may fetch out of the analysis directory.
# The lipid pair is listed explicitly rather than matched by pattern. The whole
# point of this allowlist is that the route cannot be talked into serving
# something else out of that directory, and a pattern is how that protection
# gets lost later.
VIEWER_SUFFIXES = {"pdb": "chemical/x-pdb", "xtc": "application/octet-stream",
                   "lipids.pdb": "chemical/x-pdb",
                   "lipids.xtc": "application/octet-stream"}


# Section 6.2 asks the Examples cards for a thumbnail. They live in
# `data/figures` with the rest of the generated figures rather than in the
# static tree, because they are output and not an asset, so they need a route.
# Same discipline as `viewer_file`: a fixed set of names, not a pattern.
THUMBNAILS = {"abl1", "b2ar"}


@refstates_bp.get("/thumbnail/<system_id>.png")
def thumbnail(system_id: str):
    """One system's card thumbnail, drawn by `tools/make_figures.py`."""
    if system_id not in THUMBNAILS:
        abort(404)
    directory = Path(current_app.root_path).parent / "data" / "figures"
    name = f"{system_id}_thumb.png"
    if not (directory / name).is_file():
        abort(404)
    return send_from_directory(directory, name, mimetype="image/png")


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
