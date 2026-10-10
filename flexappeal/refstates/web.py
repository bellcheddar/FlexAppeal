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

from flask import Blueprint, abort, current_app, jsonify, render_template

refstates_bp = Blueprint("refstates", __name__, url_prefix="/refstates")

# Section 6.1's state tokens. Shape as well as colour, because the palette has
# to work for a colour-blind reader and a state badge that is only a colour is
# not a label.
STATE_SHAPES = {"active": "circle", "inactive": "square",
                "intermediate": "diamond", "conflict": "triangle"}

SYSTEM_TITLES = {
    "abl1": ("ABL1 kinase in water", "Tyrosine kinase, DFG and activation loop"),
    "b2ar": ("β2AR in a POPC bilayer", "Class A GPCR, TM6 and NPxxY"),
}


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
    systems = _rows("SELECT * FROM systems ORDER BY system_id")
    cards = []
    for system in systems:
        runs = _rows("SELECT * FROM runs WHERE system_id = ? ORDER BY run_id",
                     (system["system_id"],))
        title, subtitle = SYSTEM_TITLES.get(
            system["system_id"], (system.get("name") or system["system_id"], ""))
        cards.append({
            "system": system, "runs": runs, "title": title, "subtitle": subtitle,
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
            motif["n_residues"] = len(json.loads(motif.get("residues_json") or "[]"))
        except json.JSONDecodeError:
            motif["n_residues"] = 0
    membrane = _rows(
        "SELECT * FROM membrane WHERE run_id = ? ORDER BY time_ns", (run_id,))

    from . import store

    summary = store.summarise(_database(), run_id) if _database() else {}
    return render_template(
        "refstates_run.html", run=detail, system=system, motifs=motifs,
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
    series: dict[str, dict[str, list]] = {}
    for row in rows:
        motif = series.setdefault(row["motif_id"], {
            "time_ns": [], "s": [], "label": [], "metrics": {}})
        if row["metric"] == "d_act_a":
            motif["time_ns"].append(row["time_ns"])
            motif["s"].append(row["s_score"])
            motif["label"].append(row["label"])
        motif["metrics"].setdefault(row["metric"], []).append(row["value"])
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
