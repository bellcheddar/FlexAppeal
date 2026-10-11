#!/usr/bin/env python3
"""Static figures for FINDINGS.md, from the explorer database.

Section 11 asks for FINDINGS "with figures". The explorer draws its charts in
the browser, which is right for a page a user scrubs through and no use at all
in a document, so these are generated separately and committed: a figure that
has to be regenerated to be read is a figure nobody reads.

Three kinds, one per thing the build actually claims:

  <run>_motifs.png      the per-motif state score against time, which is the
                        headline result for a row
  <system>_map.png      the reference structures on their two state axes with
                        the fitted threshold drawn, which is what P3 tests, and
                        the run's own path across the same axes
  <run>_bilayer.png     area per lipid and thickness per chunk against section
                        5.8's windows

Idempotent and safe to run on a partial run: it draws whatever frames are in
the database, so it can be run during production and again at the end. Figures
for which there is no data are skipped with a line saying so rather than
written empty, because an axis with nothing on it looks like a result of zero.

Units are angstrom and nanoseconds throughout, as everywhere else in this
build.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The brand palette, read off flexappeal/static/brand.css rather than guessed,
# so a figure in the document and a chart on the page are the same colours.
INK = "#1a1a2e"
MUTED = "#5a6b80"
BORDER = "#dde4ed"
GREEN = "#00d084"
ORANGE = "#ff6900"
PURPLE = "#9b51e0"
AMBER = "#fcb900"
RED = "#d63638"
BLUE = "#2b6cb0"

def motif_palette() -> dict[str, str]:
    """Section 6.1's motif tokens, read out of the stylesheet that defines them.

    Not a copy. A figure in the document and a line on the page have to be the
    same colour for the same motif or the token system is pointless, and this
    file assigning its own colours by index was a fourth place for one fact to
    live. `brand.css` defines the variables; the `motifs` table says which one
    each motif has; this joins them.
    """
    import re

    css = (ROOT / "flexappeal" / "static" / "brand.css").read_text()
    return {name: value for name, value in
            re.findall(r"(--motif-\d+):\s*(#[0-9a-fA-F]{3,8})\s*;", css)}

# Section 5.8's windows, imported rather than restated so a figure cannot show
# a window the checks do not use.
sys.path.insert(0, str(ROOT))


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def shown(path: Path) -> str:
    """A path for the log, relative to the project where that makes sense.

    `relative_to` raises rather than giving up when the path is elsewhere, and
    `--out` is allowed to point anywhere, so the short form is an attempt and
    not an assumption.
    """
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def rows(connection, sql: str, params: tuple = ()) -> list[dict]:
    connection.row_factory = sqlite3.Row
    return [dict(r) for r in connection.execute(sql, params)]


def _style(axis) -> None:
    """The house look: no top or right spine, muted ticks, a light grid."""
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(BORDER)
    axis.tick_params(colors=MUTED, labelsize=9)
    axis.grid(True, color=BORDER, linewidth=0.6, alpha=0.8)
    axis.set_axisbelow(True)
    axis.xaxis.label.set_color(INK)
    axis.yaxis.label.set_color(INK)


def motif_figure(connection, run_id: str, out: Path) -> bool:
    """The per-motif score against time, which is what a row is for."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frames = rows(connection,
                  "SELECT DISTINCT motif_id, time_ns, s_score FROM motif_frames "
                  "WHERE run_id = ? AND s_score IS NOT NULL ORDER BY time_ns",
                  (run_id,))
    if not frames:
        log(f"  {run_id}: no scored frames, so no motif figure")
        return False

    series: dict[str, list[tuple[float, float]]] = {}
    for row in frames:
        series.setdefault(row["motif_id"], []).append(
            (row["time_ns"], row["s_score"]))

    palette = motif_palette()
    tokens = {row["motif_id"]: row["token"] for row in
              rows(connection, "SELECT motif_id, token FROM motifs")}

    figure, axis = plt.subplots(figsize=(8.0, 4.4), dpi=130)
    # The bands first, so the lines sit on top of them. Above zero is
    # active-like and below is inactive-like; the sign is the whole point of
    # the score, so the figure says it rather than the caption having to.
    axis.axhspan(0, 1, color=GREEN, alpha=0.05)
    axis.axhspan(-1, 0, color=RED, alpha=0.05)
    axis.axhline(0, color=MUTED, linewidth=1.0)

    for motif, points in sorted(series.items()):
        points.sort()
        colour = palette.get(tokens.get(motif) or "", MUTED)
        axis.plot([p[0] for p in points], [p[1] for p in points],
                  color=colour, linewidth=1.6, label=motif.split(":", 1)[-1])

    axis.set_xlabel("time (ns)")
    axis.set_ylabel("state score  s")
    axis.set_ylim(-1.05, 1.05)
    axis.set_title(f"{run_id}: per-motif state score", color=INK,
                   fontsize=11, loc="left")
    # Right of the axes, because a six-line legend inside a plot this size
    # covers the data it is labelling.
    axis.legend(frameon=False, fontsize=9, labelcolor=INK,
                loc="center left", bbox_to_anchor=(1.01, 0.5))
    _style(axis)
    axis.text(1.01, -0.02, "above zero: active-like\nbelow zero: inactive-like",
              transform=axis.transAxes, fontsize=8, color=MUTED,
              va="top", ha="left")
    figure.tight_layout()
    figure.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    log(f"  wrote {shown(out)} ({len(series)} motifs, "
        f"{len({f['time_ns'] for f in frames})} frames)")
    return True


def map_figure(client, system_id: str, out: Path,
               run_id: str | None = None) -> bool:
    """The reference structures on their own state axes, with the threshold.

    This is the figure behind section 1's corrected constants: the references
    are the evidence, and a reader can see the gap a threshold sits in rather
    than taking the number on trust. Section 6.2 also asks for the run's path
    across the same plane, which is drawn here when a run is named.

    Read from the explorer's own `/map/<system>.json` endpoint rather than from
    the database directly. The axes, the state preference and the published
    threshold line are decisions that endpoint already makes, and a figure that
    made them again could disagree with the page about what it is showing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    query = f"?run={run_id}" if run_id else ""
    response = client.get(f"/refstates/map/{system_id}.json{query}")
    if response.status_code != 200:
        log(f"  {system_id}: the map endpoint returned "
            f"{response.status_code}, so no map")
        return False
    data = response.get_json()
    points, path = data.get("points") or [], data.get("path") or []
    if not points:
        log(f"  {system_id}: the map endpoint returned no points")
        return False

    # Whatever states the manifest actually uses, rather than a fixed list: the
    # kinase labels them in and out, the receptor active and inactive.
    by_state: dict[str, list[dict]] = {}
    for point in points:
        if point.get("y") is None:
            continue
        by_state.setdefault(str(point.get("state")), []).append(point)
    colours = {"in": GREEN, "active": GREEN, "out": RED, "inactive": RED,
               "None": MUTED, "intermediate": AMBER}

    figure, axis = plt.subplots(figsize=(6.6, 5.4), dpi=130)
    for state, got in sorted(by_state.items()):
        axis.scatter([p["x"] for p in got], [p["y"] for p in got],
                     s=38, color=colours.get(state, MUTED), alpha=0.75,
                     linewidth=0.6, edgecolor="white", zorder=3,
                     label=f"{state} ({len(got)})")

    threshold = data.get("published_threshold")
    if threshold is not None:
        # On both axes, because the published criterion is a two-distance one
        # and applies the same cut to each. Drawing only the vertical line
        # hides the half of it that D010 is about: section 1.1's six
        # disagreements are all DFG-out structures sitting between 11.1 and
        # 11.3 A from the catalytic Lys, a tenth of an angstrom past the
        # horizontal line, and on the alphaC axis they are nowhere near it.
        axis.axvline(threshold, color=INK, linewidth=1.2, linestyle="--",
                     zorder=2,
                     label=f"published threshold {threshold:g} A, both axes")
        axis.axhline(threshold, color=INK, linewidth=1.2, linestyle="--",
                     zorder=2)

    if len(path) > 1:
        # Thin and faint, with the frames themselves marked over it and shaded
        # by time. Drawn at full weight the connecting line is a scribble that
        # covers the references it is supposed to be placed against, and it
        # shows no direction either.
        axis.plot([p["x"] for p in path], [p["y"] for p in path],
                  color=PURPLE, linewidth=0.7, alpha=0.35, zorder=4)
        times = [p.get("time_ns") or 0.0 for p in path]
        marks = axis.scatter([p["x"] for p in path], [p["y"] for p in path],
                             c=times, cmap="BuPu", s=16, zorder=5,
                             linewidth=0.3, edgecolor="white",
                             label=f"{run_id} ({len(path)} frames)")
        bar = figure.colorbar(marks, ax=axis, pad=0.02, fraction=0.045)
        bar.set_label("time (ns)", color=MUTED, fontsize=9)
        bar.ax.tick_params(colors=MUTED, labelsize=8)
        bar.outline.set_edgecolor(BORDER)
        axis.scatter([path[0]["x"]], [path[0]["y"]], s=90, marker="o",
                     facecolor="none", edgecolor=INK, linewidth=1.6,
                     zorder=6, label="start")

    axes_labels = data.get("axes") or ("", "")
    axis.set_xlabel(axes_labels[0])
    axis.set_ylabel(axes_labels[1])
    axis.set_title(f"{system_id}: reference state map", color=INK,
                   fontsize=11, loc="left")
    axis.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="best")
    _style(axis)
    figure.tight_layout()
    figure.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    counted = ", ".join(f"{len(v)} {k}" for k, v in sorted(by_state.items()))
    log(f"  wrote {shown(out)} ({counted}; "
        f"{len(path)} path points)")
    return True


def bilayer_figure(connection, run_id: str, out: Path) -> bool:
    """Area per lipid and thickness per chunk, against section 5.8's windows."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from flexappeal.refstates import membrane

    chunks = rows(connection,
                  "SELECT time_ns, apl, thickness, tilt FROM membrane "
                  "WHERE run_id = ? ORDER BY time_ns", (run_id,))
    if not chunks:
        log(f"  {run_id}: no membrane chunks, so no bilayer figure")
        return False

    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), dpi=130)
    times = [c["time_ns"] for c in chunks]

    windows = ((membrane.APL_EXPECTED_A2, membrane.APL_ACCEPTABLE_A2, "apl",
                "area per lipid (A^2)"),
               (membrane.THICKNESS_EXPECTED_A, None, "thickness",
                "P-P thickness (A)"))
    for axis, (expected, acceptable, key, label) in zip(axes, windows):
        if acceptable:
            axis.axhspan(*acceptable, color=AMBER, alpha=0.10,
                         label=f"acceptable {acceptable[0]:.0f} to "
                               f"{acceptable[1]:.0f}")
        axis.axhspan(*expected, color=GREEN, alpha=0.14,
                     label=f"expected {expected[0]:.0f} to {expected[1]:.0f}")
        values = [c[key] for c in chunks]
        # A single chunk has no line to draw, so mark it instead of drawing a
        # one-point series that renders as nothing.
        axis.plot(times, values, color=INK, linewidth=1.6, marker="o",
                  markersize=5, zorder=3)
        axis.set_xlabel("chunk start (ns)")
        axis.set_ylabel(label)
        axis.legend(frameon=False, fontsize=8, labelcolor=INK, loc="best")
        _style(axis)

    tilts = [c["tilt"] for c in chunks if c["tilt"] is not None]
    mean_tilt = (f", receptor tilt {sum(tilts) / len(tilts):.1f} degrees"
                 if tilts else "")
    figure.suptitle(f"{run_id}: bilayer health{mean_tilt}", color=INK,
                    fontsize=11, x=0.01, ha="left")
    figure.tight_layout()
    figure.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    log(f"  wrote {shown(out)} ({len(chunks)} chunks)")
    return True


# Which run's structure stands in for each system's card. Any run of that system
# would do, so the first alphabetically is chosen for stability rather than for
# meaning: the thumbnail is of the system, not of the run.
def thumbnail(analysis_dir: Path, run_id: str, system_id: str,
              out: Path) -> bool:
    """Section 6.2's card thumbnail, drawn from the system that was built.

    A Calpha trace and, where there is one, the bilayer phosphorus, projected
    onto the plane that shows them. Not a rendering: there is no GPU in this
    pipeline and a cartoon render would need one, so this is an honest line
    drawing of real coordinates rather than a picture of something else.

    The receptor is drawn looking along the membrane, x against z, because that
    is the view in which a 7TM bundle sitting in a bilayer looks like one. The
    kinase has no preferred direction, so it is drawn in its own widest plane.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mdtraj as md
    import numpy as np

    structure = analysis_dir / f"{run_id}.pdb"
    if not structure.is_file():
        log(f"  {system_id}: no viewer structure, so no thumbnail")
        return False

    frame = md.load(str(structure))
    ca = frame.topology.select("protein and name CA")
    if len(ca) < 3:
        log(f"  {system_id}: too few Calpha for a thumbnail")
        return False
    protein = frame.xyz[0, ca, :] * 10.0

    lipid_file = analysis_dir / f"{run_id}.lipids.pdb"
    lipids = None
    if lipid_file.is_file():
        from flexappeal.refstates import membrane

        loaded = md.load(str(lipid_file))
        lipids = loaded.xyz[0] * 10.0
        # Centred before drawing, and the protein moved by the same amount.
        # These are wrapped coordinates, and under a membrane barostat the
        # bilayer drifts until it straddles the periodic boundary: drawn raw,
        # the two leaflets land at opposite ends of the box and the receptor
        # appears to be sitting half out of the membrane. That is the fault
        # `centre_bilayer` exists for, and it returns its shift precisely so a
        # caller can put other atoms into the same frame.
        box_z = (float(loaded.unitcell_lengths[0][2]) * 10.0
                 if loaded.unitcell_lengths is not None else 0.0)
        if box_z:
            centred, shift = membrane.centre_bilayer(lipids[:, 2], box_z)
            lipids = lipids.copy()
            lipids[:, 2] = centred
            protein = protein.copy()
            protein[:, 2] = np.mod(protein[:, 2] + shift, box_z)
            # The protein can still land the long way round after the shift:
            # anything more than half a box from the bilayer centre is nearer
            # the other way, which is what the periodicity means.
            centre = float(np.median(centred))
            far = protein[:, 2] - centre > 0.5 * box_z
            protein[far, 2] -= box_z
            near = centre - protein[:, 2] > 0.5 * box_z
            protein[near, 2] += box_z

    if lipids is not None:
        # Looking along the bilayer: x across, z through.
        axes = (0, 2)
    else:
        # The two widest directions of the molecule, so it fills the frame.
        spread = protein.std(axis=0)
        axes = tuple(int(i) for i in np.argsort(spread)[::-1][:2])

    figure, axis = plt.subplots(figsize=(2.6, 2.0), dpi=150)
    if lipids is not None:
        axis.scatter(lipids[:, axes[0]], lipids[:, axes[1]], s=5,
                     color=AMBER, alpha=0.55, linewidth=0, zorder=1)
    # The backbone in order, so it reads as a chain rather than a cloud.
    axis.plot(protein[:, axes[0]], protein[:, axes[1]],
              color=PURPLE, linewidth=0.9, alpha=0.9, zorder=2)
    axis.set_aspect("equal")
    axis.axis("off")
    if lipids is not None:
        # Framed on the receptor. The bilayer is a periodic slab several times
        # wider than the protein, so letting it set the extent puts the subject
        # in a corner of its own thumbnail.
        centre_x = float(np.mean(protein[:, axes[0]]))
        half = max(float(np.ptp(protein[:, axes[0]])),
                   float(np.ptp(protein[:, axes[1]]))) * 0.75
        axis.set_xlim(centre_x - half, centre_x + half)
    figure.tight_layout(pad=0.1)
    figure.savefig(out, bbox_inches="tight", facecolor="white",
                   transparent=False)
    plt.close(figure)
    log(f"  wrote {shown(out)} ({len(ca)} Calpha"
        + (f", {len(lipids)} lipid phosphorus" if lipids is not None else "")
        + ")")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=ROOT / "data" / "refstates.db")
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "figures")
    parser.add_argument("--only", default=None,
                        help="one run id, for a quick check")
    args = parser.parse_args()

    if not args.db.is_file():
        log(f"no database at {args.db}")
        return 1
    args.out.mkdir(parents=True, exist_ok=True)

    written = 0
    # The state map comes from the explorer's own endpoint, so the app is built
    # here and pointed at the same database the figures are drawn from.
    from flexappeal.webapp import create_app

    app = create_app()
    app.config["REFSTATES_DB"] = str(args.db)
    client = app.test_client()

    with sqlite3.connect(f"file:{args.db}?mode=ro", uri=True) as connection:
        runs = rows(connection, "SELECT run_id, system_id FROM runs ORDER BY run_id")
        if args.only:
            runs = [r for r in runs if r["run_id"] == args.only]
        if not runs:
            log("no runs to draw")
            return 1

        for run in runs:
            log(f"{run['run_id']}:")
            written += motif_figure(connection, run["run_id"],
                                    args.out / f"{run['run_id']}_motifs.png")
            written += bilayer_figure(connection, run["run_id"],
                                      args.out / f"{run['run_id']}_bilayer.png")

        # One map per system, with the first of that system's runs drawn over
        # it. A map per run would repeat the same 116 references four times.
        for system in sorted({r["system_id"] for r in runs}):
            log(f"{system}:")
            first = next(r["run_id"] for r in runs if r["system_id"] == system)
            written += map_figure(client, system,
                                  args.out / f"{system}_map.png", run_id=first)
            written += thumbnail(ROOT / "data" / "analysis", first, system,
                                 args.out / f"{system}_thumb.png")

    log(f"{written} figures in {shown(args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
