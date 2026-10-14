"""The reference-state panels, in the Analysis tab's house style.

These used to be a tab of their own: five pages that read a SQLite database and
built their charts in the browser out of four JSON endpoints. The runs are now
rows in the Runs library like any other run, so the measurements that are
particular to them became what they always should have been, which is panels at
the end of a run page.

Everything here returns the same shape the rest of the app's panels do, so the
run page renders them with the macro it already has:

    {"id", "title", "blurb", "wide", "kind", "figure" | "table"}

Two kinds, because two of these are tables rather than charts and always were.
`kind` is explicit on every panel rather than inferred from which key is
present: a template that switches on the absence of a key fails silently the
first time a builder returns an empty figure.

Where the numbers come from
---------------------------
One member of the results file, `refstates.json`, written by
`tools/pack_run_fxa.py`. A run in the library is one file, so these panels read
what is inside it rather than a database beside it: the deployed site serves
four results files and nothing else.

Colour
------
Section 6.1 fixes the motif and state colours as CSS custom properties, and
CLAUDE.md says not to pick a colour in a chart. These figures are built on the
server, where there is no computed style to read, so the tokens are parsed out
of `brand.css` itself rather than restated here. One definition, in the file the
browser also reads, and a missing token raises rather than quietly defaulting to
something a reader would take for a measurement.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .. import plots
from . import store
from .state import DECISIVE_S

BRAND_CSS = Path(__file__).resolve().parents[1] / "static" / "brand.css"

# The DSSP codes mdtraj's simplified set does not use, mapped onto the four
# categories the Analysis tab's own secondary-structure panel draws. The eight
# state codes come from the reference-state analysis, which keeps them: a run
# page that showed helix, sheet and coil in one place and H, G, I, E, B, T, S, C
# in another would be two vocabularies for one measurement.
_DSSP_CATEGORY = {
    "H": 0, "G": 0, "I": 0,          # helix, including 3-10 and pi
    "E": 1, "B": 1,                  # sheet, including an isolated bridge
    "T": 2, "S": 2, "C": 2, " ": 2,  # turn, bend and coil all read as coil
}


def _tokens() -> dict[str, str]:
    """The `--motif-*` and `--state-*` values, read from `brand.css`.

    Cached on the function, because the file does not change while the process
    runs and every panel on every request would otherwise re-read it.
    """
    cached = getattr(_tokens, "_cache", None)
    if cached is not None:
        return cached
    text = BRAND_CSS.read_text(encoding="utf-8")
    found = dict(re.findall(
        r"(--(?:motif|state)-[a-z0-9]+)\s*:\s*(#[0-9a-fA-F]{3,8})", text))
    missing = [name for name in
               ("--state-active", "--state-inactive", "--state-intermediate",
                "--state-conflict", "--motif-1", "--motif-2", "--motif-3",
                "--motif-4", "--motif-5", "--motif-6")
               if name not in found]
    if missing:
        raise RuntimeError(
            f"brand.css is missing the tokens section 6.1 fixes: {missing}. "
            f"These panels read the colours from there rather than choosing "
            f"their own, so there is nothing sensible to fall back to.")
    _tokens._cache = found                                   # type: ignore[attr-defined]
    return found


# What each catalogue calls a state, mapped onto the two section 6.1 has
# tokens for. KLIFS says "in" and "out" of the DFG motif and GPCRdb says
# "active" and "inactive", and the state map coloured a kinase cloud entirely
# grey because neither of KLIFS's words is a token name: both clusters came out
# the same colour on the one panel whose whole point is that they separate.
STATE_SYNONYMS = {"in": "active", "dfg-in": "active", "active": "active",
                  "out": "inactive", "dfg-out": "inactive",
                  "inactive": "inactive", "intermediate": "intermediate",
                  "conflict": "conflict"}


def state_colour(label: str) -> str:
    """The token for a state label, so a badge and a chart agree."""
    name = STATE_SYNONYMS.get(str(label).strip().lower(), str(label).lower())
    return _tokens().get(f"--state-{name}", plots.INK_MUTED)


def motif_colour(kind: str, motif: str) -> str:
    """Section 6.1's motif colour for one motif of one family."""
    token = (store.MOTIF_TOKENS.get(kind) or {}).get(motif)
    return _tokens().get(token or "", plots.CATEGORICAL[0])


def _decisive_bands(y_title: str = "State score s") -> dict[str, Any]:
    """The bands and labels every state-score chart is read against.

    Drawn at `state.DECISIVE_S` rather than at zero, which is 4ad: a score of
    0.05 is not "slightly active", it is undecided, and a line at zero invites
    exactly that reading. Section 10 requires every threshold to be drawn and
    labelled, so the label is part of the figure rather than left to the blurb.
    """
    return {
        "yaxis": {
            "title": {"text": y_title,
                      "font": {"size": 12, "color": plots.INK_MUTED}},
            "gridcolor": plots.GRID, "zeroline": True,
            "zerolinecolor": plots.AXIS, "linecolor": plots.AXIS,
            "ticks": "outside", "tickcolor": plots.AXIS,
            "tickfont": {"size": 11}, "range": [-1.05, 1.05],
        },
        "shapes": [
            {"type": "rect", "xref": "paper", "x0": 0, "x1": 1,
             "yref": "y", "y0": DECISIVE_S, "y1": 1.05,
             "fillcolor": state_colour("active"), "opacity": 0.06,
             "line": {"width": 0}, "layer": "below"},
            {"type": "rect", "xref": "paper", "x0": 0, "x1": 1,
             "yref": "y", "y0": -1.05, "y1": -DECISIVE_S,
             "fillcolor": state_colour("inactive"), "opacity": 0.06,
             "line": {"width": 0}, "layer": "below"},
        ],
        # On a panel background rather than bare. A label at the edge of the
        # lower band sits exactly where an inactive run's trace does, and
        # `b2ar_carazolol` runs at about -0.3 for its whole 50 ns: without
        # this the threshold a reader is meant to read the trace against is
        # printed through the trace.
        "annotations": [
            {"xref": "paper", "x": 1, "xanchor": "right", "yref": "y",
             "y": DECISIVE_S, "yanchor": "bottom", "showarrow": False,
             "text": f"active beyond s = +{DECISIVE_S}",
             "bgcolor": "rgba(255,255,255,0.82)", "borderpad": 2,
             "font": {"size": 10, "color": state_colour("active")}},
            {"xref": "paper", "x": 1, "xanchor": "right", "yref": "y",
             "y": -DECISIVE_S, "yanchor": "top", "showarrow": False,
             "text": f"inactive beyond s = -{DECISIVE_S}",
             "bgcolor": "rgba(255,255,255,0.82)", "borderpad": 2,
             "font": {"size": 10, "color": state_colour("inactive")}},
        ],
    }


# ===========================================================================
#  The charts
# ===========================================================================


def state_score(payload: dict[str, Any]) -> dict[str, Any] | None:
    """The whole-protein score over time, raw and smoothed.

    Both, because they answer different questions and the smoothed one alone
    would hide how much of the scatter the smoothing removed. The window is the
    analysis's own, stated in the blurb rather than assumed by a reader.
    """
    raw = payload.get("whole_s_raw") or []
    smoothed = payload.get("whole_s_smoothed") or []
    frames = payload.get("frames") or []
    if not raw or not frames:
        return None
    time = [float(frame.get("time_ns", index))
            for index, frame in enumerate(frames)]

    stride = plots.wire_stride(len(time))
    time_thin = time[::stride]
    traces = [{
        "type": "scatter", "mode": "lines", "name": "per frame",
        "x": time_thin, "y": plots.finite(raw[::stride]),
        "line": {"color": plots.INK_MUTED, "width": 1},
        "opacity": 0.45,
        "hovertemplate": "%{y:+.3f} per frame<extra></extra>",
    }]
    if smoothed:
        traces.append({
            "type": "scatter", "mode": "lines", "name": "smoothed",
            "x": time_thin, "y": plots.finite(smoothed[::stride]),
            "line": {"color": plots.CATEGORICAL[0], "width": 2},
            "hovertemplate": "%{y:+.3f} smoothed<extra></extra>",
        })
    return plots.figure(traces, plots.layout(
        "State score over time", "Time (ns)", "State score s",
        showlegend=True,
        legend={"orientation": "h", "y": -0.2, "x": 0,
                "font": {"size": 11, "color": plots.INK_MUTED}},
        margin={"l": 62, "r": 18, "t": 44, "b": 70},
        **_decisive_bands()))


def motif_timelines(payload: dict[str, Any]) -> dict[str, Any] | None:
    """One line per motif, in section 6.1's colours.

    The specified chart was a small multiple per motif. This is one panel with a
    line each, which is 4u: the question a reader has is whether the motifs
    agree with each other, and six panels side by side answer it worse than six
    lines on one pair of axes.
    """
    frames = payload.get("frames") or []
    if not frames:
        return None
    kind = payload.get("kind") or "kinase"
    names = [motif["motif"] for motif in (frames[0].get("motifs") or [])]
    if not names:
        return None

    time = [float(frame.get("time_ns", index))
            for index, frame in enumerate(frames)]
    stride = plots.wire_stride(len(time))
    series: dict[str, list[float | None]] = {name: [] for name in names}
    for frame in frames:
        found = {motif["motif"]: motif.get("s")
                 for motif in (frame.get("motifs") or [])}
        for name in names:
            series[name].append(found.get(name))

    traces = []
    for name in sorted(names):
        traces.append({
            "type": "scatter", "mode": "lines", "name": name,
            "x": time[::stride], "y": plots.finite(series[name][::stride]),
            "line": {"color": motif_colour(kind, name), "width": 2},
            "hovertemplate": f"{name} %{{y:+.3f}}<extra></extra>",
        })
    return plots.figure(traces, plots.layout(
        "Motif timelines", "Time (ns)", "State score s",
        showlegend=True,
        legend={"orientation": "h", "y": -0.2, "x": 0,
                "font": {"size": 11, "color": plots.INK_MUTED}},
        margin={"l": 62, "r": 18, "t": 44, "b": 80},
        **_decisive_bands()))


def _motif_lanes(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Each motif's residue span, packed into lanes so none hides another.

    The spans overlap: NPxxY sits inside TM7 on this receptor, and a band drawn
    over another band hides it. Each motif takes the first lane whose previous
    band ended before this one starts, which is the usual interval packing and
    needs two lanes here.
    """
    kind = payload.get("kind") or "kinase"
    spans = []
    for motif, numbers in (payload.get("motif_residues") or {}).items():
        if numbers:
            spans.append((min(numbers), max(numbers), motif))
    lane_ends: list[int] = []
    packed = []
    for first, last, motif in sorted(spans):
        lane = next((index for index, end in enumerate(lane_ends)
                     if end < first), len(lane_ends))
        if lane == len(lane_ends):
            lane_ends.append(last)
        else:
            lane_ends[lane] = last
        packed.append({"motif": motif, "first": first, "last": last,
                       "lane": lane, "colour": motif_colour(kind, motif)})
    return packed


def dssp_ribbon(payload: dict[str, Any]) -> dict[str, Any] | None:
    """Secondary structure, with the motifs marked and the references beside it.

    The same heatmap the Analysis tab draws, in the same four colours, with the
    two things section 6.2 asks for added: a band per motif down the residue
    axis, so a reader can see which rows belong to which, and two columns
    holding the modal secondary structure of the active and inactive reference
    ensembles at that residue. The columns are what make the panel answer a
    question about state rather than only about this run.

    The eight DSSP codes the reference-state analysis keeps are folded onto the
    four categories the Analysis tab's own panel draws. That is deliberate and
    it is the one thing here that is not a port: this page already has a
    secondary-structure panel, and showing alpha, 3-10, pi, strand, bridge,
    turn, bend and coil in one chart and helix, sheet, coil and none in the one
    above it is two vocabularies for a single measurement.
    """
    block = payload.get("dssp") or {}
    frames = block.get("frames") or []
    residues = block.get("residues") or []
    if not frames or not residues:
        return None

    times = [float(frame.get("time_ns", index))
             for index, frame in enumerate(payload.get("frames") or [])]
    if len(times) != len(frames):
        times = list(range(len(frames)))

    # Half the columns the Analysis tab's own secondary-structure panel keeps.
    # That panel is on the same page at full resolution, and this one is the
    # annotated reading of the same measurement: what it adds is read along a
    # row rather than along time. At 400 columns the two heatmaps came to
    # 1.7 MB of the page between them, nearly all of it the per-cell state
    # names the hover needs. 200 columns of a 50 ns run is a column every
    # 0.25 ns, finer than any secondary-structure change worth reading off it.
    stride = plots.wire_stride(len(frames), plots.WIRE_HEATMAP_COLUMNS // 2)
    kept = frames[::stride]
    time_axis = times[::stride]

    states = plots.dssp_states
    n = len(states)
    colourscale: list[list[Any]] = []
    for value, _, colour in states:
        colourscale.append([value / n, colour])
        colourscale.append([(value + 1) / n, colour])
    labels = {value: label for value, label, _ in states}

    def category(code: str | None) -> int | None:
        if code is None or code in ("-", ""):
            return None
        return _DSSP_CATEGORY.get(code, 2)

    # Residues down y, time across x, which is the orientation the other DSSP
    # panel uses. Transposed here rather than in the browser.
    grid = [[category(row[index] if index < len(row) else None)
             for row in kept] for index in range(len(residues))]
    text = [[labels.get(cell, "no assignment") if cell is not None
             else "no assignment" for cell in row] for row in grid]

    traces: list[dict[str, Any]] = [{
        "type": "heatmap", "z": grid, "x": time_axis, "y": list(residues),
        "colorscale": colourscale, "zmin": 0, "zmax": n, "showscale": False,
        "text": text,
        "hovertemplate": "Residue %{y}<br>%{x:.2f} ns<br>%{text}<extra></extra>",
        "xgap": 0, "ygap": 0,
    }]

    span = (time_axis[-1] - time_axis[0]) if len(time_axis) > 1 else 1.0
    # A heatmap takes each cell's width from the gap between neighbouring x
    # values, so a column given one x draws at Plotly's default width of one
    # unit: a 1 ns sliver against a 50 ns axis. Each consensus column is
    # therefore two x values holding the same value, which gives the gap
    # explicitly and keeps the two columns the same width on any run length.
    width = max(span * 0.022, 0.2)
    consensus = sse_consensus(payload, residues)
    shown: list[tuple[str, float]] = []
    for offset, state in enumerate(("active", "inactive")):
        column = consensus.get(state) or []
        if not any(column):
            continue
        centre = time_axis[-1] + width * (1.3 + offset * 1.45)
        cells = [[category((item or {}).get("code")) if item else None] * 2
                 for item in column]
        hover = []
        for item in column:
            if not item:
                hover.append([""] * 2)
                continue
            agreement = item.get("agreement")
            said = (f"{state} references: "
                    f"{labels.get(category(item.get('code')), 'no assignment')}"
                    + (f", {agreement:.0%} agree" if agreement is not None else ""))
            hover.append([said] * 2)
        traces.append({
            "type": "heatmap", "z": cells,
            "x": [centre - width / 2, centre + width / 2],
            "y": list(residues),
            "colorscale": colourscale, "zmin": 0, "zmax": n,
            "showscale": False, "text": hover,
            "hovertemplate": "Residue %{y}<br>%{text}<extra></extra>",
            "xgap": 0, "ygap": 0,
        })
        shown.append((state, centre))

    # The motif bands, left of t = 0 in the same residue coordinates, so the
    # chart reads left to right as motif, then trajectory, then the two
    # reference columns.
    lanes = _motif_lanes(payload)
    n_lanes = (max(item["lane"] for item in lanes) + 1) if lanes else 0
    shapes = []
    annotations = []
    for item in lanes:
        x1 = time_axis[0] - width * (0.4 + 0.75 * item["lane"])
        x0 = x1 - width * 0.6
        shapes.append({
            "type": "rect", "xref": "x", "yref": "y",
            "x0": x0, "x1": x1,
            "y0": item["first"] - 0.5, "y1": item["last"] + 0.5,
            "fillcolor": item["colour"], "line": {"width": 0}, "layer": "above",
        })
        # A short motif has no room for a label turned on its side: NPxxY is
        # five residues and the word is longer than the bar is tall. Those are
        # written beside the stack in the motif's own colour, so a band is
        # never a colour with no key.
        tall = (item["last"] - item["first"]) >= 12
        annotations.append({
            "xref": "x", "yref": "y",
            "x": (x0 + x1) / 2 if tall else x0 - width * 0.15,
            "y": (item["first"] + item["last"]) / 2,
            "text": item["motif"], "showarrow": False,
            "textangle": -90 if tall else 0,
            "xanchor": "center" if tall else "right",
            "font": {"size": 10,
                     "color": "#ffffff" if tall else item["colour"]},
        })
    for state, centre in shown:
        # Turned on its side, because the column is about 25 px wide on a
        # desktop panel and the word "inactive" is nearer 45: written flat, the
        # two labels printed over each other and read as one word.
        annotations.append({
            "xref": "x", "yref": "paper", "x": centre, "y": 1.01,
            "yanchor": "bottom", "xanchor": "center", "showarrow": False,
            "textangle": -90, "text": state,
            "font": {"size": 10, "color": state_colour(state)},
        })
    if shown:
        # The columns are not frames, so the boundary is drawn rather than left
        # to read as the last few nanoseconds of the run.
        shapes.append({
            "type": "line", "yref": "paper", "y0": 0, "y1": 1,
            "xref": "x", "x0": time_axis[-1] + width * 0.45,
            "x1": time_axis[-1] + width * 0.45,
            "line": {"color": plots.INK, "width": 1.5},
        })

    legend = [{
        "type": "scatter", "mode": "markers", "name": label,
        "x": [None], "y": [None], "showlegend": True,
        "marker": {"color": colour, "size": 10, "symbol": "square"},
        "hoverinfo": "skip",
    } for _, label, colour in states]

    # Set explicitly, because Plotly's autorange looks at traces and not at
    # shapes: the motif bands sit left of every data column, so left to itself
    # the axis would start at the first frame and clip them entirely.
    left = time_axis[0] - width * (0.6 + 0.75 * max(n_lanes, 1)) - width
    right = ((shown[-1][1] + width) if shown
             else time_axis[-1] + span * 0.01)

    return plots.figure(traces + legend, plots.layout(
        "Secondary structure, motifs and the reference consensus",
        "Time (ns)", "Residue",
        showlegend=True,
        legend={"orientation": "h", "y": -0.16, "x": 0,
                "font": {"size": 11, "color": plots.INK_MUTED}},
        hovermode="closest",
        margin={"l": 72, "r": 18, "t": 74, "b": 70},
        shapes=shapes, annotations=annotations,
        xaxis={"title": {"text": "Time (ns)",
                         "font": {"size": 12, "color": plots.INK_MUTED}},
               "gridcolor": plots.GRID, "zeroline": False,
               "linecolor": plots.AXIS, "ticks": "outside",
               "tickcolor": plots.AXIS, "tickfont": {"size": 11},
               "range": [left, right]},
        yaxis={"title": {"text": "Residue",
                         "font": {"size": 12, "color": plots.INK_MUTED}},
               "gridcolor": plots.GRID, "zeroline": False,
               "linecolor": plots.AXIS, "ticks": "outside",
               "tickcolor": plots.AXIS, "tickfont": {"size": 11},
               "autorange": "reversed"}))


def state_map(payload: dict[str, Any]) -> dict[str, Any] | None:
    """The reference cloud on its published plane, with this run's path on it.

    For the kinase this is the Dunbrack plane, the two DFG-Phe distances, which
    is the plane the published criterion is drawn on and where D010's threshold
    failure is visible. For the receptor the x axis is the TM3-TM6 cleft, the
    axis the Dror 2011 maps use.

    The references are coloured by the catalogue's own label rather than by this
    build's, so the cloud is the thing being agreed with rather than a picture
    of this build agreeing with itself.
    """
    system = payload.get("system_id") or ""
    references = payload.get("reference_structures") or []
    if not references:
        return None

    if system == "abl1":
        keys = ("dfg_d_alphac", "dfg_d_lys")
        # The unit symbol and the Greek letters, not ASCII transliterations.
        # The same convention as every other label on the page, which is what
        # the guard in the suite is for: the labels this panel replaced were
        # written "(A)" with "Calpha" spelled out, so one plot carried a
        # different convention from the panel beside it.
        axes = ("DFG-Phe to \u03b1C-Glu+4 C\u03b1 (\u00c5)",
                "DFG-Phe to catalytic Lys C\u03b1 (\u00c5)")
    else:
        keys = ("tm3_tm6", "npxxy_to_inactive")
        axes = ("TM3-TM6 C\u03b1 distance, 3.50 to 6.34 (\u00c5)",
                "NPxxY RMSD to the nearest inactive reference (\u00c5)")

    clouds: dict[str, dict[str, list]] = {}
    for reference in references:
        metrics = reference.get("metrics") or {}
        x, y = metrics.get(keys[0]), metrics.get(keys[1])
        if x is None or y is None:
            continue
        state = (reference.get("source_state") or reference.get("our_state")
                 or "intermediate")
        bucket = clouds.setdefault(state, {"x": [], "y": [], "pdb": []})
        bucket["x"].append(x)
        bucket["y"].append(y)
        bucket["pdb"].append(reference.get("pdb", ""))
    if not clouds:
        return None

    traces = []
    for state, bucket in sorted(clouds.items()):
        # Shape as well as colour: section 6.1 asks that no state be readable
        # by colour alone, and this is the one panel where two clouds have to
        # be told apart at a glance.
        resolved = STATE_SYNONYMS.get(state.lower(), state.lower())
        traces.append({
            "type": "scatter", "mode": "markers",
            "name": f"{state} references",
            "x": bucket["x"], "y": bucket["y"], "text": bucket["pdb"],
            "marker": {"color": state_colour(state), "size": 9,
                       "symbol": {"active": "circle",
                                  "inactive": "square"}.get(resolved, "diamond"),
                       "line": {"width": 1, "color": "#ffffff"}},
            "hovertemplate": "%{text}<br>%{x:.2f}, %{y:.2f} \u00c5"
                             f"<br>{state} in the catalogue<extra></extra>",
        })

    # This run's path across the same plane. Each point carries the frame the
    # viewer should show, which is what lets a click on the chart move it there.
    #
    # The viewer's frame, not the analysis's. They are the same for these four
    # rows, where both are every frame of a 500-frame trajectory, and they are
    # not the same in general: a results file is packed at a stride, so a
    # 2,000-frame run carries 500. Scaled here rather than in the click
    # handler, so what the point says is what the slider wants and the handler
    # has nothing to get wrong.
    packed = int(payload.get("frames_packed") or 0)
    analysed = len(payload.get("frames") or [])

    def viewer_frame(index: int) -> int:
        if packed < 2 or analysed < 2:
            return 0
        return round(index / (analysed - 1) * (packed - 1))

    def coordinates(frame: dict[str, Any]) -> tuple[Any, Any]:
        """Where this frame sits on the plane the references are drawn on.

        The frame's own measurements are nested and the reference manifest's
        are flat, which is why this is a function rather than two `get`s: the
        first version read the reference keys out of a frame, found nothing,
        and drew a reference cloud with no trajectory across it. The panel
        looked complete, which is what made it worth a function with a name.
        """
        found = frame.get("measurements") or {}
        if system == "abl1":
            plane = found.get("dfg_plane") or {}
            return plane.get("d_alphac_a"), plane.get("d_lys_a")
        x = (found.get("micro_switches") or {}).get("tm3_tm6_a")
        # The NPxxY distance to the nearest inactive reference is already
        # computed as that motif's score, so it is read from there rather than
        # measured twice.
        y = next((motif.get("d_inact_a") for motif in frame.get("motifs") or []
                  if motif.get("motif") == "NPxxY"), None)
        return x, y

    path_x, path_y, frame_index, hover = [], [], [], []
    for index, frame in enumerate(payload.get("frames") or []):
        x, y = coordinates(frame)
        if x is None or y is None:
            continue
        path_x.append(x)
        path_y.append(y)
        frame_index.append(viewer_frame(index))
        hover.append(f"frame {index}, {frame.get('time_ns', 0):.1f} ns")
    if path_x:
        stride = plots.wire_stride(len(path_x), 300)
        # Coloured by time, so the direction of travel is visible rather than
        # implied. Five hundred points joined by a line is a scribble, and a
        # reader looking at one cannot tell where it started: the question this
        # panel is asked is whether a run left one cluster for the other, and
        # that is a question about direction. The ramp is the house sequential
        # one, so it means the same thing here as on a contact map.
        times = [float((payload.get("frames") or [])[index].get("time_ns", 0.0))
                 for index in range(len(payload.get("frames") or []))]
        along = [times[index] if index < len(times) else 0.0
                 for index in frame_index][::stride]
        traces.append({
            "type": "scatter", "mode": "lines+markers", "name": "this run",
            "x": path_x[::stride], "y": path_y[::stride],
            "text": hover[::stride],
            "customdata": frame_index[::stride],
            "line": {"color": "rgba(26,26,46,0.25)", "width": 1},
            "marker": {"size": 5, "color": along,
                       "colorscale": [list(stop) for stop in plots.SEQUENTIAL],
                       "showscale": True,
                       "colorbar": {"title": {"text": "ns",
                                              "font": {"size": 11}},
                                    "thickness": 10, "len": 0.6,
                                    "outlinewidth": 0,
                                    "tickfont": {"size": 10}}},
            "hovertemplate": "%{text}<br>%{x:.2f}, %{y:.2f} \u00c5"
                             "<extra></extra>",
        })

    # The line the state is called on, drawn and labelled, which is section
    # 10's rule for every threshold on every plot. The two systems get it from
    # different places and the label says which: the kinase's is published and
    # the receptor's is fitted on this very reference set, and a reader is
    # entitled to know that one of them was chosen to fit the cloud it is
    # drawn over.
    shapes: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = []
    if system == "abl1":
        from .kinase import DFG_PUBLISHED_A

        value, said = DFG_PUBLISHED_A, f"published {DFG_PUBLISHED_A:.1f} Å"
        axis = "x"
    else:
        from .gpcr import fit_threshold

        by_state: dict[str, list[float]] = {"active": [], "inactive": []}
        for reference in references:
            state = (reference.get("source_state") or "").lower()
            distance = (reference.get("metrics") or {}).get("tm3_tm6")
            if state in by_state and distance is not None:
                by_state[state].append(float(distance))
        fitted = fit_threshold(by_state["active"], by_state["inactive"])
        value = fitted.get("threshold")
        accuracy = fitted.get("accuracy")
        said = (f"fitted {value:.2f} Å"
                + (f", {accuracy:.0%} leave-one-out" if accuracy else "")
                if value is not None else "")
        axis = "x"
    if value is not None:
        shapes.append(
            {"type": "line", "xref": axis, "x0": value, "x1": value,
             "yref": "paper", "y0": 0, "y1": 1,
             "line": {"color": state_colour("inactive"), "width": 1.5,
                      "dash": "dash"}})
        annotations.append(
            {"xref": axis, "x": value, "yref": "paper", "y": 1,
             "yanchor": "bottom", "xanchor": "left", "showarrow": False,
             "text": said, "font": {"size": 10, "color": plots.INK_MUTED}})

    return plots.figure(traces, plots.layout(
        "Reference state map", axes[0], axes[1],
        showlegend=True,
        legend={"orientation": "h", "y": -0.22, "x": 0,
                "font": {"size": 11, "color": plots.INK_MUTED}},
        hovermode="closest",
        margin={"l": 70, "r": 18, "t": 52, "b": 86},
        shapes=shapes, annotations=annotations))


def bilayer(payload: dict[str, Any]) -> dict[str, Any] | None:
    """Area per lipid and thickness per chunk, against section 5.8's two windows.

    The Analysis tab already plots both per frame. This panel is the question
    the bilayer QC asks instead: whether each 10 ns chunk sat inside the
    expected window, inside the wider acceptable one, or outside both, which is
    what G6 is decided on. Both windows are drawn, because a chunk between them
    is a real and reportable result rather than a failure.
    """
    chunks = payload.get("membrane_chunks") or []
    scored = [chunk for chunk in chunks if chunk.get("apl") is not None]
    if not scored:
        return None

    from .membrane import (APL_ACCEPTABLE_A2, APL_EXPECTED_A2,
                           THICKNESS_EXPECTED_A)

    time = [float(chunk.get("time_ns") or 0.0) for chunk in scored]
    apl = [float(chunk["apl"]) for chunk in scored]
    thickness = [chunk.get("thickness") for chunk in scored]

    traces = [{
        "type": "scatter", "mode": "lines+markers", "name": "area per lipid",
        "x": time, "y": apl,
        "line": {"color": plots.CATEGORICAL[0], "width": 2},
        "marker": {"size": 8, "color": [
            state_colour("active") if APL_EXPECTED_A2[0] <= value <= APL_EXPECTED_A2[1]
            else state_colour("intermediate") for value in apl]},
        "hovertemplate": "%{x:.0f} ns<br>%{y:.2f} \u00c5\u00b2<extra></extra>",
    }]
    if any(value is not None for value in thickness):
        traces.append({
            "type": "scatter", "mode": "lines+markers", "name": "thickness",
            "x": time, "y": plots.finite([v if v is not None else float("nan")
                                          for v in thickness]),
            "yaxis": "y2",
            "line": {"color": plots.CATEGORICAL[1], "width": 2, "dash": "dot"},
            "marker": {"size": 6},
            "hovertemplate": "%{x:.0f} ns<br>%{y:.2f} \u00c5 thick<extra></extra>",
        })

    shapes = [
        {"type": "rect", "xref": "paper", "x0": 0, "x1": 1, "yref": "y",
         "y0": APL_EXPECTED_A2[0], "y1": APL_EXPECTED_A2[1],
         "fillcolor": state_colour("active"), "opacity": 0.07,
         "line": {"width": 0}, "layer": "below"},
        {"type": "line", "xref": "paper", "x0": 0, "x1": 1, "yref": "y",
         "y0": APL_ACCEPTABLE_A2[0], "y1": APL_ACCEPTABLE_A2[0],
         "line": {"color": plots.AXIS, "width": 1, "dash": "dash"}},
        {"type": "line", "xref": "paper", "x0": 0, "x1": 1, "yref": "y",
         "y0": APL_ACCEPTABLE_A2[1], "y1": APL_ACCEPTABLE_A2[1],
         "line": {"color": plots.AXIS, "width": 1, "dash": "dash"}},
    ]
    annotations = [
        {"xref": "paper", "x": 0, "xanchor": "left", "yref": "y",
         "y": APL_EXPECTED_A2[1], "yanchor": "bottom", "showarrow": False,
         "text": f"expected {APL_EXPECTED_A2[0]:.0f} to {APL_EXPECTED_A2[1]:.0f} \u00c5\u00b2",
         "font": {"size": 10, "color": state_colour("active")}},
        {"xref": "paper", "x": 0, "xanchor": "left", "yref": "y",
         "y": APL_ACCEPTABLE_A2[1], "yanchor": "bottom", "showarrow": False,
         "text": f"acceptable to {APL_ACCEPTABLE_A2[1]:.0f}",
         "font": {"size": 10, "color": plots.INK_MUTED}},
    ]

    return plots.figure(traces, plots.layout(
        "Bilayer quality, chunk by chunk", "Time (ns)",
        "Area per lipid (\u00c5\u00b2)",
        showlegend=True,
        legend={"orientation": "h", "y": -0.22, "x": 0,
                "font": {"size": 11, "color": plots.INK_MUTED}},
        margin={"l": 70, "r": 62, "t": 44, "b": 80},
        shapes=shapes, annotations=annotations,
        yaxis2={"title": {"text": f"Thickness (\u00c5), expected "
                                  f"{THICKNESS_EXPECTED_A[0]:.0f} to "
                                  f"{THICKNESS_EXPECTED_A[1]:.0f}",
                          "font": {"size": 11, "color": plots.INK_MUTED}},
                "overlaying": "y", "side": "right",
                "gridcolor": "rgba(0,0,0,0)", "zeroline": False,
                "linecolor": plots.AXIS, "tickfont": {"size": 10}}))


# ===========================================================================
#  The tables
# ===========================================================================


def sse_consensus(payload: dict[str, Any],
                  residues: list[int]) -> dict[str, list]:
    """The reference ensembles' secondary structure, on the ribbon's own axis.

    Section 6.2 asks the ribbon for an "active/inactive reference consensus row
    above it", so a reader can see whether this run's secondary structure looks
    like the active references or the inactive ones rather than only what it is.

    `report.sse_profiles` holds a count of DSSP codes per state, per motif, per
    position within that motif. Two things have to happen to put it on the
    ribbon: the modal code has to be picked, and the motif position resolved to
    a residue number through `motif_residues`, which is the same index space.

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
        out[state] = [found.get(residue) for residue in residues]
    return out


# Section 5.4 and 5.5's measurements, as rows: a label, where the value is in
# a frame's `measurements`, its unit, and where a state word for it lives. The
# fifth element is where the analysis's own resolved residue label is, which is
# preferred over the written label because it names the actual residues: "R3.50
# (R131) - E6.30 (E268)" for the ionic lock rather than "Ionic lock". The sixth
# suffixes the label, because one measurement can give two angles.
#
# Ported from the page this replaced, unchanged. It is a list of where things
# are in a payload, and rewriting it would be rewriting the payload's shape
# from memory.
MEASUREMENT_ROWS = [
    ("alphaC salt bridge", ("alphac_salt_bridge", "distance_a"), "\u00c5",
     ("alphac_salt_bridge", "state"), ("alphac_salt_bridge", "pair")),
    ("alphaC axis rotation", ("alphac_displacement", "axis_rotation_deg"), "\u00b0", None),
    ("alphaC centroid shift", ("alphac_displacement", "centroid_shift_a"), "\u00c5", None),
    ("P-loop tip to hinge", ("p_loop", "tip_to_hinge_a"), "\u00c5", None),
    ("R-spine", ("r_spine", "state"), "", None),
    ("TM3\u2013TM6 cleft", ("micro_switches", "tm3_tm6_a"), "\u00c5", None,
     ("micro_switches", "tm3_tm6_pair")),
    ("Ionic lock", ("micro_switches", "ionic_lock", "distance_a"), "\u00c5",
     ("micro_switches", "ionic_lock", "state"),
     ("micro_switches", "ionic_lock", "pair")),
    ("Y5.58\u2013Y7.53", ("micro_switches", "y5x58_y7x53_oh_a"), "\u00c5", None),
    ("TM6 axis rotation", ("tm6_rotation", "axis_rotation_deg"), "\u00b0", None),
    ("TM6 centroid shift", ("tm6_rotation", "centroid_shift_a"), "\u00c5", None),
    ("PIF connector RMSD", ("connector", "rmsd_a"), "\u00c5", None),
    ("Helix 8 to the membrane plane", ("helix8", "angle_to_membrane_plane_deg"),
     "\u00b0", None),
    # Section 5.5's sodium. The two pocket distances were always here; the ion
    # itself could not be reported until the analysis was given ion positions
    # to look at, and a row that is absent for a system with no sodium is
    # dropped by `frame_measurements` rather than shown empty.
    ("D2.50 to S3.39", ("micro_switches", "sodium_pocket", "d_s3x39_a"),
     "\u00c5", None),
    ("D2.50 to N7.49", ("micro_switches", "sodium_pocket", "d_n7x49_a"),
     "\u00c5", None),
    ("Nearest sodium to D2.50",
     ("micro_switches", "sodium_pocket", "nearest_sodium_a"), "\u00c5", None),
    # Section 5.5 lists the toggle switch as "W6.48 chi1 and chi2" and NPxxY
    # as "Y7.53 chi1/chi2 ...". Both were computed, both carried their own
    # resolved label, and neither had a row, so four of the section's angles
    # existed in the payload and on no page. A sixth element suffixes the
    # label, because one measurement gives two of them.
    ("Toggle switch", ("micro_switches", "toggle_w6x48", "chi1_deg"),
     "\u00b0", None, ("micro_switches", "toggle_w6x48", "position"),
     " \u03c7\u2081"),
    ("Toggle switch", ("micro_switches", "toggle_w6x48", "chi2_deg"),
     "\u00b0", None, ("micro_switches", "toggle_w6x48", "position"),
     " \u03c7\u2082"),
    ("NPxxY tyrosine", ("micro_switches", "npxxy_y7x53", "chi1_deg"),
     "\u00b0", None, ("micro_switches", "npxxy_y7x53", "position"),
     " \u03c7\u2081"),
    ("NPxxY tyrosine", ("micro_switches", "npxxy_y7x53", "chi2_deg"),
     "\u00b0", None, ("micro_switches", "npxxy_y7x53", "position"),
     " \u03c7\u2082"),
]


def _dig(payload: dict, path: tuple) -> Any:
    """Follow a key path, returning None rather than raising on a gap."""
    found: Any = payload
    for key in path:
        if not isinstance(found, dict):
            return None
        found = found.get(key)
    return found


def measurements_table(payload: dict[str, Any]) -> dict[str, Any] | None:
    """What every published criterion measured, at the start and at the end.

    These are the measurements the state calls are made from: the DFG
    distances and the R-spine for a kinase, the TM3-TM6 cleft, the ionic lock,
    the toggle switch and the NPxxY tyrosine for a receptor. The score panels
    above are built out of them, and this is where a reader checks that the
    score is saying what the geometry says.

    First and last rather than every frame: the table says what moved over the
    run, and the per-frame series belongs in a plot. A reading absent from both
    ends is left out rather than shown blank, so the table is what was measured
    rather than a list of what might have been.
    """
    frames = [frame for frame in (payload.get("frames") or [])
              if frame.get("measurements")]
    if not frames:
        return None
    first, last = frames[0]["measurements"], frames[-1]["measurements"]

    rows = []
    for entry in MEASUREMENT_ROWS:
        label, value_path, unit, state_path = entry[:4]
        label_path = entry[4] if len(entry) > 4 else None
        suffix = entry[5] if len(entry) > 5 else ""
        start, end = _dig(first, value_path), _dig(last, value_path)
        if start is None and end is None:
            continue
        # The resolved label where there is one, falling back to the written
        # one: a row whose analysis could not name its residues still has to
        # say what it measured.
        resolved = _dig(last, label_path) if label_path else None
        numeric = isinstance(end if end is not None else start, (int, float))
        state = _dig(last, state_path) if state_path else None

        def fmt(value: Any) -> str:
            if value is None:
                return ""
            if isinstance(value, (int, float)):
                return f"{value:.2f}{unit}" if unit else f"{value:.2f}"
            return str(value)

        change = ""
        if numeric and start is not None and end is not None:
            change = f"{end - start:+.2f}"
        rows.append({"cells": [
            {"text": (resolved or label) + suffix},
            {"text": fmt(start)},
            {"text": fmt(end)},
            {"text": change},
            {"text": state or "",
             "colour": state_colour(str(state).lower()) if state else None},
        ]})
    if not rows:
        return None
    return {
        "columns": ["Measurement", "At the start", "At the end", "Change",
                    "State at the end"],
        "rows": rows,
        "note": ("Every residue here was resolved at run time from KLIFS "
                 "pocket positions or GPCRdb generic numbers, which is why the "
                 "labels name actual residues. A measurement with no reading "
                 "at either end is left out rather than shown empty."),
    }


def motif_table(payload: dict[str, Any]) -> dict[str, Any] | None:
    """Every motif at the last frame, with the reference it is nearest to.

    Section 6.2's six columns. The distances are what the score is computed
    from, so they are shown beside it: a motif at s = +0.4 that is 1 angstrom
    from one ensemble and 2.3 from the other is a different claim from one at
    the same score that is 7 and 16 angstrom from them, and the score alone
    cannot tell them apart.
    """
    frames = payload.get("frames") or []
    if not frames:
        return None
    last = frames[-1]
    motifs = last.get("motifs") or []
    if not motifs:
        return None
    kind = payload.get("kind") or "kinase"
    residues = payload.get("motif_residues") or {}

    halves = _half_means(payload)
    agreement = {name: (block or {}).get("reference_agreement") or {}
                 for name, block in
                 ((payload.get("dssp") or {}).get("motifs") or {}).items()}

    rows = []
    for motif in sorted(motifs, key=lambda item: item.get("motif", "")):
        name = motif.get("motif", "")
        numbers = residues.get(name) or []
        label = motif.get("label") or "intermediate"
        first, second = halves.get(name, (None, None))
        leans = agreement.get(name) or {}
        # Withheld unless the analysis called it decisive. A lean read off two
        # residues is not a finding, and the payload says so in its own note:
        # printing "leans inactive" beside a margin that small would put a
        # conclusion on the page that the measurement does not support.
        said = (leans.get("leans") or "") if leans.get("decisive") else ""
        rows.append({
            "cells": [
                {"text": name, "colour": motif_colour(kind, name)},
                {"text": (f"{len(numbers)}, {min(numbers)} to {max(numbers)}"
                          if numbers else "")},
                {"text": (f"{motif['d_act_a']:.2f}"
                          if motif.get("d_act_a") is not None else "")},
                {"text": (f"{motif['d_inact_a']:.2f}"
                          if motif.get("d_inact_a") is not None else "")},
                {"text": (f"{motif['s']:+.3f}"
                          if motif.get("s") is not None else "")},
                {"text": (f"{first:+.3f} to {second:+.3f}"
                          if first is not None and second is not None else "")},
                {"text": label, "colour": state_colour(label)},
                {"text": said,
                 "colour": state_colour(said) if said in
                 ("active", "inactive") else None},
                {"text": " / ".join(
                    part for part in (motif.get("nearest_active"),
                                      motif.get("nearest_inactive")) if part)},
            ],
        })

    return {
        "columns": ["Motif", "Residues", "To active (Å)", "To inactive (Å)",
                    "s", "First half to second", "Label",
                    "Secondary structure leans", "Nearest active / inactive"],
        "rows": rows,
        "note": (f"At {last.get('time_ns', 0):.1f} ns, the last frame analysed. "
                 f"A label is decided beyond s = ±{DECISIVE_S}; inside that the "
                 f"motif reads as intermediate, which is a measurement and not "
                 f"a missing answer. The two half means say whether a motif was "
                 f"still moving: one number at the end cannot distinguish a "
                 f"motif that settled from one that was on its way somewhere. "
                 f"The last column is left blank where too few residues "
                 f"distinguish the reference sets to read a state from their "
                 f"secondary structure."),
    }


def _half_means(payload: dict[str, Any]) -> dict[str, tuple[float, float]]:
    """Each motif's mean score over the first and second half of the run.

    The committed figures carry these and the page did too, for a reason worth
    keeping: a single number at the end cannot tell a motif that settled from
    one that was still moving when the run stopped. Computed here from the same
    frames the timelines are drawn from, so the page and the figure cannot
    disagree about them.
    """
    frames = payload.get("frames") or []
    if len(frames) < 4:
        return {}
    cut = len(frames) // 2
    totals: dict[str, list[list[float]]] = {}
    for index, frame in enumerate(frames):
        half = 0 if index < cut else 1
        for motif in frame.get("motifs") or []:
            value = motif.get("s")
            if value is None:
                continue
            bucket = totals.setdefault(motif.get("motif", ""), [[], []])
            bucket[half].append(float(value))
    out = {}
    for name, (first, second) in totals.items():
        if first and second:
            out[name] = (sum(first) / len(first), sum(second) / len(second))
    return out


def ligand_table(payload: dict[str, Any]) -> dict[str, Any] | None:
    """The bound ligand, as the reference-state analysis measured it.

    The Analysis tab's own two ligand panels plot the pose RMSD and the contact
    occupancy, and they are on this page already. Three things it does not say
    are worth saying, and each is here because leaving it out would make a
    number read as more than it is:

    * whether the RMSD was symmetry corrected. An uncorrected pose RMSD on a
      ligand with equivalent atoms overstates drift, sometimes by angstroms.
    * what cutoff a contact was counted at, because "15 contacts" means nothing
      without it.
    * each contact's closest approach to two decimals. At one, a contact at the
      cutoff prints as the cutoff and reads as no contact at all.
    """
    ligand = payload.get("ligand") or {}
    if not ligand:
        return None

    rows = []

    def add(label: str, text: str) -> None:
        if text:
            rows.append({"cells": [{"text": label}, {"text": text}]})

    add("Residue name", str(ligand.get("resname") or ""))
    if ligand.get("n_heavy_atoms"):
        add("Heavy atoms", str(ligand["n_heavy_atoms"]))
    if ligand.get("final_pose_rmsd_a") is not None:
        add("Pose RMSD at the end", f"{ligand['final_pose_rmsd_a']:.2f} Å")
    if ligand.get("max_pose_rmsd_a") is not None:
        add("Widest pose RMSD", f"{ligand['max_pose_rmsd_a']:.2f} Å")
    add("Symmetry corrected",
        "yes" if ligand.get("symmetry_corrected") else "no")
    if ligand.get("contact_cutoff_a") is not None:
        add("Contact cutoff", f"{ligand['contact_cutoff_a']:.1f} Å")
    for contact in (ligand.get("persistent_contacts") or [])[:12]:
        fraction = contact.get("fraction_of_frames")
        closest = contact.get("closest_a")
        add(f"{contact.get('name', '')}{contact.get('residue', '')}",
            (f"{fraction:.0%} of frames" if fraction is not None else "")
            + (f", closest {closest:.2f} Å" if closest is not None else ""))
    if not rows:
        return None
    return {
        "columns": ["", ""],
        "rows": rows,
        "note": ("A contact is counted when any heavy atom of the ligand comes "
                 "within the cutoff of any heavy atom of a residue, and the "
                 "closest approach is printed to two decimals: at one, a "
                 "contact at the cutoff prints as the cutoff and reads as no "
                 "contact."),
    }


def references_table(payload: dict[str, Any]) -> dict[str, Any] | None:
    """The reference ensemble, with what the catalogue calls each structure and
    what this build computed from its coordinates.

    One table rather than two. The page this replaced had a reference browser
    on a page of its own and a separate agreement table, and the first attempt
    at moving them both made two panels that were 80 per cent the same columns.
    They are one question: here is the ensemble this run was scored against,
    and here is whether our spatial label agrees with the annotation it came
    with. That agreement, counted, is P1 and P2 for the kinase and P3 for the
    receptor, which is why the panel carries the counts in its note.

    This run's own system only. A reader on a kinase page has no use for 116
    receptor structures, and the page that showed both was reached from a tab
    that no longer exists.
    """
    references = payload.get("reference_structures") or []
    if not references:
        return None
    system = payload.get("system_id") or ""
    ensembles = payload.get("ensembles") or {}
    in_ensemble = {str(pdb).upper()
                   for state in ("active", "inactive")
                   for pdb in (ensembles.get(state) or [])}

    rows = []
    agreed = compared = 0
    for reference in sorted(references, key=lambda item: item.get("pdb", "")):
        metrics = reference.get("metrics") or {}
        pdb = reference.get("pdb", "")
        catalogue = reference.get("source_state") or ""
        ours = reference.get("our_state") or ""
        if system == "abl1":
            # KLIFS annotates the DFG and the alphaC separately, and P1 and P2
            # are those two labels. The spatial label is in the same column so
            # the two sit side by side rather than in two tables.
            extra = [{"text": metrics.get("their_dfg") or ""},
                     {"text": metrics.get("our_dfg") or ""}]
            agrees = metrics.get("dfg_agrees")
        else:
            distance = metrics.get("tm3_tm6")
            extra = [{"text": f"{distance:.2f}" if distance is not None else ""}]
            agrees = (None if not catalogue or not ours
                      else catalogue.lower() == ours.lower())
        if agrees is not None:
            compared += 1
            agreed += bool(agrees)
        verdict = ("agrees" if agrees else
                   ("differs" if agrees is False else "no opinion"))
        rows.append({"cells": [
            {"text": pdb},
            {"text": reference.get("chain", "")},
            {"text": (f"{reference['resolution']:.2f}"
                      if reference.get("resolution") is not None else "")},
            {"text": catalogue,
             "colour": state_colour(catalogue.lower()) if catalogue else None},
            {"text": ours, "colour": state_colour(ours.lower()) if ours else None},
            *extra,
            {"text": verdict,
             "colour": (state_colour("active") if agrees else
                        (state_colour("inactive") if agrees is False else None))},
            {"text": "scored" if pdb.upper() in in_ensemble else "not used"},
        ]})

    if system == "abl1":
        columns = ["PDB", "Chain", "Resolution (Å)", "KLIFS state", "Our state",
                   "KLIFS DFG", "Our DFG", "Agreement", "In the ensemble"]
        control = (f"P1 and P2: our spatial DFG and alphaC labels agree with "
                   f"KLIFS on {agreed} of {compared} structures that carry one. "
                   f"Both labels are measurements here and annotations there.")
    else:
        columns = ["PDB", "Chain", "Resolution (Å)", "GPCRdb state", "Our state",
                   "TM3-TM6 (Å)", "Agreement", "In the ensemble"]
        control = (f"P3: the threshold on the TM3-TM6 distance, fitted on this "
                   f"set under leave-one-out, reproduces GPCRdb's own state "
                   f"call on {agreed} of {compared} structures. No structure "
                   f"helped decide the line it is judged by.")
    return {
        "columns": columns,
        "rows": rows,
        "note": (f"{len(rows)} structures, resolved at run time from the "
                 f"catalogue rather than listed in the code. {control} A "
                 f"structure reads \"not used\" when it has no usable state "
                 f"label or the alignment dropped it; the core panel says how "
                 f"many that was."),
    }


def core_table(payload: dict[str, Any]) -> dict[str, Any] | None:
    """The superposition core, and how well the references agree on it.

    Every distance on this page is measured after superposing on this core, so
    how tight it is bounds what any of them can mean. The numbers are the
    alignment's own report rather than a restatement of it.
    """
    core = payload.get("core") or {}
    report = payload.get("report") or {}
    alignment = report.get("alignment") or {}
    rigidity = report.get("rigidity") or {}
    if not core and not alignment:
        return None

    rows = []

    def add(label: str, value: Any, unit: str = "") -> None:
        if value in (None, "", []):
            return
        text = (f"{value:.2f}{unit}" if isinstance(value, float)
                else f"{value}{unit}")
        rows.append({"cells": [{"text": label}, {"text": text}]})

    add("Core residues", core.get("n_residues"))
    add("Core defined by", core.get("source"))
    add("Alignment target", alignment.get("target"))
    add("References aligned", alignment.get("n_aligned"))
    add("References dropped", alignment.get("n_dropped"))
    add("Mean core RMSD", alignment.get("mean_core_rmsd_a"), " Å")
    add("Worst core RMSD", alignment.get("max_core_rmsd_a"), " Å")
    add("Rigidity reference", rigidity.get("reference"))
    add("Mean rigidity RMSD", rigidity.get("mean_rmsd_a"), " Å")
    if not rows:
        return None
    return {
        "columns": ["", ""],
        "rows": rows,
        "note": ("The core is the part of the fold that does not move between "
                 "states, so superposing on it leaves the motifs free to."),
    }


# ===========================================================================
#  The set
# ===========================================================================


def viewer_controls(results) -> dict[str, Any]:
    """What the viewer card offers beyond the scrubber, for this results file.

    `analysis.js` has carried all three of these since the page that is now
    gone: the cartoon coloured by motif, a ghost per reference state, and the
    bilayer on a toggle. They are driven entirely by markup, so moving the runs
    into the library needed the data attributes rather than new viewer code.

    Returns an empty dict for a file with no reference-state analysis, so the
    card renders exactly as it does for every other run.
    """
    member = (results.members or {}).get("refstates.json")
    if not member:
        return {}
    try:
        payload = json.loads(member)
    except (TypeError, ValueError):
        return {}
    viewer = payload.get("viewer") or {}
    members = results.members or {}

    # Offered only where the member it needs is actually in the file. A toggle
    # whose structure 404s is worse than no toggle: it reads as a broken page
    # rather than as a run that was packed without a ghost.
    overlay = {}
    for state, detail in (viewer.get("overlay") or {}).items():
        name = (detail or {}).get("member")
        if name and name in members:
            overlay[state] = detail
    lipids = viewer.get("lipids") or {}
    # The topology is required and the trajectory is optional: the bilayer is
    # the whole lipid now rather than its phosphorus, which is 2 MB for one
    # frame and 48 for five hundred, so it is a single frame and has no
    # trajectory to go with it. An older file that carries both still works.
    if lipids.get("topology") not in members:
        lipids = {}
    elif lipids.get("trajectory") and lipids["trajectory"] not in members:
        lipids = dict(lipids)
        lipids.pop("trajectory")
    return {
        "motif_painted": bool(viewer.get("motif_painted")),
        "overlay": overlay,
        "lipids": lipids,
        # Offered for any run that carries one, membrane or not. It was gated
        # on the bilayer when the only thing that could define a view was the
        # membrane normal; a water run's view comes from the shape of its own
        # core, and the point of it is that the two runs of a system open
        # facing the same way.
        "orientation": viewer.get("orientation") or {},
        # The key to the cartoon's colours, for a legend under the viewport.
        # Section 6.1 gives each motif a token and the structure carries the
        # token's number in its B-factor column, so the viewer draws colours
        # that nothing on the page named: a reader saw an orange helix and a
        # blue one with no way to tell which motif either was.
        "motifs": ([{"name": motif,
                     "colour": motif_colour(payload.get("kind") or "kinase", motif),
                     "residues": (f"{min(numbers)} to {max(numbers)}"
                                  if numbers else "")}
                    for motif, numbers in
                    sorted((payload.get("motif_residues") or {}).items())]
                   if viewer.get("motif_painted") else []),
        "kind": payload.get("kind") or "",
    }


def from_results(results) -> list[dict[str, Any]]:
    """Every reference-state panel a results file can support, in reading order.

    Returns an empty list for a results file with no `refstates.json`, which is
    every run but these four, so the run page simply has nothing extra to show.
    """
    member = (results.members or {}).get("refstates.json")
    if not member:
        return []
    try:
        payload = json.loads(member)
    except (TypeError, ValueError):
        return []
    # How many frames the viewer actually holds, which only the manifest knows.
    # The state map needs it to say which frame a click should move to.
    payload["frames_packed"] = (results.trajectory or {}).get("frames_packed")
    return build(payload)


def build(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """The panels, from the payload. Separate from `from_results` for the tests."""
    panels: list[dict[str, Any]] = []

    def chart(panel_id: str, figure: dict | None, blurb: str,
              wide: bool = True) -> None:
        if figure is None:
            return
        panels.append({"id": f"refstates_{panel_id}", "kind": "figure",
                       "title": figure["layout"]["title"]["text"],
                       "figure": figure, "blurb": blurb, "wide": wide})

    def table(panel_id: str, title: str, built: dict | None, blurb: str,
              wide: bool = True) -> None:
        if built is None:
            return
        panels.append({"id": f"refstates_{panel_id}", "kind": "table",
                       "title": title, "table": built, "blurb": blurb,
                       "wide": wide})

    window = payload.get("smoothing_window_ns")
    chart("state", state_score(payload),
          "The whole protein as one number, against the experimental active and "
          "inactive structures of this same protein: +1 is the active "
          "references, -1 the inactive ones. The thin line is per frame and the "
          "bold one is smoothed"
          + (f" over {window} ns" if window else "")
          + f". The shaded bands start at s = ±{DECISIVE_S}, which is where a "
            "label is decided; inside them the protein reads as intermediate.")
    chart("motifs", motif_timelines(payload),
          "Each structural motif scored separately, in the colours section 6.1 "
          "fixes for them. The whole-protein score above is these weighted by "
          "residue count, which is why a small motif moving a long way can "
          "leave it almost unchanged. Motifs of unequal discriminating power "
          "carry equal residue weight, so read this panel rather than the "
          "headline number.")
    table("motif_table", "Motifs at the end of the run", motif_table(payload),
          "The last frame, motif by motif, with the distances the score is "
          "computed from and the nearest reference structure in each ensemble.")
    table("measurements", "What the published criteria measured",
          measurements_table(payload),
          "The measurements every state call on this page is made from, at the "
          "first frame and the last. Each residue was resolved at run time from "
          "a KLIFS pocket position or a GPCRdb generic number, so the labels "
          "name the actual residues rather than a motif in the abstract. This "
          "is where a reader checks that the score agrees with the geometry.")
    table("ligand", "The ligand, as the reference-state analysis saw it",
          ligand_table(payload),
          "The pose and the contacts are plotted above by the Analysis tab's "
          "own ligand panels. This is what those panels do not say: whether "
          "the RMSD was symmetry corrected, what cutoff a contact was counted "
          "at, and how close each persistent contact actually came.",
          wide=False)
    chart("ribbon", dssp_ribbon(payload),
          "Secondary structure per residue per frame, with each motif marked "
          "down the left edge and the reference ensembles' own consensus in the "
          "two columns on the right. A motif whose secondary structure matches "
          "the inactive references and not the active ones is the clearest "
          "evidence this analysis produces, and it is visible along a row.")
    chart("map", state_map(payload),
          "The reference structures on the plane their published criterion is "
          "drawn on, coloured by what the catalogue calls them, with this run's "
          "path across it in dark grey. This is the panel where a threshold "
          "that does not separate the clusters is obvious rather than argued "
          "about.")
    chart("bilayer", bilayer(payload),
          "The bilayer QC, chunk by chunk, against both of section 5.8's "
          "windows: the expected one it should sit in and the wider acceptable "
          "one outside which the membrane invalidates the run. A chunk between "
          "them is reported rather than failed.")
    table("core", "Superposition core", core_table(payload),
          "What every distance on this page was measured after superposing on, "
          "and how well the reference structures agree about it.", wide=False)
    table("references", "Reference structures", references_table(payload),
          "The ensemble this run was scored against, and whether our spatial "
          "label agrees with the annotation each structure came with. Counted, "
          "that agreement is the positive control for the whole method: if the "
          "labels did not reproduce the catalogue on structures whose state is "
          "known, nothing this page says about a trajectory would mean "
          "anything.")
    return panels
