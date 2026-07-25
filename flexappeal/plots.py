"""Plotly figure builders for the Analysis tab.

Returns plain dicts (Plotly figure JSON) rather than HTML: the vendored
plotly.js draws them client-side, so the server never renders an image and the
payload stays small.

Colour
------
Every colour here is validated, not chosen by eye. The categorical set is
derived from the marcdeller.com brand hues, darkened where the brand values
failed:

    brand #fcb900 (amber)  L 0.826 -- above the 0.43-0.77 categorical band
    brand #00d084 (green)  1.74:1 -- below the 3:1 contrast floor on white
    brand #ff6900 (orange) 2.89:1 -- below the 3:1 contrast floor on white

The chart steps below pass all six checks at light/#ffffff (lightness band,
chroma floor, CVD separation, normal-vision floor, contrast). The brand values
stay in brand.css for UI chrome, where the thresholds do not apply.

Two hard constraints that shape the figures:

* **Adjacent-pair validation caps categorical use at six series; all-pairs
  caps it at three.** Any form where two arbitrary marks can touch -- a
  category heatmap, a scatter -- uses at most the first three. DSSP has exactly
  three real states; "no assignment" is neutral grey, which is what an absence
  should look like anyway.

* **No dual-axis charts, ever.** The convergence panel plots five measures on
  wildly different scales (kJ/mol, K, bar, nm3, g/mL). Two y-axes would invent
  a correlation that is not in the data, so it is small multiples instead.
"""

from __future__ import annotations

import csv
import io
import math
from typing import Any

import numpy as np

# --- validated categorical set (adjacent pairs, light mode, white surface) ---
# node scripts/validate_palette.js "#1e73be,#d95f00,#00915c,#8a3fd1,#a87c00,#c02c2e"
#   --mode light --surface "#ffffff"   ->  ALL CHECKS PASS
CATEGORICAL = ("#1e73be", "#d95f00", "#00915c", "#8a3fd1", "#a87c00", "#c02c2e")

# The first three also clear the harder all-pairs test (worst ΔE 8.1 protan).
CATEGORICAL_ALL_PAIRS = CATEGORICAL[:3]

# Single-hue sequential ramp for magnitude. The lightest step recedes toward the
# surface on purpose: on a contact map, "near zero" should read as empty.
SEQUENTIAL = (
    (0.00, "#f2f7fd"), (0.15, "#dbe9f7"), (0.30, "#b3d1ec"), (0.45, "#84b3de"),
    (0.60, "#5b9bd5"), (0.75, "#2f7cc0"), (0.90, "#1e5c92"), (1.00, "#123c61"),
)

# Ink, never series colour: text stays in these regardless of what it labels.
INK = "#1a1a2e"
INK_MUTED = "#6b7c93"
GRID = "#eef2f7"
AXIS = "#dde4ed"
NEUTRAL = "#c8d0da"          # "no data" / unassigned, not a categorical slot

FONT = "Inter, -apple-system, BlinkMacSystemFont, sans-serif"
MONO = "Roboto Mono, monospace"


def _layout(title: str, x: str, y: str, **extra: Any) -> dict[str, Any]:
    """Shared layout: recessive axes, unified hover, no chart junk."""
    layout: dict[str, Any] = {
        "title": {"text": title, "font": {"family": FONT, "size": 15, "color": INK},
                  "x": 0, "xanchor": "left", "pad": {"l": 4, "b": 8}},
        "autosize": True,
        "margin": {"l": 62, "r": 18, "t": 44, "b": 48},
        "font": {"family": FONT, "size": 12, "color": INK_MUTED},
        "paper_bgcolor": "rgba(0,0,0,0)",
        "plot_bgcolor": "rgba(0,0,0,0)",
        "xaxis": {
            "title": {"text": x, "font": {"size": 12, "color": INK_MUTED}},
            "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
            "ticks": "outside", "tickcolor": AXIS, "tickfont": {"size": 11},
        },
        "yaxis": {
            "title": {"text": y, "font": {"size": 12, "color": INK_MUTED}},
            "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
            "ticks": "outside", "tickcolor": AXIS, "tickfont": {"size": 11},
        },
        # A crosshair plus one tooltip for the whole x position: the default for
        # any time series, so a reader never has to hit a 2px line exactly.
        "hovermode": "x unified",
        "hoverlabel": {"bgcolor": "#ffffff", "bordercolor": AXIS,
                       "font": {"family": MONO, "size": 12, "color": INK}},
        "showlegend": False,
    }
    layout.update(extra)
    return layout


def _line(x, y, name: str, colour: str = CATEGORICAL[0], hover: str = "") -> dict[str, Any]:
    return {
        "type": "scatter", "mode": "lines", "name": name,
        "x": list(x), "y": list(y),
        "line": {"color": colour, "width": 2},   # thin marks
        "hovertemplate": hover or "%{y:.3f}<extra></extra>",
    }


def _figure(data: list[dict], layout: dict) -> dict[str, Any]:
    return {"data": data, "layout": layout}


def _finite(values) -> list[float]:
    """Replace NaN/inf with None so Plotly draws a gap rather than a spike."""
    return [None if (v is None or not math.isfinite(float(v))) else float(v) for v in values]


# ===========================================================================
#  Convergence -- from state_data.csv
# ===========================================================================

# Column name in the CSV -> (label, unit, whether higher-is-not-meaningful)
_STATE_COLUMNS = {
    "Potential Energy (kJ/mole)": ("Potential energy", "kJ/mol"),
    "Kinetic Energy (kJ/mole)": ("Kinetic energy", "kJ/mol"),
    "Total Energy (kJ/mole)": ("Total energy", "kJ/mol"),
    "Temperature (K)": ("Temperature", "K"),
    "Box Volume (nm^3)": ("Box volume", "nm³"),
    "Density (g/mL)": ("Density", "g/mL"),
}


def parse_state_data(raw: bytes) -> dict[str, list[float]]:
    """Read OpenMM's StateDataReporter CSV into columns."""
    text = raw.decode("utf-8", errors="replace")
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        return {}
    # OpenMM writes the header with a leading '#' and quotes each name.
    header = [h.strip().lstrip("#").strip().strip('"') for h in header]

    columns: dict[str, list[float]] = {name: [] for name in header}
    for row in reader:
        if len(row) != len(header):
            continue
        for name, value in zip(header, row):
            try:
                columns[name].append(float(value))
            except ValueError:
                columns[name].append(float("nan"))
    return columns


def convergence(columns: dict[str, list[float]]) -> list[dict[str, Any]]:
    """One figure per measure -- small multiples, never a shared axis.

    Potential energy is ~10^5 kJ/mol, temperature ~10^2 K, density ~1 g/mL.
    Putting any two of those on one plot with two y-scales would make their
    alignment arbitrary and the apparent correlation fictional.
    """
    if not columns:
        return []

    time = columns.get("Time (ps)") or columns.get("Step") or []
    x_label = "Time (ps)" if "Time (ps)" in columns else "Step"
    if not time:
        return []

    figures = []
    for index, (key, (label, unit)) in enumerate(_STATE_COLUMNS.items()):
        values = columns.get(key)
        if not values or all(not math.isfinite(v) for v in values):
            continue
        colour = CATEGORICAL[index % len(CATEGORICAL)]
        figures.append({
            "id": f"conv-{key.split(' ')[0].lower()}",
            "title": label,
            "figure": _figure(
                [_line(time, _finite(values), label, colour,
                       hover=f"%{{y:.4g}} {unit}<extra></extra>")],
                _layout(f"{label}  ({unit})", x_label, unit,
                        margin={"l": 66, "r": 14, "t": 38, "b": 44}),
            ),
        })
    return figures


# ===========================================================================
#  Geometry
# ===========================================================================


def rmsd(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("rmsd_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    # Ångström reads better than nanometres for RMSD; the axis says so.
    angstrom = [v * 10.0 for v in values]
    return _figure(
        [_line(time, _finite(angstrom), "RMSD",
               hover="%{y:.2f} Å<extra></extra>")],
        _layout("Backbone RMSD from the reference", "Time (ns)", "RMSD (Å)"),
    )


def radius_of_gyration(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("rgyr_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    angstrom = [v * 10.0 for v in values]
    return _figure(
        [_line(time, _finite(angstrom), "Rg", CATEGORICAL[1],
               hover="%{y:.2f} Å<extra></extra>")],
        _layout("Radius of gyration", "Time (ns)", "R<sub>g</sub> (Å)"),
    )


def rmsf(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Per-residue fluctuation, with the most mobile residues labelled.

    Selective direct labels rather than a value on every point: the question a
    reader has is "which bits move", and five labels answer it.
    """
    values = metrics.get("rmsf_nm")
    resids = metrics.get("rmsf_resids")
    if not values or not resids:
        return None

    angstrom = [v * 10.0 for v in values]
    names = metrics.get("rmsf_resnames") or [""] * len(resids)

    trace = {
        "type": "scatter", "mode": "lines", "name": "RMSF",
        "x": list(resids), "y": _finite(angstrom),
        "line": {"color": CATEGORICAL[0], "width": 2},
        "fill": "tozeroy", "fillcolor": "rgba(30,115,190,0.10)",
        "customdata": [[n] for n in names],
        "hovertemplate": "%{customdata[0]} %{x}<br>%{y:.2f} Å<extra></extra>",
    }

    # Label the five most mobile residues, spaced so the text does not collide.
    order = sorted(range(len(angstrom)), key=lambda i: angstrom[i], reverse=True)
    labelled: list[int] = []
    for i in order:
        if len(labelled) >= 5:
            break
        if all(abs(resids[i] - resids[j]) > max(4, len(resids) // 20) for j in labelled):
            labelled.append(i)

    annotations = [
        {
            "x": resids[i], "y": angstrom[i], "text": f"{names[i]}{resids[i]}",
            "showarrow": True, "arrowhead": 0, "arrowwidth": 1, "arrowcolor": AXIS,
            "ax": 0, "ay": -18,
            "font": {"family": MONO, "size": 10, "color": INK_MUTED},
        }
        for i in labelled
    ]

    return _figure([trace], _layout(
        "Per-residue fluctuation", "Residue", "RMSF (Å)",
        annotations=annotations, hovermode="closest",
    ))


def sasa(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("sasa_total_nm2")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    return _figure(
        [_line(time, _finite(values), "SASA", CATEGORICAL[2],
               hover="%{y:.1f} nm²<extra></extra>")],
        _layout("Solvent-accessible surface area", "Time (ns)", "SASA (nm²)"),
    )


def native_contacts(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("native_contacts_q")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    return _figure(
        [_line(time, _finite(values), "Q", CATEGORICAL[3],
               hover="Q = %{y:.3f}<extra></extra>")],
        _layout("Native contacts retained", "Time (ns)", "Q",
                yaxis={"title": {"text": "Q", "font": {"size": 12, "color": INK_MUTED}},
                       "range": [0, 1.02], "gridcolor": GRID, "zeroline": False,
                       "linecolor": AXIS, "ticks": "outside", "tickcolor": AXIS}),
    )


# ===========================================================================
#  Secondary structure
# ===========================================================================

# Exactly three real states, which is the all-pairs ceiling. "No assignment"
# takes neutral grey rather than a fourth hue -- an absence should look absent.
_DSSP_STATES = [
    (0, "Helix", CATEGORICAL_ALL_PAIRS[0]),
    (1, "Sheet", CATEGORICAL_ALL_PAIRS[1]),
    (2, "Coil", CATEGORICAL_ALL_PAIRS[2]),
    (3, "None", NEUTRAL),
]


def dssp(metrics: dict[str, Any], arrays: dict[str, np.ndarray]) -> dict[str, Any] | None:
    grid = arrays.get("dssp")
    time = metrics.get("time_ns")
    if grid is None or not time:
        return None

    resids = metrics.get("contact_residues") or list(range(grid.shape[1]))
    present = sorted({int(v) for v in np.unique(grid)})

    # A discrete colourscale: hard stops so no cell is drawn an interpolated
    # colour that means nothing.
    colourscale = []
    n = len(_DSSP_STATES)
    for value, _, colour in _DSSP_STATES:
        colourscale.append([value / n, colour])
        colourscale.append([(value + 1) / n, colour])

    label_lookup = {v: label for v, label, _ in _DSSP_STATES}
    text = [[label_lookup.get(int(v), "?") for v in row] for row in grid]

    trace = {
        "type": "heatmap",
        "z": grid.T.tolist(),                     # residues on y, time on x
        "x": list(time), "y": list(resids[: grid.shape[1]]),
        "colorscale": colourscale, "zmin": 0, "zmax": n,
        "showscale": False,
        "text": list(map(list, zip(*text))),
        "hovertemplate": "Residue %{y}<br>%{x:.2f} ns<br>%{text}<extra></extra>",
        "xgap": 0, "ygap": 0,
    }

    # Identity is never colour-alone: a legend is always present for >= 2
    # categories. Zero-width scatter traces are how a heatmap gets one.
    legend = [
        {
            "type": "scatter", "mode": "markers", "name": label,
            "x": [None], "y": [None], "showlegend": True,
            "marker": {"color": colour, "size": 10, "symbol": "square"},
            "hoverinfo": "skip",
        }
        for value, label, colour in _DSSP_STATES if value in present
    ]

    return _figure([trace, *legend], _layout(
        "Secondary structure over time", "Time (ns)", "Residue",
        showlegend=True,
        legend={"orientation": "h", "y": -0.18, "x": 0,
                "font": {"size": 11, "color": INK_MUTED}},
        hovermode="closest",
        margin={"l": 62, "r": 18, "t": 44, "b": 70},
    ))


def secondary_structure_fractions(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Helix and sheet content over time -- two series, so a legend is required."""
    helix = metrics.get("dssp_helix_fraction")
    sheet = metrics.get("dssp_sheet_fraction")
    time = metrics.get("time_ns")
    if not time or (not helix and not sheet):
        return None

    data = []
    if helix:
        data.append(_line(time, _finite(helix), "Helix", CATEGORICAL_ALL_PAIRS[0],
                          hover="Helix %{y:.1%}<extra></extra>"))
    if sheet:
        data.append(_line(time, _finite(sheet), "Sheet", CATEGORICAL_ALL_PAIRS[1],
                          hover="Sheet %{y:.1%}<extra></extra>"))

    return _figure(data, _layout(
        "Secondary structure content", "Time (ns)", "Fraction of residues",
        showlegend=len(data) > 1,
        legend={"orientation": "h", "y": 1.02, "yanchor": "bottom", "x": 1,
                "xanchor": "right", "font": {"size": 11, "color": INK_MUTED}},
        yaxis={"title": {"text": "Fraction of residues",
                         "font": {"size": 12, "color": INK_MUTED}},
               "range": [0, 1], "tickformat": ".0%", "gridcolor": GRID,
               "zeroline": False, "linecolor": AXIS, "ticks": "outside",
               "tickcolor": AXIS},
    ))


# ===========================================================================
#  Contacts and collective motion
# ===========================================================================


def contact_map(metrics: dict[str, Any], arrays: dict[str, np.ndarray]) -> dict[str, Any] | None:
    matrix = arrays.get("contact_map")
    if matrix is None or matrix.size == 0:
        return None
    resids = metrics.get("contact_residues") or list(range(matrix.shape[0]))

    return _figure([{
        "type": "heatmap",
        "z": matrix.tolist(),
        "x": list(resids[: matrix.shape[0]]), "y": list(resids[: matrix.shape[0]]),
        "colorscale": [[stop, colour] for stop, colour in SEQUENTIAL],
        "zmin": 0, "zmax": 1,
        "colorbar": {
            "title": {"text": "Occupancy", "font": {"size": 11, "color": INK_MUTED},
                      "side": "right"},
            "tickformat": ".0%", "thickness": 12, "len": 0.85,
            "outlinewidth": 0, "tickfont": {"size": 10, "color": INK_MUTED},
        },
        "hovertemplate": "%{y} – %{x}<br>%{z:.0%} of frames<extra></extra>",
    }], _layout(
        "Residue contact occupancy", "Residue", "Residue",
        hovermode="closest",
        yaxis={"title": {"text": "Residue", "font": {"size": 12, "color": INK_MUTED}},
               "scaleanchor": "x", "scaleratio": 1, "gridcolor": GRID,
               "zeroline": False, "linecolor": AXIS, "ticks": "outside",
               "tickcolor": AXIS},
    ))


def pca(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """PC1 vs PC2, coloured by simulated time.

    Time is a magnitude, so it takes the sequential ramp rather than categorical
    hues -- and it turns the scatter into a trajectory you can follow.
    """
    projection = metrics.get("pca_projection")
    time = metrics.get("time_ns")
    if not projection or not time or len(projection[0]) < 2:
        return None

    pc1 = [row[0] for row in projection]
    pc2 = [row[1] for row in projection]
    ratios = metrics.get("pca_variance_ratio") or [0, 0]

    return _figure([{
        "type": "scatter", "mode": "markers", "name": "frames",
        "x": pc1, "y": pc2,
        "marker": {
            "size": 8,                                   # >= 8px markers
            "color": list(time),
            "colorscale": [[stop, colour] for stop, colour in SEQUENTIAL],
            "showscale": True,
            "colorbar": {"title": {"text": "ns", "font": {"size": 11, "color": INK_MUTED},
                                   "side": "right"},
                         "thickness": 12, "len": 0.85, "outlinewidth": 0,
                         "tickfont": {"size": 10, "color": INK_MUTED}},
            "line": {"color": "#ffffff", "width": 1},    # 2px surface ring on overlap
        },
        "customdata": [[t] for t in time],
        "hovertemplate": "PC1 %{x:.2f}<br>PC2 %{y:.2f}<br>%{customdata[0]:.2f} ns<extra></extra>",
    }], _layout(
        "Essential dynamics",
        f"PC1 ({ratios[0]:.0%} of variance)",
        f"PC2 ({ratios[1]:.0%} of variance)" if len(ratios) > 1 else "PC2",
        hovermode="closest",
    ))


def clusters(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Cluster populations. One series, so every bar takes slot 1."""
    found = metrics.get("clusters")
    if not found:
        return None

    labels = [f"{i + 1}" for i in range(len(found))]
    fractions = [c.get("fraction", 0) for c in found]
    frames = [c.get("centre_frame", 0) for c in found]

    return _figure([{
        "type": "bar", "x": labels, "y": fractions,
        # One colour for every bar: colouring by value would re-encode what the
        # bar length already shows and spend the identity channel for nothing.
        "marker": {"color": CATEGORICAL[0], "line": {"width": 0}},
        "width": 0.62,
        "customdata": [[f, c.get("size", 0)] for f, c in zip(frames, found)],
        "hovertemplate": ("Cluster %{x}<br>%{y:.1%} of frames<br>"
                          "%{customdata[1]} members<br>centre: frame %{customdata[0]}"
                          "<extra></extra>"),
        "text": [f"{f:.0%}" for f in fractions],
        "textposition": "outside",
        "textfont": {"family": MONO, "size": 11, "color": INK_MUTED},
        "cliponaxis": False,
    }], _layout(
        "Conformational clusters", "Cluster", "Fraction of frames",
        hovermode="closest",
        yaxis={"title": {"text": "Fraction of frames",
                         "font": {"size": 12, "color": INK_MUTED}},
               "tickformat": ".0%", "range": [0, min(1.0, max(fractions) * 1.18)],
               "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
               "ticks": "outside", "tickcolor": AXIS},
    ))


# ===========================================================================
#  Membrane
# ===========================================================================


def area_per_lipid(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("area_per_lipid_nm2")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    return _figure(
        [_line(time, _finite(values), "APL", CATEGORICAL[4],
               hover="%{y:.3f} nm²<extra></extra>")],
        _layout("Area per lipid", "Time (ns)", "Area (nm²)"),
    )


def bilayer_thickness(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("bilayer_thickness_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    return _figure(
        [_line(time, _finite(values), "Thickness", CATEGORICAL[5],
               hover="%{y:.2f} nm<extra></extra>")],
        _layout("Bilayer thickness", "Time (ns)", "P–P distance (nm)"),
    )


# ===========================================================================
#  Stat tiles -- where the story is one number, it is not a chart
# ===========================================================================


def tiles(summary: dict[str, Any], metrics: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def add(label, value, unit="", note=""):
        if value is not None:
            out.append({"label": label, "value": value, "unit": unit, "note": note})

    duration = summary.get("duration_ns")
    if duration:
        add("Simulated", f"{duration:,.1f}" if duration >= 1 else f"{duration:.3f}", "ns")

    frames = summary.get("frames_analysed")
    if frames:
        note = ""
        packed = summary.get("frames_packed")
        if packed and packed < frames:
            note = f"{packed:,} packed for viewing"
        add("Frames", f"{frames:,}", "", note)

    rmsd_values = metrics.get("rmsd_nm")
    if rmsd_values:
        # The last decile, not the final frame: one frame is noise.
        tail = rmsd_values[max(0, len(rmsd_values) - max(1, len(rmsd_values) // 10)):]
        mean = sum(tail) / len(tail) * 10
        add("Final RMSD", f"{mean:.2f}", "Å", "mean of the last 10% of frames")

    rgyr = metrics.get("rgyr_nm")
    if rgyr:
        add("Radius of gyration", f"{sum(rgyr) / len(rgyr) * 10:.2f}", "Å", "trajectory mean")

    if metrics.get("hbond_count") is not None:
        add("Hydrogen bonds", f"{metrics['hbond_count']:,}", "", "occupancy > 10%")

    helix = metrics.get("dssp_helix_fraction")
    if helix:
        add("Helix content", f"{sum(helix) / len(helix):.0%}", "", "trajectory mean")

    ns_per_day = summary.get("ns_per_day")
    if ns_per_day:
        add("Throughput", f"{ns_per_day:,.0f}", "ns/day", summary.get("platform", ""))

    wall = summary.get("wall_hours")
    if wall:
        add("Wall time", f"{wall:,.1f}", "h")

    return out


# ===========================================================================
#  Assembly
# ===========================================================================


def build_all(results) -> dict[str, Any]:
    """Every panel the results file can support, in reading order."""
    metrics = results.metrics
    arrays = results.arrays

    panels: list[dict[str, Any]] = []

    def add(key, title, figure, blurb="", wide=False):
        if figure is None:
            return
        # The card's own <h3> names the panel, so an in-figure title would say
        # the same thing twice and steal vertical space from the plot. Blanked
        # here rather than in each builder so the figures stay self-describing
        # when used on their own (the convergence multiples keep theirs, because
        # in a grid the title is the only label a reader gets).
        figure["layout"]["title"] = {"text": ""}
        figure["layout"]["margin"] = dict(figure["layout"].get("margin", {}), t=12)
        panels.append({"id": key, "title": title, "figure": figure,
                       "blurb": blurb, "wide": wide})

    add("rmsd", "RMSD", rmsd(metrics),
        "How far the structure has moved from its reference. A plateau means "
        "the trajectory has settled; a steady climb means it has not.")
    add("rgyr", "Radius of gyration", radius_of_gyration(metrics),
        "Global compactness. A rise suggests unfolding or expansion.")
    add("rmsf", "Per-residue fluctuation", rmsf(metrics),
        "Which parts of the molecule move. Peaks are usually loops and termini.",
        wide=True)
    add("sasa", "Solvent exposure", sasa(metrics),
        "Total solvent-accessible surface area over time.")
    add("nativecontacts", "Native contacts", native_contacts(metrics),
        "The fraction of the starting structure's contacts still present.")
    add("ss_fractions", "Secondary structure content", secondary_structure_fractions(metrics),
        "Helix and sheet content over the trajectory.")
    add("dssp", "Secondary structure timeline", dssp(metrics, arrays),
        "Per-residue assignment at every frame. Vertical stripes are transient "
        "melting; horizontal bands are stable elements.", wide=True)
    add("contacts", "Contact map", contact_map(metrics, arrays),
        "How often each residue pair is in contact across the whole trajectory.",
        wide=True)
    add("pca", "Essential dynamics", pca(metrics),
        "The two dominant collective motions, with each frame coloured by time. "
        "Distinct clouds mean distinct conformational states.")
    add("clusters", "Conformational clusters", clusters(metrics),
        "Representative states found by RMSD clustering.")
    add("ligand_rmsd", "Ligand pose stability", ligand_rmsd(metrics),
        "How far the ligand moves within the binding site, with the protein "
        "held fixed. A flat trace means the pose held; a jump means it moved "
        "or left.")
    add("ligand_contacts", "Protein-ligand contacts", ligand_contacts(metrics),
        "Which residues the ligand touches, and for what fraction of the run.")
    add("membrane_apl", "Area per lipid", area_per_lipid(metrics),
        "The standard bilayer equilibration diagnostic: it should fall and then "
        "plateau. POPC sits near 0.63 nm² at 310 K. Two caveats on the absolute "
        "value: the protein's own cross-section is counted as membrane area, "
        "which inflates it for a large complex, and a freshly built bilayer "
        "starts loose and needs nanoseconds to condense.")
    add("membrane_scd", "Lipid chain order", lipid_order(metrics),
        "Order parameters down the acyl chain. A plateau near the headgroup "
        "falling away toward the tail is what a fluid bilayer looks like.")
    add("membrane_thickness", "Bilayer thickness", bilayer_thickness(metrics),
        "Phosphate-to-phosphate distance across the bilayer.")

    return {"panels": panels}


def ligand_rmsd(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Ligand pose stability, measured after aligning on the protein.

    The trajectory is superposed on the protein first, so this is movement
    within the binding site rather than tumbling of the whole complex -- which
    is the question anyone asks about a bound ligand.
    """
    values = metrics.get("ligand_rmsd_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    angstrom = [v * 10.0 for v in values]
    return _figure(
        [_line(time, _finite(angstrom), "Ligand RMSD", CATEGORICAL[3],
               hover="%{y:.2f} Å<extra></extra>")],
        _layout("Ligand pose stability", "Time (ns)", "Ligand RMSD (Å)"),
    )


def ligand_contacts(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Which residues the ligand actually touches, and how often.

    A horizontal bar chart: residue labels are long and there are up to twenty
    of them, so vertical bars would need rotated tick labels to fit.
    """
    contacts = metrics.get("ligand_contacts")
    if not contacts:
        return None

    top = contacts[:20][::-1]           # highest occupancy nearest the top
    labels = [c["residue"] for c in top]
    values = [c["occupancy"] for c in top]

    return _figure([{
        "type": "bar", "orientation": "h",
        "x": values, "y": labels,
        # One series, so every bar takes slot 1. Colouring by value would
        # re-encode the bar length and spend the identity channel for nothing.
        "marker": {"color": CATEGORICAL[0], "line": {"width": 0}},
        "hovertemplate": "%{y}<br>in contact for %{x:.0%} of frames<extra></extra>",
        "text": [f"{v:.0%}" for v in values],
        "textposition": "outside",
        "textfont": {"family": MONO, "size": 10, "color": INK_MUTED},
        "cliponaxis": False,
    }], _layout(
        "Protein–ligand contacts", "Occupancy", "Residue",
        hovermode="closest",
        margin={"l": 86, "r": 46, "t": 44, "b": 48},
        height=max(240, 22 * len(top) + 90),
        xaxis={"title": {"text": "Fraction of frames",
                         "font": {"size": 12, "color": INK_MUTED}},
               "tickformat": ".0%", "range": [0, min(1.05, max(values) * 1.2)],
               "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
               "ticks": "outside", "tickcolor": AXIS},
        yaxis={"title": {"text": "", "font": {"size": 12, "color": INK_MUTED}},
               "gridcolor": "rgba(0,0,0,0)", "zeroline": False, "linecolor": AXIS,
               "ticks": "outside", "tickcolor": AXIS,
               "tickfont": {"family": MONO, "size": 11}},
    ))


def lipid_order(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Deuterium order parameters down the acyl chain.

    Plotted as -S_CD, the literature convention, so an ordered chain reads high.
    The shape is diagnostic: a plateau near the headgroup falling away toward
    the tail end is what a fluid bilayer looks like, and a POPC plateau well
    below about 0.2 means the membrane has not finished equilibrating.
    """
    profile = metrics.get("lipid_order_parameters")
    if not profile:
        return None

    carbons = [p["carbon"] for p in profile]
    values = [p["scd"] for p in profile]

    return _figure([{
        "type": "scatter", "mode": "lines+markers", "name": "-S(CD)",
        "x": carbons, "y": values,
        "line": {"color": CATEGORICAL[2], "width": 2},
        "marker": {"size": 8, "color": CATEGORICAL[2],
                   "line": {"color": "#ffffff", "width": 1}},
        "hovertemplate": "C%{x}<br>-S<sub>CD</sub> = %{y:.3f}<extra></extra>",
    }], _layout(
        "Lipid chain order", "Acyl carbon", "-S<sub>CD</sub>",
        hovermode="closest",
    ))
