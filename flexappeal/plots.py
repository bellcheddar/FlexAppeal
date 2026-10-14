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
import hashlib
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


# --- decimation for the wire ------------------------------------------------
# What a panel is reduced to before it is written into the page. These are about
# the page, not about the data: the file keeps everything the run measured, and
# the full-resolution arrays are what `/reanalyse` works from.
#
# The numbers come from a browser rather than from taste. A matrix is stored at
# up to 500 on a side, which is right for a payload, but 500 by 500 is a quarter
# of a million numbers in the HTML: for a 162-residue receptor that measured
# 4.5 MB of JSON on one page, against a panel a few hundred pixels wide that
# cannot show a quarter of it.
WIRE_MATRIX_SIDE = 200
# Points per line panel. DeepFlex decimates to about this for the same reason.
WIRE_SERIES_POINTS = 600
# Columns a time-against-residue heat map keeps. Lower than the line panels
# because every column carries a cell per residue and, on the DSSP panel, a
# state name per cell for the hover: the 100 ns example came to 18.7 MB on one
# panel undecimated, which is the whole page budget spent on one chart.
WIRE_HEATMAP_COLUMNS = 400


def _wire_matrix(matrix) -> np.ndarray:
    """Block-reduce a matrix to something a page can carry."""
    side = matrix.shape[0]
    if side <= WIRE_MATRIX_SIDE:
        return matrix
    step = int(np.ceil(side / WIRE_MATRIX_SIDE))
    return matrix[::step, ::step]


def _wire_stride(length: int, limit: int = WIRE_SERIES_POINTS) -> int:
    """Frames to skip so a panel stays within `limit` points."""
    if length <= limit:
        return 1
    return int(np.ceil(length / limit))


def _thin(time, *series, limit: int = WIRE_SERIES_POINTS):
    """One stride applied to a time axis and everything plotted against it.

    Returned together on purpose. Striding the axis and the series separately
    is how a chart ends up drawn against the wrong x, and the two calls are far
    enough apart in a builder for that to look right while reading.
    """
    stride = _wire_stride(len(time), limit)
    if stride == 1:
        return (list(time), *[list(s) for s in series])
    return (list(time)[::stride], *[list(s)[::stride] for s in series])


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


# Answers already computed, keyed by the series they were computed from. The
# detection is superlinear in the point count and runs on every view of a run
# page: the committed example writes a state row every 10 ps for 100 ns, and
# 10,000 points cost 1.2 seconds here and about 3 on the droplet, which took
# that page from 0.44 seconds to 3.
#
# Thinning the series was the first answer and it is the wrong one: at one
# point in five the example's own answer moved from 0.52 ns to nothing at all,
# so the cheaper number is a different number. The work is done once per file
# instead, and `deploy.sh` warms every run page, so a reader never pays it.
_EQUIL_CACHE: dict[str, dict[str, Any]] = {}


def equilibration(columns: dict[str, list[float]]) -> dict[str, Any]:
    """Where the production run became stationary, and what said so.

    pymbar's `detect_equilibration` on the potential energy, which is the
    series this is conventionally read off: it is the one that falls as the
    system relaxes and then fluctuates about a mean, and it is recorded for
    every run whether or not there is a ligand in it.

    This is not the restrained equilibration that ran before production. That
    phase is reported separately, from the manifest, because it is a protocol
    setting rather than a measurement. This is the measurement: how much of the
    production trajectory is still relaxing.

    Returns an empty dict where there is no series to read or no answer to
    give, which is the rule for this file: a key that cannot be produced is
    omitted rather than filled with a zero that reads as "none needed".
    """
    series = columns.get("Potential Energy (kJ/mole)") or []
    time = columns.get("Time (ps)") or []
    finite = [v for v in series if math.isfinite(v)]
    if len(finite) < 8 or len(time) != len(series):
        return {}

    from . import measure  # noqa: PLC0415

    values = np.asarray(finite, dtype=float)
    key = hashlib.sha256(values.tobytes()).hexdigest()
    if key in _EQUIL_CACHE:
        return _EQUIL_CACHE[key]

    start, method = measure.detect_equilibration(values)
    if start <= 0:
        # Stationary from the first frame. That is an answer and it is worth
        # saying, so it is reported rather than omitted.
        answer = {"start_index": 0, "ns": 0.0, "fraction": 0.0,
                  "method": method}
    else:
        index = min(int(start), len(time) - 1)
        span = float(time[-1]) - float(time[0])
        answer = {
            "start_index": index,
            "ns": (float(time[index]) - float(time[0])) / 1000.0,
            "fraction": ((float(time[index]) - float(time[0])) / span
                         if span > 0 else 0.0),
            "method": method,
        }
    _EQUIL_CACHE[key] = answer
    return answer


def convergence(columns: dict[str, list[float]]) -> list[dict[str, Any]]:
    """One figure per measure -- small multiples, never a shared axis.

    Potential energy is ~10^5 kJ/mol, temperature ~10^2 K, density ~1 g/mL.
    Putting any two of those on one plot with two y-scales would make their
    alignment arbitrary and the apparent correlation fictional.

    Each carries the same shaded head: the part of production that `pymbar`
    says is still relaxing, read off the potential energy and drawn on all six
    so that a reader comparing them is comparing the same window.
    """
    if not columns:
        return []

    time = columns.get("Time (ps)") or columns.get("Step") or []
    x_label = "Time (ps)" if "Time (ps)" in columns else "Step"
    if not time:
        return []

    settled = equilibration(columns)
    band = _equil_band(settled.get("start_index", 0), time) if settled else []

    figures = []
    for index, (key, (label, unit)) in enumerate(_STATE_COLUMNS.items()):
        values = columns.get(key)
        if not values or all(not math.isfinite(v) for v in values):
            continue
        colour = CATEGORICAL[index % len(CATEGORICAL)]
        # Six of these on one page, each with a point per reporter interval:
        # undecimated, a 100 ns run spent 2.2 MB of the page on them between
        # them. The state file itself is packed whole and keeps every row.
        thinned_time, thinned = _thin(time, _finite(values))
        figures.append({
            "id": f"conv-{key.split(' ')[0].lower()}",
            "title": label,
            "figure": _figure(
                [_line(thinned_time, thinned, label, colour,
                       hover=f"%{{y:.4g}} {unit}<extra></extra>")],
                _layout(f"{label}  ({unit})", x_label, unit,
                        margin={"l": 66, "r": 14, "t": 38, "b": 44},
                        shapes=band),
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
    time, angstrom = _thin(time, _finite(angstrom))
    return _figure(
        [_line(time, angstrom, "RMSD",
               hover="%{y:.2f} Å<extra></extra>")],
        _layout("Backbone RMSD from the reference", "Time (ns)", "RMSD (Å)"),
    )


def radius_of_gyration(metrics: dict[str, Any]) -> dict[str, Any] | None:
    values = metrics.get("rgyr_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    angstrom = [v * 10.0 for v in values]
    time, angstrom = _thin(time, _finite(angstrom))
    return _figure(
        [_line(time, angstrom, "Rg", CATEGORICAL[1],
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
    time, values = _thin(time, _finite(values))
    return _figure(
        [_line(time, values, "SASA", CATEGORICAL[2],
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

    # Decimate along time before anything is built from the grid. Every column
    # costs a cell and a state name per residue, so a long run spends the whole
    # page budget here: 10,000 frames of a 162-residue receptor came to 18.7 MB
    # of JSON for a panel a few hundred pixels wide. The residue axis is left
    # alone -- it is the axis a reader reads, and there are never many of them.
    stride = _wire_stride(len(time), WIRE_HEATMAP_COLUMNS)
    if stride > 1:
        grid = grid[::stride]
        time = list(time)[::stride]

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

    # One stride for both series and the axis, so the two curves stay comparable
    # point for point as well as sharing an x.
    time, helix, sheet = _thin(time, _finite(helix or []), _finite(sheet or []))

    data = []
    if helix:
        data.append(_line(time, helix, "Helix", CATEGORICAL_ALL_PAIRS[0],
                          hover="Helix %{y:.1%}<extra></extra>"))
    if sheet:
        data.append(_line(time, sheet, "Sheet", CATEGORICAL_ALL_PAIRS[1],
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

    # Eight-pixel markers on a panel a few hundred pixels wide: past about this
    # many the scatter is solid colour, so the extra points cost page weight and
    # show nothing. The variance ratios still come from every frame.
    time, pc1, pc2 = _thin(time, pc1, pc2)

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


def tiles(summary: dict[str, Any], metrics: dict[str, Any],
          settled: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def add(label, value, unit="", note=""):
        if value is not None:
            out.append({"label": label, "value": value, "unit": unit, "note": note})

    duration = summary.get("duration_ns")
    if duration:
        # Whole nanoseconds for a run of a nanosecond or more. A 50 ns run
        # whose last frame lands at 49.9 is a 50 ns run, and a tile that reads
        # "49.9" invites a reader to wonder what happened to the other 0.1.
        # Under a nanosecond the decimals are the measurement, so they stay.
        add("Simulated", f"{duration:,.0f}" if duration >= 1 else f"{duration:.3f}",
            "ns")

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

    # Two different facts, and only one of them is a measurement.
    #
    # `equilibration_ns` is the restrained phase the protocol ran before
    # production, which is a setting; the tile beside it is where pymbar says
    # the production trajectory itself stopped relaxing, read off the potential
    # energy. A run can have a long restrained phase and still take a while to
    # settle, which is exactly what the apo receptor does, so reporting one as
    # though it were the other would hide the thing worth knowing.
    equilibrated = summary.get("equilibration_ns")
    if equilibrated:
        add("Equilibrated", f"{equilibrated:,.0f}" if equilibrated >= 1
            else f"{equilibrated:.2f}", "ns", "restrained, before production")

    if settled:
        by = settled["method"].split("->")[0]
        # Zero is the measurement here, not a missing answer: the series was
        # stationary from the first production frame. It says so in words,
        # because a tile reading "0.0 ns" beside a note reading "0%" invites a
        # reader to think the test did not run.
        add("Settles after", f"{settled['ns']:,.1f}", "ns",
            f"stationary from the first frame, by {by}" if not settled["ns"]
            else f"{settled['fraction']:.0%} of production, by {by}")

    # What a run was, rather than what it measured. These were a subtitle under
    # the page's title, in a line that also repeated the duration the first
    # tile carries: a reader meeting a run page wants to know which force field
    # and which structure as much as they want the RMSD, and a subtitle is not
    # where the page puts a fact worth reading.
    #
    # The platform is already the note under Throughput, where it says what
    # produced that number, and it earns a box of its own as well: what the run
    # was computed on is a property of the run.
    add("Force field", summary.get("force_field") or None)
    add("Platform", summary.get("platform") or None)
    add("Structure", (summary.get("structure") or {}).get("filename") or None)

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


def distance(metrics: dict[str, Any]) -> dict[str, Any] | None:
    """Centre-of-geometry separation between two selections over time."""
    values = metrics.get("distance_nm")
    time = metrics.get("time_ns")
    if not values or not time:
        return None
    return _figure(
        [_line(time, _finite([v * 10 for v in values]), "Distance", CATEGORICAL[5],
               hover="%{y:.2f} Å<extra></extra>")],
        _layout("Distance between selections", "Time (ns)", "Distance (Å)"),
    )


def build_reanalysis(metrics: dict[str, Any]) -> list[dict[str, Any]]:
    """Panels for a server-side re-analysis result.

    Reuses the same builders as the main page, so a re-analysed RMSD looks
    identical to one computed in the bundle -- which is the point: the user
    should not have to work out which panel came from where.
    """
    panels: list[dict[str, Any]] = []

    def add(key, title, figure):
        if figure is not None:
            figure["layout"]["title"] = {"text": ""}
            figure["layout"]["margin"] = dict(figure["layout"].get("margin", {}), t=12)
            panels.append({"id": f"re-{key}", "title": title, "figure": figure})

    add("rmsd", "RMSD", rmsd(metrics))
    add("rgyr", "Radius of gyration", radius_of_gyration(metrics))
    add("rmsf", "Per-residue fluctuation", rmsf(metrics))
    add("sasa", "Solvent exposure", sasa(metrics))
    add("ss", "Secondary structure content", secondary_structure_fractions(metrics))
    add("pca", "Essential dynamics", pca(metrics))
    add("distance", "Distance", distance(metrics))

    contact_map_values = metrics.get("contact_map")
    if contact_map_values:
        arrays = {"contact_map": np.asarray(contact_map_values, dtype=np.float32)}
        add("contacts", "Contact map", contact_map(metrics, arrays))

    return panels


# ===========================================================================
#  Binding-site panels, from flexappeal.measure
# ===========================================================================
#
# The same measurements an `.fxa` version 2 package carries, drawn for the one
# complex an upload contains. The grade that goes above them is built separately
# in `flexappeal.grade`, which is a transcription of DeepFlex's own weights and
# bands; these builders stay about what was measured.
#
# What a single complex still cannot have is anything needing the cohort: the
# energy sub-score is a percentile among the other ligands run against the same
# receptor, and the league table and both rank columns are positions in it. A
# package is what produces those.

# Interaction types, in the order they are stacked and coloured. Polar and ionic
# are the directional, specific ones and get the first two slots; hydrophobic is
# the bulk term and sits third.
_ITYPE_COLOUR = {
    "HBond": CATEGORICAL[0],
    "Ionic": CATEGORICAL[1],
    "Hydrophobic": CATEGORICAL[2],
    "contact": NEUTRAL,
}


def _equil_band(start_frame: int, time_ns) -> list[dict[str, Any]]:
    """A shaded region over the frames that were measured and not scored.

    Every scored number on these panels is cut from `equil_start_frame`, so the
    discarded head has to be visible: a reader looking at a chart that includes
    it would otherwise reasonably assume the number beside it does too.
    """
    if start_frame <= 0 or len(time_ns) == 0:
        return []
    end = float(time_ns[min(int(start_frame), len(time_ns) - 1)])
    return [{
        "type": "rect", "xref": "x", "yref": "paper",
        "x0": float(time_ns[0]), "x1": end, "y0": 0, "y1": 1,
        "fillcolor": GRID, "opacity": 0.85, "line": {"width": 0}, "layer": "below",
    }]


def pose_stability(measured) -> dict[str, Any] | None:
    """Ligand RMSD two ways: fitted on the pocket, and fitted on itself.

    The pair is the point. Fitted on the pocket, the curve is the ligand moving
    inside a receptor held still. Fitted on the ligand, it is the ligand
    changing shape. A rise in the first with a flat second is a pose sliding; a
    rise in both is a ligand doing something more interesting.
    """
    series = measured.series or {}
    if "pose_rmsd_a" not in series:
        return None
    stride = _wire_stride(len(series["time_ns"]))
    time_ns = series["time_ns"][::stride]
    layout = _layout("Ligand pose", "Time (ns)", "RMSD (Å)",
                     showlegend=True,
                     legend={"orientation": "h", "y": 1.12, "x": 0})
    layout["shapes"] = _equil_band(measured.metrics.get("equil_start_frame", 0), time_ns)
    return _figure([
        _line(time_ns, _finite(series["pose_rmsd_a"][::stride]), "pocket fit",
              CATEGORICAL[0], "pocket fit %{y:.2f} Å<extra></extra>"),
        _line(time_ns, _finite(series["inplace_rmsd_a"][::stride]), "ligand fit",
              CATEGORICAL[1], "ligand fit %{y:.2f} Å<extra></extra>"),
    ], layout)


def centroid_drift(measured) -> dict[str, Any] | None:
    """How far the ligand centroid has gone, and how far it sits from the pocket.

    The drift is measured from the ligand's own starting centroid. The offset is
    to the pocket centroid, which sits several angstrom away even in a perfectly
    held pose, so it is a diagnostic rather than a second drift measure.
    """
    series = measured.series or {}
    if "com_drift_a" not in series:
        return None
    stride = _wire_stride(len(series["time_ns"]))
    time_ns = series["time_ns"][::stride]
    layout = _layout("Centroid drift", "Time (ns)", "Distance (Å)",
                     showlegend=True,
                     legend={"orientation": "h", "y": 1.12, "x": 0})
    layout["shapes"] = _equil_band(
        measured.metrics.get("equil_start_frame", 0), time_ns)
    # The 5 A line is DeepFlex's hard gate: a mean drift above it, or a drift
    # above it in more than a fifth of frames, caps the grade at C.
    layout["shapes"].append({
        "type": "line", "xref": "x", "yref": "y",
        "x0": float(time_ns[0]), "x1": float(time_ns[-1]), "y0": 5.0, "y1": 5.0,
        "line": {"color": CATEGORICAL[5], "width": 1, "dash": "dot"},
    })
    return _figure([
        _line(time_ns, _finite(series["com_drift_a"][::stride]),
              "from its own start", CATEGORICAL[0],
              "drift %{y:.2f} Å<extra></extra>"),
        _line(time_ns, _finite(series["pocket_offset_a"][::stride]),
              "from the pocket centre", CATEGORICAL[3],
              "offset %{y:.2f} Å<extra></extra>"),
    ], layout)


def pocket_occupancy(measured) -> dict[str, Any] | None:
    """Native contacts retained per frame, against the quarter-retained threshold.

    The threshold is what `in_pocket_fraction` counts: the fraction of scored
    frames holding at least a quarter of the contacts present at the start.
    """
    series = measured.series or {}
    if "contacts" not in series:
        return None
    stride = _wire_stride(len(series["time_ns"]))
    time_ns = series["time_ns"][::stride]
    threshold = float(series["native_contact_threshold"][0]) \
        if "native_contact_threshold" in series else None
    layout = _layout("Pocket occupancy", "Time (ns)", "Residues in contact")
    layout["shapes"] = _equil_band(
        measured.metrics.get("equil_start_frame", 0), time_ns)
    if threshold is not None:
        layout["shapes"].append({
            "type": "line", "xref": "x", "yref": "y",
            "x0": float(time_ns[0]), "x1": float(time_ns[-1]),
            "y0": threshold, "y1": threshold,
            "line": {"color": CATEGORICAL[5], "width": 1, "dash": "dot"},
        })
    return _figure([
        _line(time_ns, [int(v) for v in series["contacts"][::stride]], "contacts",
              CATEGORICAL[2], "%{y} residues<extra></extra>"),
    ], layout)


def interaction_occupancy(measured) -> dict[str, Any] | None:
    """Which contacts held, and for how much of the scored window.

    Sorted by occupancy rather than by residue number: the question a reader
    brings to this panel is which interactions carried the pose, and a plot
    ordered by sequence position buries the answer.
    """
    rows = measured.summary or []
    if not rows:
        return None
    rows = sorted(rows, key=lambda r: r["occupancy"])[-24:]
    return _figure([{
        "type": "bar", "orientation": "h",
        "x": [float(r["occupancy"]) for r in rows],
        "y": [f"{r['residue']} {r['itype']}" for r in rows],
        "marker": {"color": [_ITYPE_COLOUR.get(r["itype"], NEUTRAL) for r in rows]},
        "hovertemplate": "%{y}<br>%{x:.0%} of scored frames<extra></extra>",
    }], _layout("Interaction occupancy", "Fraction of scored frames", "",
                hovermode="closest",
                xaxis={"range": [0, 1], "tickformat": ".0%",
                       "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
                       "title": {"text": "Fraction of scored frames",
                                 "font": {"size": 12, "color": INK_MUTED}}},
                margin={"l": 150, "r": 18, "t": 44, "b": 48}))


def contact_barcode(measured) -> dict[str, Any] | None:
    """Which contact was present when, busiest first.

    A barcode rather than a set of lines: twenty-plus on-off series stacked as
    lines is unreadable, and the question is when each was present rather than
    its exact value.
    """
    barcode = measured.barcode or {}
    if "z" not in barcode or not len(barcode.get("labels", [])):
        return None
    labels = [str(v) for v in barcode["labels"]]
    return _figure([{
        "type": "heatmap",
        "z": [[int(v) for v in row] for row in barcode["z"]],
        "x": [float(v) for v in barcode["time_ns"]],
        "y": labels,
        "colorscale": [(0.0, "#f2f7fd"), (1.0, CATEGORICAL[0])],
        "showscale": False,
        "hovertemplate": "%{y}<br>%{x:.1f} ns<extra></extra>",
    }], _layout("Contact barcode", "Time (ns)", "", hovermode="closest",
                margin={"l": 150, "r": 18, "t": 44, "b": 48}))


def pose_clusters(measured) -> dict[str, Any] | None:
    """How the scored window divides into distinct poses.

    Each bar is a cluster's share of the window; the x position is how far that
    cluster's centroid sits from the starting pose. One tall bar near zero is a
    ligand that held; two bars wide apart is one that moved and stayed.
    """
    clusters = measured.metrics.get("pose_clusters") or []
    if not clusters:
        return None
    return _figure([{
        "type": "bar",
        "x": [float(c["rmsd_to_t0_a"]) for c in clusters],
        "y": [float(c["population"]) for c in clusters],
        "marker": {"color": CATEGORICAL[0]},
        "width": 0.18,
        "customdata": [[int(c["cluster_id"]), int(c["n_frames"])] for c in clusters],
        "hovertemplate": ("cluster %{customdata[0]}<br>%{y:.0%} of the window"
                          "<br>%{customdata[1]} frames<br>%{x:.2f} Å from the "
                          "start<extra></extra>"),
    }], _layout("Pose clusters", "Centroid RMSD from the starting pose (Å)",
                "Share of the scored window", hovermode="closest",
                yaxis={"range": [0, 1], "tickformat": ".0%",
                       "gridcolor": GRID, "zeroline": False, "linecolor": AXIS,
                       "title": {"text": "Share of the scored window",
                                 "font": {"size": 12, "color": INK_MUTED}}}))


def ligand_strain(measured) -> dict[str, Any] | None:
    """MMFF94 strain in excess of the starting pose.

    The excess, never the absolute energy: these coordinates come from an OpenFF
    or GAFF force field and are being scored by MMFF94, and both offsets are
    constant for a given ligand and cancel when the first frame is subtracted.
    Below zero means the ligand relaxed out of the pose it started in, which is
    a good outcome rather than a neutral one.
    """
    conformer = measured.conformer or {}
    if "strain_excess_kcal" not in conformer:
        return None
    layout = _layout("Ligand strain", "Time (ns)", "Excess strain (kcal/mol)")
    layout["shapes"] = [{
        "type": "line", "xref": "paper", "yref": "y",
        "x0": 0, "x1": 1, "y0": 0, "y1": 0,
        "line": {"color": AXIS, "width": 1},
    }]
    return _figure([
        _line(_finite(conformer["time_ns"]),
              _finite(conformer["strain_excess_kcal"]),
              "excess strain", CATEGORICAL[3],
              "%{y:+.1f} kcal/mol<extra></extra>"),
    ], layout)


def ligand_torsions(measured) -> dict[str, Any] | None:
    """Each rotatable torsion over the scored window.

    Flat traces are a rigid ligand. A trace jumping between bands is a rotamer
    flip, which is what the flip rate counts.
    """
    conformer = measured.conformer or {}
    angles = conformer.get("torsion_deg")
    if angles is None or not angles.size:
        return None
    time_ns = _finite(conformer["time_ns"])
    data = [
        _line(time_ns, _finite(angles[:, column]),
              f"torsion {column + 1}", CATEGORICAL[column % len(CATEGORICAL)],
              "%{y:.0f}°<extra></extra>")
        for column in range(angles.shape[1])
    ]
    return _figure(data, _layout(
        "Ligand torsions", "Time (ns)", "Dihedral (°)",
        showlegend=angles.shape[1] > 1,
        legend={"orientation": "h", "y": 1.12, "x": 0},
        yaxis={"range": [-180, 180], "dtick": 90, "gridcolor": GRID,
               "zeroline": False, "linecolor": AXIS,
               "title": {"text": "Dihedral (°)",
                         "font": {"size": 12, "color": INK_MUTED}}}))


def ligand_atom_fluctuation(measured) -> dict[str, Any] | None:
    """Per-atom ligand movement about its mean position, after the pocket fit."""
    series = measured.series or {}
    values = series.get("ligand_rmsf_a")
    if values is None or not len(values):
        return None
    return _figure([{
        "type": "bar",
        "x": list(range(1, len(values) + 1)),
        "y": _finite(values),
        "marker": {"color": CATEGORICAL[0]},
        "hovertemplate": "atom %{x}: %{y:.2f} Å<extra></extra>",
    }], _layout("Ligand atom fluctuation", "Heavy atom", "RMSF (Å)",
                hovermode="closest"))


def pocket_flexibility(measured) -> dict[str, Any] | None:
    """Per-residue receptor fluctuation, with the pocket residues marked.

    The pocket is picked out rather than left to a reader to look up: the
    question this panel answers is whether the site itself is rigid, and that is
    invisible on an unmarked trace of two hundred residues.
    """
    panels = measured.panels or {}
    if "rmsf_resid" not in panels:
        return None
    resid = [int(v) for v in panels["rmsf_resid"]]
    values = _finite(panels["rmsf_a"])
    pocket = set(measured.pocket.get("pocket_residues_filtered") or [])
    marked = [(r, v) for r, v in zip(resid, values) if r in pocket]
    data = [{
        "type": "scatter", "mode": "lines", "name": "all residues",
        "x": resid, "y": values,
        "line": {"color": CATEGORICAL[0], "width": 1.5},
        "hovertemplate": "residue %{x}: %{y:.2f} Å<extra></extra>",
    }]
    if marked:
        data.append({
            "type": "scatter", "mode": "markers", "name": "pocket",
            "x": [r for r, _ in marked], "y": [v for _, v in marked],
            "marker": {"color": CATEGORICAL[1], "size": 7},
            "hovertemplate": "pocket residue %{x}: %{y:.2f} Å<extra></extra>",
        })
    return _figure(data, _layout(
        "Receptor flexibility", "Residue", "Cα RMSF (Å)", showlegend=True,
        legend={"orientation": "h", "y": 1.12, "x": 0}, hovermode="closest"))


# The MM/GBSA terms worth drawing, in reading order: the four physical
# contributions, then the two sums.
#
# The bonded terms are deliberately absent. BOND, ANGLE, DIHED and both 1-4
# terms cancel to zero in a single-trajectory calculation, because the receptor
# and the ligand are read out of the complex's own coordinates, so five bars at
# zero would say nothing about the binding. That they cancel is a check on the
# calculation rather than a finding, and `energy._check_cancellation` refuses the
# result when they do not.
#
# DELTA TOTAL is left out too, which is what DeepFlex's own panel does: it is the
# sum of everything else and is printed beside the chart instead of competing
# with the parts for width.
ENERGY_TERM_ORDER = ("VDWAALS", "EEL", "EGB", "ESURF",
                     "DELTA G gas", "DELTA G solv")


def energy_terms(measured) -> dict[str, Any] | None:
    """The MM/GBSA difference block, term by term.

    The only part of an MM/GBSA calculation that is a binding quantity: the
    complex, receptor and ligand blocks are each thousands of kcal/mol and only
    their difference means anything.

    Error bars are the standard deviation over the frames MMPBSA.py sampled, not
    the standard error. The spread is the interesting number here -- a van der
    Waals term of -11 with a deviation of 8 is a different claim from the same
    mean with a deviation of 1 -- and a standard error would shrink with frame
    count and hide it.
    """
    terms = (measured.metrics or {}).get("mmpbsa_gb_terms")
    if not isinstance(terms, dict) or not terms:
        return None

    rows = [(name, terms[name]) for name in ENERGY_TERM_ORDER if name in terms]
    if not rows:
        return None

    labels = [name for name, _ in rows]
    means = [float(values[0]) for _, values in rows]
    spread = [float(values[1]) if len(values) > 1 else 0.0 for _, values in rows]

    return _figure([{
        "type": "bar", "orientation": "h",
        "y": labels, "x": means,
        "error_x": {"type": "data", "array": spread, "color": INK_MUTED,
                    "thickness": 1, "width": 3},
        # Favourable and unfavourable read differently at a glance, and the sign
        # is the first thing anybody asks of an energy term.
        "marker": {"color": [CATEGORICAL[0] if value < 0 else CATEGORICAL[5]
                             for value in means]},
        "hovertemplate": "%{y}: %{x:.2f} kcal/mol<extra></extra>",
    }], _layout(
        "MM/GBSA, term by term", "\u0394G contribution (kcal/mol)", "",
        bargap=0.35, hovermode="closest",
        margin={"l": 96, "r": 18, "t": 12, "b": 44},
        xaxis={"title": {"text": "\u0394G contribution (kcal/mol)",
                         "font": {"size": 12, "color": INK_MUTED}},
               "zeroline": True, "zerolinecolor": AXIS, "gridcolor": GRID,
               "linecolor": AXIS, "ticks": "outside", "tickcolor": AXIS},
        yaxis={"automargin": True, "title": {"text": ""},
               "linecolor": AXIS, "ticks": "outside", "tickcolor": AXIS},
    ))


def rmsd_matrix(measured) -> dict[str, Any] | None:
    """Frame against frame Cα RMSD: a map of which parts of the run resemble which.

    Block structure along the diagonal is a trajectory that visited distinct
    states; a uniform field is one that stayed in one.
    """
    panels = measured.panels or {}
    matrix = panels.get("rmsd_2d")
    if matrix is None or matrix.ndim != 2 or matrix.shape[0] < 2:
        return None
    time_ns = panels["time_ns"]
    step = max(1, len(time_ns) // matrix.shape[0])
    axis = [float(v) for v in time_ns[::step][: matrix.shape[0]]]
    matrix = _wire_matrix(matrix)
    axis = axis[::max(1, int(np.ceil(len(axis) / matrix.shape[0])))][: matrix.shape[0]]
    return _figure([{
        "type": "heatmap",
        "z": [[round(float(v), 2) for v in row] for row in matrix],
        "x": axis, "y": axis,
        "colorscale": [list(stop) for stop in SEQUENTIAL],
        "colorbar": {"title": {"text": "Å", "font": {"size": 11}},
                     "thickness": 12, "outlinewidth": 0},
        "hovertemplate": "%{x:.1f} ns vs %{y:.1f} ns: %{z:.2f} Å<extra></extra>",
    }], _layout("Frame-to-frame RMSD", "Time (ns)", "Time (ns)",
                hovermode="closest"))


def cross_correlation(measured) -> dict[str, Any] | None:
    """Residue motion correlation: which parts of the receptor move together.

    Diverging, because the sign is the whole point: +1 is moving together, -1 is
    moving oppositely, and a sequential ramp would hide the difference.
    """
    panels = measured.panels or {}
    matrix = panels.get("cross_correlation")
    if matrix is None or matrix.ndim != 2 or matrix.shape[0] < 2:
        return None
    resid = [int(v) for v in panels["dccm_resid"]]
    reduced = _wire_matrix(matrix)
    if reduced.shape[0] != matrix.shape[0]:
        step = int(np.ceil(len(resid) / reduced.shape[0]))
        resid = resid[::step][: reduced.shape[0]]
    matrix = reduced
    return _figure([{
        "type": "heatmap",
        "z": [[round(float(v), 2) for v in row] for row in matrix],
        "x": resid, "y": resid, "zmin": -1, "zmax": 1,
        "colorscale": [(0.0, "#c02c2e"), (0.5, "#f7f7f7"), (1.0, "#1e5c92")],
        "colorbar": {"thickness": 12, "outlinewidth": 0,
                     "tickvals": [-1, 0, 1]},
        "hovertemplate": "%{x} vs %{y}: %{z:.2f}<extra></extra>",
    }], _layout("Residue cross-correlation", "Residue", "Residue",
                hovermode="closest"))


def _plain(value: Any) -> Any:
    """numpy scalars and arrays reduced to what `json.dumps` accepts.

    Every series on these panels came out of numpy, and a stray `float32`
    reaches the template as an unserialisable object rather than as a number:
    the panel then fails to render with a traceback about the whole page. One
    sweep where the panels are assembled beats remembering it in twelve
    builders.
    """
    if isinstance(value, np.generic):
        value = value.item()
        return round(value, 3) if isinstance(value, float) else value
    if isinstance(value, float):
        return round(value, 3)
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {key: _plain(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def build_binding(measured) -> list[dict[str, Any]]:
    """Every binding-site panel this measurement can support, in reading order."""
    panels: list[dict[str, Any]] = []

    def add(key, title, figure, blurb="", wide=False):
        if figure is None:
            return
        figure = _plain(figure)
        figure["layout"]["title"] = {"text": ""}
        figure["layout"]["margin"] = dict(figure["layout"].get("margin", {}), t=12)
        # Namespaced, because these panels share a page with `build_all`'s and
        # the id becomes a DOM id. `clusters` collided: the whole-system panels
        # already use it for conformational clustering, so a ligand run rendered
        # two elements with id="plot-clusters". Both still drew, because the
        # bootstrap walks nodes rather than ids, but the document was invalid and
        # anything reaching for one by id found whichever came first. Prefixing
        # here makes the clash impossible rather than fixing the one instance of
        # it that happened to exist.
        panels.append({"id": f"binding-{key}", "title": title, "figure": figure,
                       "blurb": blurb, "wide": wide})

    add("pose", "Ligand pose", pose_stability(measured),
        "Two fits of the same ligand. Against the pocket, the curve is the "
        "ligand moving in a receptor held still; against itself, it is the "
        "ligand changing shape. The shaded head is the equilibration that was "
        "measured and not scored.")
    add("centroid", "Centroid drift", centroid_drift(measured),
        "Displacement of the ligand centroid from where it started. The dotted "
        "line at 5 Å is the drift beyond which DeepFlex caps a grade.")
    add("occupancy", "Pocket occupancy", pocket_occupancy(measured),
        "Receptor residues in contact with the ligand. The dotted line is a "
        "quarter of the contacts present at the start, which is the threshold "
        "for counting a frame as still bound.")
    add("interactions", "Interaction occupancy", interaction_occupancy(measured),
        "Which contacts carried the pose, and for how much of the scored "
        "window. Typed from heavy-atom geometry, so these are polar, ionic and "
        "hydrophobic contacts rather than a full interaction fingerprint.",
        wide=True)
    add("barcode", "Contact barcode", contact_barcode(measured),
        "When each contact was present, busiest first. Gaps are the contact "
        "breaking and reforming.", wide=True)
    add("clusters", "Pose clusters", pose_clusters(measured),
        "How the scored window divides into distinct poses. One tall bar near "
        "zero is a ligand that held its pose.")
    add("strain", "Ligand strain", ligand_strain(measured),
        "MMFF94 strain in excess of the starting pose, never the absolute "
        "energy. Below zero means the ligand relaxed out of the pose it "
        "started in, which is a good outcome rather than a neutral one.")
    add("torsions", "Ligand torsions", ligand_torsions(measured),
        "Each rotatable torsion over the run. Flat traces are a rigid ligand; "
        "a trace jumping between bands is a rotamer flip.")
    add("ligand_rmsf", "Ligand atom fluctuation", ligand_atom_fluctuation(measured),
        "Which part of the ligand moves, about its own mean position.")
    add("pocket_rmsf", "Receptor flexibility", pocket_flexibility(measured),
        "Per-residue backbone fluctuation, with the binding-site residues "
        "marked. A mobile pocket makes every pose measurement noisier.",
        wide=True)
    # Wide, like every other heat map on the page. A square matrix in a
    # half-width card is drawn 340 pixels tall and then letterboxed to fit,
    # which wastes the width it was given and makes a block pattern harder to
    # read than the same data in the wide card the version 1 maps get.
    add("rmsd_2d", "Frame-to-frame RMSD", rmsd_matrix(measured),
        "Every frame against every other. Blocks along the diagonal are "
        "distinct states; a uniform field is a run that stayed in one.",
        wide=True)
    add("dccm", "Residue cross-correlation", cross_correlation(measured),
        "Which parts of the receptor move together. Blue is correlated, red is "
        "anticorrelated.", wide=True)
    add("energy_terms", "MM/GBSA, term by term", energy_terms(measured),
        "Where the binding energy comes from. The bars are the difference block, "
        "complex minus receptor minus ligand, which is the only part of an "
        "MM/GBSA calculation that is a binding quantity. Error bars are the "
        "standard deviation over the sampled frames, not the standard error. "
        "igb=2, single trajectory, no entropy term: good for ranking ligands "
        "against each other, not an affinity.")
    return panels


# ===========================================================================
#  The house style, for a sibling module that builds a panel
# ===========================================================================
#
# `flexappeal/refstates/panels.py` builds the reference-state panels that sit at
# the end of a run page. They have to look like every panel above them: the same
# recessive axes, the same unified hover, the same ink, the same decimation
# limits. These aliases are how it gets them.
#
# Public names for things this module has always had, rather than a second copy
# of the layout over there. A panel built from a different `_layout` is a panel
# that drifts the first time this one is touched, which is the thing section 6.1
# exists to prevent.
layout = _layout
line = _line
figure = _figure
finite = _finite
wire_stride = _wire_stride
dssp_states = _DSSP_STATES
