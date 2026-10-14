"""The reference-state panels, which were a tab of their own until 13 October.

Every claim in here was a test of that tab: four templates, four JSON endpoints
and a SQLite database behind them. The runs are rows in the Runs library now and
their measurements are panels at the end of a run page, so these follow the
panels rather than being deleted with the pages. Where a claim could not follow,
it is gone and `FINDINGS.md` 4bu says which and why.

They are unit tests on the builders and not browser tests. `tests/test_pages.py`
opens both of the fixture runs in Chromium at two viewports, runs axe over them
and asserts that each gated panel is drawn by one of them, which is the half
that needs a browser.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

import pytest

from flexappeal import fxa
from flexappeal.refstates import panels
from flexappeal.refstates.state import DECISIVE_S


# ---------------------------------------------------------------- fixtures
#
# Generated here, never committed, which is the rule for every fixture in this
# build. A payload is the shape `tools/pack_run_fxa.py` writes into
# `refstates.json`, with only the keys the panel under test reads: a fixture
# carrying the whole of a real payload would be a 2 MB file whose contents
# nobody could check against the thing being asserted.


def _frames(n: int = 12, motifs=("xDFG", "hinge"), *, kind: str = "kinase"):
    out = []
    for index in range(n):
        built = []
        for position, motif in enumerate(motifs):
            # A drift rather than a constant, so a half mean differs from the
            # value at the end and a test on one cannot pass on the other.
            score = 0.1 * position + index * 0.02
            built.append({
                "motif": motif, "n_residues": 4,
                "d_act_a": 1.0 + index * 0.01,
                "d_inact_a": 4.0 - index * 0.01,
                "s": round(score, 4),
                "label": "active" if score > DECISIVE_S else "intermediate",
                "nearest_active": "2F4J", "nearest_inactive": "3PYY",
            })
        out.append({"time_ns": index / 2.0, "n_conflicts": 0,
                    "whole_s": round(index * 0.02, 4),
                    "whole_label": "intermediate", "motifs": built,
                    "measurements": {
                        # The real shape, which is nested: the reference
                        # manifest's measurements are flat and a frame's are
                        # not, and a fixture that flattened them would let a
                        # panel reading the wrong path pass.
                        "dfg_plane": {"status": "ok",
                                      "d_alphac_a": 6.5 + index * 0.05,
                                      "d_lys_a": 14.0 - index * 0.04,
                                      "state": "in"},
                        "r_spine": {"state": "assembled"},
                        "alphac_salt_bridge": {
                            "distance_a": 3.10 + index * 0.02,
                            "state": "formed", "pair": "E286 - K271"},
                    }})
    return out


def _payload(**overrides):
    payload = {
        "name": "fixture", "kind": "kinase", "system_id": "abl1",
        "n_frames": 12, "frames_packed": 6, "smoothing_window_ns": 2.0,
        "core": {"n_residues": 43, "source": "fixture core"},
        "motif_residues": {"xDFG": [380, 381, 382, 383],
                           "hinge": [316, 317, 318]},
        "ensembles": {"active": ["2F4J"], "inactive": ["3PYY"]},
        "report": {
            "alignment": {"target": "2GQG", "n_aligned": 4, "n_dropped": 1,
                          "mean_core_rmsd_a": 0.49, "max_core_rmsd_a": 0.92},
            "rigidity": {"reference": "2GQG", "mean_rmsd_a": 0.31},
            "sse_profiles": {
                "active": {"xDFG": {"0": {"E": 4}, "1": {"E": 3, "C": 1}}},
                "inactive": {"xDFG": {"0": {"C": 4}, "1": {"C": 4}}}},
        },
        "dssp": {"n_frames": 12, "n_residues": 5,
                 "residues": [316, 317, 318, 380, 381],
                 "frames": ["HHHEE"] * 12,
                 "motifs": {
                     "xDFG": {"reference_agreement": {
                         "decisive": True, "leans": "active", "margin": 0.61,
                         "n_discriminating": 4, "n_compared": 4}},
                     "hinge": {"reference_agreement": {
                         "decisive": False, "leans": "inactive",
                         "margin": 0.04, "n_discriminating": 1,
                         "n_compared": 3,
                         "note": "only 1 residue(s) distinguish the "
                                 "reference sets here"}}}},
        "whole_s_raw": [index * 0.02 for index in range(12)],
        "whole_s_smoothed": [0.01 + index * 0.02 for index in range(12)],
        "frames": _frames(),
        "reference_structures": [
            {"pdb": "2F4J", "chain": "A", "resolution": 1.9, "source": "KLIFS",
             "source_state": "in", "our_state": "in",
             "metrics": {"dfg_d_alphac": 6.83, "dfg_d_lys": 13.63,
                         "their_dfg": "in", "our_dfg": "in",
                         "dfg_agrees": True}},
            {"pdb": "2HYY", "chain": "A", "resolution": 2.4, "source": "KLIFS",
             "source_state": "out", "our_state": "out",
             "metrics": {"dfg_d_alphac": 13.82, "dfg_d_lys": 11.14,
                         "their_dfg": "out", "our_dfg": "out",
                         "dfg_agrees": True}},
            {"pdb": "2V7A", "chain": "B", "resolution": 2.2, "source": "KLIFS",
             "source_state": "in", "our_state": None,
             "metrics": {"dfg_d_alphac": 6.71, "dfg_d_lys": 13.9,
                         "their_dfg": "in", "our_dfg": None,
                         "dfg_agrees": None}},
        ],
        "membrane_chunks": [],
        "notes": [],
        "viewer": {"motif_painted": True,
                   "overlay": {"active": {"member": "refstates_active.pdb",
                                          "code": "2HZ4",
                                          "mean_motif_rmsd_a": 4.06}},
                   "lipids": {}},
    }
    payload.update(overrides)
    return payload


def _results(payload, *, members=None, frames_packed=6):
    """A `Results` carrying this payload, built through the real reader."""
    body = {
        "metrics.json": json.dumps({"time_ns": [0.0, 0.5]}).encode(),
        "refstates.json": json.dumps(payload).encode(),
        **(members or {}),
    }
    manifest = {
        "fxa_version": 1, "job_name": payload.get("name", "fixture"),
        "trajectory": {"frames_packed": frames_packed, "atoms_packed": 100,
                       "duration_ns": 5.0},
        "files": {name: {"bytes": len(blob),
                         "sha256": hashlib.sha256(blob).hexdigest()}
                  for name, blob in body.items()},
    }
    import io

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("manifest.json", json.dumps(manifest))
        for name, blob in body.items():
            archive.writestr(name, blob)
    return fxa.load(buffer.getvalue(), verify_checksums=True)


def _by_id(built):
    return {panel["id"].removeprefix("refstates_"): panel for panel in built}


# ------------------------------------------------------------- the palette


def test_the_panel_colours_come_from_the_stylesheet():
    """Section 6.1 fixes them as CSS custom properties; CLAUDE.md says do not
    pick a colour in a chart.

    These figures are built on the server, where there is no computed style to
    read, so the tokens are parsed out of `brand.css` itself. One definition,
    in the file the browser also reads.
    """
    import re

    from flexappeal.refstates.store import MOTIF_TOKENS

    css = panels.BRAND_CSS.read_text(encoding="utf-8")
    declared = dict(re.findall(r"(--[a-z0-9-]+)\s*:\s*(#[0-9a-fA-F]{6})", css))

    # The value the stylesheet gives this motif's own token, from section 6.1's
    # table in `store.py`, rather than whichever colour happened to come out.
    token = MOTIF_TOKENS["kinase"]["xDFG"]
    assert panels.motif_colour("kinase", "xDFG") == declared[token]
    assert panels.state_colour("active") == declared["--state-active"]
    assert panels.state_colour("inactive") == declared["--state-inactive"]


def test_a_missing_token_raises_rather_than_defaulting(monkeypatch, tmp_path):
    """There is nothing sensible to fall back to.

    A chart that quietly substituted its own colour for a missing token would
    put a colour on the page that section 6.1 does not define, and a reader
    comparing it with a badge would be comparing two different things.
    """
    thin = tmp_path / "brand.css"
    thin.write_text(":root { --motif-1: #0072b2; }", encoding="utf-8")
    monkeypatch.setattr(panels, "BRAND_CSS", thin)
    monkeypatch.setattr(panels._tokens, "_cache", None, raising=False)
    with pytest.raises(RuntimeError) as raised:
        panels._tokens()
    assert "--state-active" in str(raised.value)
    monkeypatch.setattr(panels._tokens, "_cache", None, raising=False)


# --------------------------------------------------------- the state score


def test_the_whole_protein_series_is_drawn_raw_and_smoothed():
    """Section 5.6 wants the smoothed line for display and the raw stored.

    Both are drawn. The smoothed one alone would hide how much of the scatter
    the smoothing removed, which is the reader's to judge.
    """
    figure = panels.state_score(_payload())
    names = [trace["name"] for trace in figure["data"]]
    assert names == ["per frame", "smoothed"]
    assert figure["data"][0]["y"][:3] == [0.0, 0.02, 0.04]
    assert figure["data"][1]["y"][0] == 0.01


def test_the_score_plot_draws_and_labels_its_threshold():
    """Section 10: threshold lines labelled on every plot.

    Drawn at `DECISIVE_S` rather than at zero, which is 4ad: a score of 0.05 is
    not slightly active, it is undecided, and a line at zero invites exactly
    that reading.
    """
    layout = panels.state_score(_payload())["layout"]
    bands = [shape for shape in layout["shapes"] if shape["type"] == "rect"]
    assert len(bands) == 2
    assert {band["y0"] for band in bands} == {DECISIVE_S, -1.05}
    said = " ".join(note["text"] for note in layout["annotations"])
    assert f"+{DECISIVE_S}" in said and f"-{DECISIVE_S}" in said


def test_the_state_threshold_has_one_definition():
    """It was written out three times and named once.

    Every number the panels draw or print for this threshold has to come from
    `state.DECISIVE_S`, so a reader is never told one thing by a band and
    another by a sentence.
    """
    moved = 0.45
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(panels, "DECISIVE_S", moved)
        built = _by_id(panels.build(_payload()))
        bands = [shape for shape in built["state"]["figure"]["layout"]["shapes"]
                 if shape["type"] == "rect"]
        assert moved in {band["y0"] for band in bands}
        said = " ".join(note["text"] for note in
                        built["state"]["figure"]["layout"]["annotations"])
        assert f"+{moved}" in said
        # And the prose, which is where it was written out by hand before: the
        # blurb said "+0.2" beside bands drawn from the constant.
        assert f"±{moved}" in built["state"]["blurb"]
        assert f"±{moved}" in built["motif_table"]["table"]["note"]
    finally:
        monkeypatch.undo()


# ------------------------------------------------------- motif timelines


def test_every_motif_is_one_line_in_its_own_colour():
    """One panel with a line each, rather than one small multiple per motif.

    The specified chart was a panel per motif and this is 4u: the question a
    reader has is whether the motifs agree with each other, and six panels side
    by side answer it worse than six lines on one pair of axes.
    """
    figure = panels.motif_timelines(_payload())
    names = [trace["name"] for trace in figure["data"]]
    assert names == ["hinge", "xDFG"]
    colours = {trace["name"]: trace["line"]["color"] for trace in figure["data"]}
    assert colours["xDFG"] == panels.motif_colour("kinase", "xDFG")
    assert colours["hinge"] == panels.motif_colour("kinase", "hinge")
    assert colours["xDFG"] != colours["hinge"]


def test_a_motif_that_resolved_one_side_still_has_a_time_axis():
    """Every series a panel draws must be plottable or absent.

    The specific case: a motif whose score could not be computed for some
    frames. The gap is drawn as a gap rather than as a zero, and the trace
    still spans the run so the other motifs stay comparable against it.
    """
    payload = _payload()
    for frame in payload["frames"][4:8]:
        for motif in frame["motifs"]:
            if motif["motif"] == "hinge":
                motif["s"] = None
    figure = panels.motif_timelines(payload)
    hinge = next(t for t in figure["data"] if t["name"] == "hinge")
    assert len(hinge["x"]) == len(hinge["y"]) == 12
    assert hinge["y"][5] is None
    assert hinge["y"][0] is not None


def test_a_panel_lets_the_container_set_its_height():
    """Setting a height in the layout as well pinned it before the grid ran.

    Every panel is `autosize`, and the stylesheet owns the height: a figure
    that carries its own is the one that comes out the wrong size on a phone.
    """
    for panel in panels.build(_payload()):
        if panel["kind"] != "figure":
            continue
        layout = panel["figure"]["layout"]
        assert layout.get("autosize") is True, panel["id"]
        assert "height" not in layout, panel["id"]
        assert "width" not in layout, panel["id"]


# ----------------------------------------------------------- the motif table


def test_the_motif_table_shows_both_half_means():
    """A single number at the end cannot tell a motif that settled from one
    that was still moving when the run stopped.

    The committed figures carry these, and the page did too. Computed from the
    same frames the timelines are drawn from, so the two cannot disagree.
    """
    table = panels.motif_table(_payload())
    columns = table["columns"]
    assert "First half to second" in columns
    cell = table["rows"][0]["cells"][columns.index("First half to second")]
    assert " to " in cell["text"], cell
    first, second = (float(part) for part in cell["text"].split(" to "))
    assert second > first, "this fixture drifts upwards, so the halves differ"


def test_an_undecided_secondary_structure_lean_is_withheld():
    """A lean on a small margin must not appear as a finding.

    The payload says when it is decisive and when it is not, and the page is
    where that distinction either survives or is quietly dropped. `xDFG` is
    decisive here and `hinge` is not, so one prints and the other is blank.
    """
    table = panels.motif_table(_payload())
    column = table["columns"].index("Secondary structure leans")
    said = {row["cells"][0]["text"]: row["cells"][column]["text"]
            for row in table["rows"]}
    assert said["xDFG"] == "active"
    assert said["hinge"] == ""


def test_the_motif_table_never_calls_a_frame_final():
    """The same slip P5 made, in the column a reader looks at.

    A running row's latest frame is not its final one, and this table cannot
    know which it has, so it says what it does know: the last frame analysed.
    """
    table = panels.motif_table(_payload())
    assert "last frame analysed" in table["note"]
    assert "final" not in table["note"].lower()


def test_a_coloured_cell_always_carries_its_own_text():
    """Section 6.1: never by colour alone.

    A swatch in a table is decoration on a label a reader can also read, so the
    page says the same thing in greyscale.
    """
    for panel in panels.build(_payload(ligand={"resname": "STI"})):
        if panel["kind"] != "table":
            continue
        for row in panel["table"]["rows"]:
            for cell in row["cells"]:
                if cell.get("colour"):
                    assert cell.get("text"), (panel["id"], row)


# ------------------------------------------------------------- the ribbon


def test_the_dssp_ribbon_is_decimated_and_keeps_its_gaps():
    """An undecimated per-frame heatmap is what made one page 24 MB before.

    And a residue with no assignment stays a gap: filling it with coil would
    invent secondary structure the analysis did not find.
    """
    payload = _payload()
    payload["dssp"]["frames"] = ["HHHE-"] * 900
    payload["frames"] = _frames(900)
    payload["dssp"]["n_frames"] = 900
    figure = panels.dssp_ribbon(payload)
    heatmap = figure["data"][0]
    assert len(heatmap["x"]) <= panels.plots.WIRE_HEATMAP_COLUMNS // 2
    assert len(heatmap["z"]) == 5, "one row per residue"
    assert heatmap["z"][4][0] is None, "a dash is a gap, not coil"


def test_the_dssp_ribbon_carries_the_reference_consensus():
    """Section 6.2 asks for "active/inactive reference consensus rows".

    Two extra columns on the same residue axis, so a reader compares along a
    row. The agreement fraction travels with the code: a residue where every
    reference agrees and one where half do are not the same claim.
    """
    payload = _payload()
    consensus = panels.sse_consensus(payload, payload["dssp"]["residues"])
    assert consensus["active"][3]["code"] == "E"
    assert consensus["inactive"][3]["code"] == "C"
    assert consensus["active"][3]["agreement"] == 1.0
    assert consensus["active"][0] is None, "no consensus outside the motifs"

    figure = panels.dssp_ribbon(payload)
    # One heatmap for the run and one per state, each two x values wide so the
    # column has a width Plotly can derive.
    extra = [trace for trace in figure["data"]
             if trace.get("type") == "heatmap"][1:]
    assert len(extra) == 2
    assert all(len(trace["x"]) == 2 for trace in extra)
    said = " ".join(note["text"] for note in figure["layout"]["annotations"])
    assert "active" in said and "inactive" in said


def test_the_ribbon_packs_overlapping_motifs_into_lanes():
    """NPxxY sits inside TM7 on this receptor, and a band drawn over another
    band hides it.

    The motif bands need residue numbers rather than a span sentence, and the
    nested pair has to end up in different lanes.
    """
    payload = _payload(kind="gpcr", motif_residues={
        "TM7": list(range(305, 331)), "NPxxY": [322, 323, 324, 325, 326]})
    lanes = panels._motif_lanes(payload)
    assert {item["motif"] for item in lanes} == {"TM7", "NPxxY"}
    assert len({item["lane"] for item in lanes}) == 2, lanes


# ----------------------------------------------------------- the state map


def test_the_state_map_never_plots_a_point_it_cannot_place():
    """A point with no y renders as nothing while the panel looks fine.

    So a reference missing either coordinate is left out of the cloud, and the
    count in the reference table is what says how many there were.
    """
    payload = _payload()
    payload["reference_structures"].append(
        {"pdb": "9XXX", "source_state": "in", "metrics": {"dfg_d_alphac": 7.0}})
    figure = panels.state_map(payload)
    plotted = {code for trace in figure["data"]
               for code in (trace.get("text") or [])}
    assert "9XXX" not in plotted
    assert "2F4J" in plotted


def test_the_published_cut_is_drawn_once_and_named_in_angstrom():
    """It was 11.0 in the comparison rule and 11.0 again in the state map.

    The line on the kinase map is `kinase.DFG_PUBLISHED_A` and it names its own
    unit, which is section 10's rule for a threshold on a plot.
    """
    from flexappeal.refstates.kinase import DFG_PUBLISHED_A

    layout = panels.state_map(_payload())["layout"]
    lines = [shape for shape in layout["shapes"] if shape["type"] == "line"]
    assert len(lines) == 1
    assert lines[0]["x0"] == DFG_PUBLISHED_A
    said = layout["annotations"][0]["text"]
    assert "published" in said and "Å" in said
    assert f"{DFG_PUBLISHED_A:.1f}" in said


def test_the_receptor_map_says_its_threshold_was_fitted():
    """And the kinase's says it was published.

    A reader is entitled to know that one of the two lines was chosen to fit
    the cloud it is drawn over. The receptor's is refitted on this very
    reference set under leave-one-out, which is P3, so the label says fitted
    and carries the accuracy.
    """
    payload = _payload(system_id="b2ar", kind="gpcr", reference_structures=[
        {"pdb": "3P0G", "source_state": "active",
         "metrics": {"tm3_tm6": 14.2, "npxxy_to_inactive": 2.95}},
        {"pdb": "3SN6", "source_state": "active",
         "metrics": {"tm3_tm6": 14.9, "npxxy_to_inactive": 2.79}},
        {"pdb": "2RH1", "source_state": "inactive",
         "metrics": {"tm3_tm6": 9.1, "npxxy_to_inactive": 0.42}},
        {"pdb": "5D5A", "source_state": "inactive",
         "metrics": {"tm3_tm6": 9.4, "npxxy_to_inactive": 0.51}}])
    said = panels.state_map(payload)["layout"]["annotations"][0]["text"]
    assert said.startswith("fitted")
    assert "leave-one-out" in said


def test_the_path_carries_the_frame_the_viewer_should_show():
    """Not the analysis's own frame index.

    They are the same for these four rows and not the same in general: a
    results file is packed at a stride, so the panel scales the index onto what
    the file actually carries and the click handler has no arithmetic to get
    wrong. Twelve analysed frames and six packed, so the last point is 5.
    """
    figure = panels.state_map(_payload())
    path = next(trace for trace in figure["data"] if trace["name"] == "this run")
    assert path["customdata"][0] == 0
    assert path["customdata"][-1] == 5


# ---------------------------------------------------------------- bilayer


def test_the_bilayer_windows_come_from_the_constants():
    """The note said "60 to 68" and the badge read `APL_EXPECTED_A2`.

    Both windows are drawn, because a chunk between them is what section 5.8
    calls acceptable and not expected, which is a real result rather than a
    failure.
    """
    from flexappeal.refstates.membrane import (APL_ACCEPTABLE_A2,
                                               APL_EXPECTED_A2)

    payload = _payload(membrane_chunks=[
        {"time_ns": 0.0, "apl": 68.05, "thickness": 38.42},
        {"time_ns": 10.0, "apl": 66.23, "thickness": 39.35}])
    figure = panels.bilayer(payload)
    band = next(shape for shape in figure["layout"]["shapes"]
                if shape["type"] == "rect")
    assert (band["y0"], band["y1"]) == APL_EXPECTED_A2
    lines = {shape["y0"] for shape in figure["layout"]["shapes"]
             if shape["type"] == "line"}
    assert lines == set(APL_ACCEPTABLE_A2)
    said = " ".join(note["text"] for note in figure["layout"]["annotations"])
    assert f"{APL_EXPECTED_A2[1]:.0f}" in said and "acceptable" in said


def test_a_chunk_on_the_boundary_is_marked_outside_the_expected_window():
    """68.05 against a window ending at 68.0 is outside it, and read as inside.

    4bn, which was a rounding in a summary document and then the same rounding
    in the control that produced it. Here the chunk's own marker colour says
    which window it is in, and the hover carries two decimals.
    """
    payload = _payload(membrane_chunks=[
        {"time_ns": 0.0, "apl": 68.05, "thickness": 38.42},
        {"time_ns": 10.0, "apl": 66.23, "thickness": 39.35}])
    figure = panels.bilayer(payload)
    apl = figure["data"][0]
    assert apl["y"] == [68.05, 66.23]
    assert ".2f" in apl["hovertemplate"]
    colours = apl["marker"]["color"]
    assert colours[0] == panels.state_colour("intermediate"), "outside expected"
    assert colours[1] == panels.state_colour("active"), "inside expected"


def test_a_water_row_has_no_bilayer_panel():
    """Omitted rather than drawn empty, which is the rule for every panel here."""
    assert panels.bilayer(_payload()) is None
    assert "bilayer" not in _by_id(panels.build(_payload()))


# ------------------------------------------------------- the measurements


def test_the_measurement_table_shows_what_moved_and_omits_the_rest():
    """Driven by a table, so adding a measurement does not mean editing a page.

    A reading absent from both ends is left out rather than shown blank, so the
    table is what was measured rather than a list of what might have been.
    """
    table = panels.measurements_table(_payload())
    labels = [row["cells"][0]["text"] for row in table["rows"]]
    # The resolved label, which names the actual residues, in preference to the
    # written one: a row whose analysis could not name its residues still has
    # to say what it measured.
    assert "E286 - K271" in labels
    assert "R-spine" in labels
    assert not any("Ionic lock" in label for label in labels), (
        "a receptor measurement has no place on a kinase row")
    change = {row["cells"][0]["text"]: row["cells"][3]["text"]
              for row in table["rows"]}
    assert change["E286 - K271"].startswith("+"), change


def test_the_measurement_table_offers_the_sodium_rows():
    """Including the nearest approach, which is the honest form of a negative.

    Section 5.5's sodium pocket: two distances that are always measurable and
    one that depends on there being an ion in range. A row that cannot be
    measured is dropped rather than shown as zero.
    """
    payload = _payload(kind="gpcr", system_id="b2ar")
    for frame in payload["frames"]:
        frame["measurements"] = {"micro_switches": {"sodium_pocket": {
            "d_s3x39_a": 4.1, "d_n7x49_a": 5.2, "nearest_sodium_a": 7.8}}}
    labels = [row["cells"][0]["text"]
              for row in panels.measurements_table(payload)["rows"]]
    assert "D2.50 to S3.39" in labels
    assert "Nearest sodium to D2.50" in labels


def test_a_run_with_no_measurements_has_no_measurement_panel():
    payload = _payload()
    for frame in payload["frames"]:
        frame.pop("measurements")
    assert panels.measurements_table(payload) is None


# ------------------------------------------------------------- the ligand


def test_the_ligand_panel_is_offered_on_a_water_row():
    """It was inside the membrane block once, so a water row with a ligand
    would never have shown it.

    `abl1_imatinib` is exactly that row, and it is one of the four.
    """
    payload = _payload(ligand={"resname": "STI", "n_heavy_atoms": 37,
                               "final_pose_rmsd_a": 0.661,
                               "symmetry_corrected": True,
                               "contact_cutoff_a": 4.0})
    built = _by_id(panels.build(payload))
    assert "ligand" in built
    assert "bilayer" not in built


def test_the_ligand_panel_says_whether_the_pose_was_symmetry_corrected():
    """An uncorrected RMSD overstates drift, so a reader has to be told which."""
    table = panels.ligand_table({"ligand": {
        "resname": "STI", "final_pose_rmsd_a": 0.66,
        "symmetry_corrected": True, "contact_cutoff_a": 4.0}})
    said = {row["cells"][0]["text"]: row["cells"][1]["text"]
            for row in table["rows"]}
    assert said["Symmetry corrected"] == "yes"
    assert said["Contact cutoff"] == "4.0 Å"


def test_a_contact_at_the_cutoff_does_not_print_as_the_cutoff():
    """The closest approach carries two decimals.

    At one, a contact measured at 3.96 angstrom against a 4.0 cutoff prints as
    4.0 and reads as no contact at all.
    """
    table = panels.ligand_table({"ligand": {
        "resname": "CAU", "contact_cutoff_a": 4.0,
        "persistent_contacts": [{"residue": 113, "name": "ASP",
                                 "fraction_of_frames": 1.0,
                                 "closest_a": 3.96}]}})
    said = " ".join(cell["text"] for row in table["rows"]
                    for cell in row["cells"])
    assert "3.96" in said
    assert "4.00 Å" not in said


# --------------------------------------------------------- the references


def test_the_reference_table_counts_its_own_rows():
    """It said 72 active and 28 inactive where the ensembles held 84 and 32.

    A count written into a sentence goes stale and nothing can catch it, so
    every number in this note is computed from the rows beside it.
    """
    table = panels.references_table(_payload())
    assert len(table["rows"]) == 3
    assert "3 structures" in table["note"]
    # Two of the three carry a label on both sides; the third is the "no
    # opinion" case and is not counted as agreement either way.
    assert "2 of 2" in table["note"]
    verdicts = [row["cells"][-2]["text"] for row in table["rows"]]
    assert verdicts.count("agrees") == 2
    assert verdicts.count("no opinion") == 1


def test_the_reference_table_says_which_structures_were_scored():
    """A structure the alignment dropped is in the catalogue and not in the
    ensemble, and the table has to be able to say so."""
    table = panels.references_table(_payload())
    column = table["columns"].index("In the ensemble")
    said = {row["cells"][0]["text"]: row["cells"][column]["text"]
            for row in table["rows"]}
    assert said["2F4J"] == "scored"
    assert said["2HYY"] == "not used"


def test_the_core_panel_reports_the_alignment_it_was_built_on():
    """Every distance on the page is measured after superposing on this core,
    so how tight it is bounds what any of them can mean."""
    table = panels.core_table(_payload())
    said = {row["cells"][0]["text"]: row["cells"][1]["text"]
            for row in table["rows"]}
    assert said["Core residues"] == "43"
    assert said["Alignment target"] == "2GQG"
    assert said["Mean core RMSD"] == "0.49 Å"
    assert said["References dropped"] == "1"


# ------------------------------------------------------- the viewer controls


def test_the_viewer_only_offers_a_ghost_whose_structure_is_in_the_file():
    """The guard the bilayer toggle had and the overlay did not.

    A toggle whose structure 404s is worse than no toggle: it reads as a broken
    page rather than as a run that was packed without a ghost.
    """
    payload = _payload()
    with_ghost = panels.viewer_controls(
        _results(payload, members={"refstates_active.pdb": b"ATOM\n"}))
    assert set(with_ghost["overlay"]) == {"active"}
    assert with_ghost["overlay"]["active"]["code"] == "2HZ4"

    without = panels.viewer_controls(_results(payload))
    assert without["overlay"] == {}


def test_the_bilayer_toggle_needs_the_structure_it_draws():
    """A toggle whose structure is not in the file would 404 rather than draw.

    It used to need a trajectory as well, which was right when the bilayer was
    476 phosphorus atoms over 500 frames and wrong now that it is 24,752 atoms
    at one frame: the newer shape has no trajectory to require, and
    `test_an_older_bilayer_that_moves_still_works` covers the older one.
    """
    payload = _payload(viewer={"motif_painted": False, "overlay": {},
                               "lipids": {"topology": "lipids.pdb",
                                          "n_atoms": 24752}})
    offered = panels.viewer_controls(
        _results(payload, members={"lipids.pdb": b"ATOM\n"}))
    assert offered["lipids"]["n_atoms"] == 24752

    assert panels.viewer_controls(_results(payload))["lipids"] == {}


def test_the_cartoon_is_only_coloured_when_the_topology_was_painted():
    """A results file packed before this existed has no motif number in its
    B-factor column, and painting from a column of zeroes would wash the whole
    cartoon to one colour."""
    painted = panels.viewer_controls(_results(_payload()))
    assert painted["motif_painted"] is True

    plain = _payload()
    plain["viewer"] = {"motif_painted": False, "overlay": {}, "lipids": {}}
    assert panels.viewer_controls(_results(plain))["motif_painted"] is False


def test_a_results_file_with_no_reference_states_gets_no_panels():
    """Which is every run but these four, so the page simply has nothing extra
    to show and the viewer card renders as it does for any other run."""
    import io

    body = {"metrics.json": json.dumps({"time_ns": [0.0]}).encode()}
    manifest = {"fxa_version": 1, "job_name": "plain",
                "trajectory": {"frames_packed": 2},
                "files": {name: {"bytes": len(blob),
                                 "sha256": hashlib.sha256(blob).hexdigest()}
                          for name, blob in body.items()}}
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("manifest.json", json.dumps(manifest))
        for name, blob in body.items():
            archive.writestr(name, blob)
    results = fxa.load(buffer.getvalue(), verify_checksums=True)
    assert panels.from_results(results) == []
    assert panels.viewer_controls(results) == {}


# ------------------------------------------------------------ the whole set


def test_every_panel_has_a_title_a_blurb_and_a_kind():
    """The template switches on `kind` rather than on which key is present: a
    template that infers it fails silently the first time a builder returns an
    empty figure."""
    built = panels.build(_payload(
        ligand={"resname": "STI", "contact_cutoff_a": 4.0},
        membrane_chunks=[{"time_ns": 0.0, "apl": 67.0, "thickness": 38.9}]))
    assert len(built) == 10
    for panel in built:
        assert panel["kind"] in ("figure", "table"), panel["id"]
        assert panel["title"] and panel["blurb"], panel["id"]
        assert panel["id"].startswith("refstates_")
        assert ("figure" in panel) == (panel["kind"] == "figure")


def test_the_panels_read_the_payload_out_of_the_results_file():
    """One member, so a run in the library is one file: the deployed site
    serves four results files and nothing else."""
    built = panels.from_results(_results(_payload()))
    assert [panel["id"] for panel in built][:2] == ["refstates_state",
                                                    "refstates_motifs"]


# ----------------------------------------------------- the page they land on
#
# These two were tests of the removed tab's four pages. Their claim is about
# markup rather than about those pages, so they follow the panels onto the run
# page that now carries them. Server-side, with Flask's test client:
# `tests/test_pages.py` is where a browser and axe-core see the same page.


@pytest.fixture
def run_page(tmp_path):
    """A run in a library of its own, rendered by the real app."""
    from flexappeal.webapp import create_app

    runs = tmp_path / "web_runs" / "abl1_apo"
    runs.mkdir(parents=True)
    payload = _payload(ligand={"resname": "STI", "contact_cutoff_a": 4.0,
                               "symmetry_corrected": True,
                               "persistent_contacts": [
                                   {"residue": 315, "name": "THR",
                                    "fraction_of_frames": 0.98,
                                    "closest_a": 2.78}]},
                       membrane_chunks=[{"time_ns": 0.0, "apl": 67.0,
                                         "thickness": 38.9}])
    results = _results(payload)
    # Written back out as the archive the route will read.
    import io

    body = dict(results.members)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, blob in body.items():
            archive.writestr(name, blob)
    (runs / "results.fxa").write_bytes(buffer.getvalue())
    (runs / "run.json").write_text(json.dumps({"label": "fixture"}))

    app = create_app()
    app.config["RUNS_ROOT"] = str(tmp_path / "web_runs")
    app.config.update(TESTING=True)
    return app.test_client()


def test_every_referenced_asset_resolves(run_page):
    """A 404 on a script or stylesheet is the usual cause of a console error."""
    import re

    response = run_page.get("/runs/abl1_apo")
    assert response.status_code == 200
    html = response.get_data(as_text=True)
    missing = [ref for ref in
               set(re.findall(r'(?:src|href)="(/static/[^"]+)"', html))
               if run_page.get(ref).status_code != 200]
    assert not missing, f"referenced but not served: {missing}"


def test_the_page_meets_the_markup_rules_axe_would_check(run_page):
    """The subset of axe-core's serious rules that live in the markup.

    Not a substitute for axe-core, and FINDINGS says so. These are the ones a
    template can get wrong and a server can see: a table without header scope,
    an image without alt text, a page without a language. The panel tables
    added nine new ones to this page, which is why it is worth re-asserting
    here rather than trusting the browser pass alone.
    """
    import re

    html = run_page.get("/runs/abl1_apo").get_data(as_text=True)
    problems = []
    if not re.search(r"<html[^>]+lang=", html):
        problems.append("<html> has no lang")
    for image in re.findall(r"<img[^>]*>", html):
        if "alt=" not in image:
            problems.append("an img has no alt")
    # `<th(?=[\s>])` so the pattern does not also match `<thead>`.
    for header in re.findall(r"<th(?=[\s>])[^>]*>", html):
        if "scope=" not in header:
            problems.append(f"a th has no scope: {header}")
    assert not problems, problems


def test_the_state_map_axes_are_labelled_in_angstrom(run_page):
    """The unit symbol and the Greek letters, not ASCII transliterations.

    The labels this panel replaced were written `(A)` with `Calpha` and
    `alphaC` spelled out while the same page wrote the symbol everywhere else,
    so one plot carried a different convention from the panel beside it. The
    first version of the new panel reintroduced exactly that, and this guard is
    what caught it.
    """
    figure = panels.state_map(_payload())
    labels = [figure["layout"]["xaxis"]["title"]["text"],
              figure["layout"]["yaxis"]["title"]["text"]]
    for label in labels:
        assert "\u00c5" in label, f"{label!r} carries no Angstrom"
        assert "(A)" not in label, f"{label!r} uses ASCII A"
        assert "Calpha" not in label and "alphaC" not in label, \
            f"{label!r} spells out a Greek letter"


# ------------------------------------------------- the bilayer and the camera


def test_a_single_frame_bilayer_needs_only_its_topology():
    """The bilayer is the whole lipid now rather than its phosphorus.

    Phosphorus alone is 476 atoms and draws as a cloud of dots; the full chains
    are 24,752 heavy atoms, which is 2 MB for one frame and 48 MB for five
    hundred. So it is one frame and there is no trajectory to require.
    """
    payload = _payload(viewer={"motif_painted": False, "overlay": {},
                               "lipids": {"topology": "lipids.pdb",
                                          "n_atoms": 24752, "n_residues": 476,
                                          "n_frames": 1}})
    controls = panels.viewer_controls(
        _results(payload, members={"lipids.pdb": b"ATOM\n"}))
    assert controls["lipids"]["n_residues"] == 476
    assert "trajectory" not in controls["lipids"]


def test_an_older_bilayer_that_moves_still_works():
    """A results file packed when the bilayer was phosphorus over 500 frames.

    The topology is what is required; a trajectory is carried when the file has
    one, so an older file keeps the moving bilayer it was packed with.
    """
    payload = _payload(viewer={"motif_painted": False, "overlay": {},
                               "lipids": {"topology": "lipids.pdb",
                                          "trajectory": "lipids.xtc",
                                          "n_atoms": 476, "n_frames": 500}})
    controls = panels.viewer_controls(_results(payload, members={
        "lipids.pdb": b"ATOM\n", "lipids.xtc": b"\x00"}))
    assert controls["lipids"]["trajectory"] == "lipids.xtc"

    # And one whose trajectory did not make it into the file is offered as the
    # single frame it has, rather than as a control that 404s.
    half = panels.viewer_controls(
        _results(payload, members={"lipids.pdb": b"ATOM\n"}))
    assert "trajectory" not in half["lipids"]
    assert half["lipids"]["topology"] == "lipids.pdb"


def test_a_run_carries_its_camera_whether_or_not_it_has_a_bilayer():
    """The two runs of a system have to open facing the same way.

    It was gated on the bilayer when the only thing that could define a view
    was the membrane normal. A water run's view comes from the shape of its own
    core instead, and both runs of a system are matched to one of them by
    superposing their cores, so a kinase pair opens as one picture rather than
    two unrelated angles.
    """
    orientation = {"up": [0, 0, 1], "along": [1, 0, 0], "normal_axis": 2}
    membrane = _payload(viewer={"motif_painted": False, "overlay": {},
                                "orientation": orientation,
                                "lipids": {"topology": "lipids.pdb",
                                           "n_residues": 476, "n_frames": 1}})
    controls = panels.viewer_controls(
        _results(membrane, members={"lipids.pdb": b"ATOM\n"}))
    assert controls["orientation"]["up"] == [0, 0, 1]

    water = _payload(viewer={"motif_painted": False, "overlay": {},
                             "orientation": orientation, "lipids": {}})
    assert panels.viewer_controls(_results(water))["orientation"] == orientation

    # And a file with no camera in it leaves the viewer alone.
    plain = _payload(viewer={"motif_painted": False, "overlay": {},
                             "lipids": {}})
    assert panels.viewer_controls(_results(plain))["orientation"] == {}


def test_the_viewer_legend_names_every_motif_the_cartoon_colours():
    """The structure carries the motif number and Mol* colours from it.

    Without a key a reader sees an orange helix and a blue one with nothing on
    the page naming either. The colours are section 6.1's tokens, the same ones
    the panel tables use, so a motif looks the same wherever it appears.
    """
    controls = panels.viewer_controls(_results(_payload()))
    legend = {item["name"]: item["colour"] for item in controls["motifs"]}
    assert set(legend) == {"xDFG", "hinge"}
    assert legend["xDFG"] == panels.motif_colour("kinase", "xDFG")
    assert legend["hinge"] == panels.motif_colour("kinase", "hinge")
    spans = {item["name"]: item["residues"] for item in controls["motifs"]}
    assert spans["xDFG"] == "380 to 383"


def test_an_unpainted_structure_gets_no_legend():
    """A key to colours that are not there would be worse than no key."""
    payload = _payload()
    payload["viewer"] = {"motif_painted": False, "overlay": {}, "lipids": {}}
    assert panels.viewer_controls(_results(payload))["motifs"] == []


# --------------------------------------------------- the packer's own staging


def _packer():
    """`tools/pack_run_fxa.py`, loaded the way the ingest tests load theirs."""
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_pack_for_tests", root / "tools" / "pack_run_fxa.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pdb(rows):
    """A minimal PDB: (serial, name, resname, resSeq, element) per atom."""
    out = []
    for serial, (name, resname, number, element) in enumerate(rows, start=1):
        # Columns matter: 13-16 the atom name, 17 altLoc, 18-20 the residue
        # name, 22 the chain, 23-26 the residue number. The first version of
        # this helper left out altLoc, which shifted the residue name a column
        # left and made the test assert on a file the parser could not read.
        out.append(
            f"HETATM{serial:>5d} {name:<4s} {resname:>3s} A{number:>4d}"
            f"    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00"
            f"          {element:>2s}\n")
    return "".join(out) + "END\n"


def test_a_ligand_split_across_two_residues_is_put_back_together(tmp_path):
    """The built systems have every ligand in two residues, heavy and hydrogen.

    Imatinib is 37 and 31, carazolol 22 and 22. The simulation does not care
    and neither does the reference-state analysis, which selects the ligand by
    name. The binding measurements do: with the halves apart they took "the
    largest fragment", which is the heavy atoms with no hydrogens, and could
    not infer bond orders from them. No two-dimensional depiction, so the
    mobility inset was missing from every bound run while the committed example
    had one; no strain or torsion, so the conformation sub-score stayed
    pending; and a pose RMSD with no symmetry correction, which overstates
    drift.
    """
    module = _packer()
    path = tmp_path / "topology.pdb"
    path.write_text(_pdb([("N", "ALA", 1, "N"), ("CA", "ALA", 1, "C")]
                         + [("C1", "STI", 2, "C"), ("C2", "STI", 2, "C")]
                         + [("H1", "STI", 3, "H"), ("H2", "STI", 3, "H")]))
    module._merge_ligand_residues(path, ["STI"])

    numbers = {int(line[22:26]) for line in path.read_text().splitlines()
               if line.startswith("HETATM") and line[17:20].strip() == "STI"}
    assert numbers == {2}, "both halves should share one residue number"


def test_two_real_copies_of_a_ligand_are_left_alone(tmp_path):
    """Merging those would join two molecules into one.

    The split this repairs is a specific shape: two residues of the same name,
    next to each other, one carrying hydrogens and the other not. Two copies of
    a ligand have hydrogens in both and are not touched.
    """
    module = _packer()
    path = tmp_path / "topology.pdb"
    before = _pdb([("C1", "STI", 1, "C"), ("H1", "STI", 1, "H"),
                   ("C1", "STI", 2, "C"), ("H1", "STI", 2, "H")])
    path.write_text(before)
    module._merge_ligand_residues(path, ["STI"])
    assert path.read_text() == before


def test_a_ligand_in_one_residue_already_is_left_alone(tmp_path):
    module = _packer()
    path = tmp_path / "topology.pdb"
    before = _pdb([("C1", "BNZ", 1, "C"), ("H1", "BNZ", 1, "H")])
    path.write_text(before)
    module._merge_ligand_residues(path, ["BNZ"])
    assert path.read_text() == before


def test_staging_never_writes_through_to_the_simulation_output(tmp_path):
    """The staged topology is a copy, because the staging rewrites it.

    `_merge_ligand_residues` puts a ligand back into one residue in the staged
    copy. The first version of the staging hard-linked that file from
    `data/md/<row>/solvated.pdb` and then wrote through the link, which edited
    the simulation's own output: both bound rows lost the heavy-and-hydrogen
    split that is the evidence for the fault the merge works around.

    Asserted on the inode rather than on the content, because content is what
    the mistake preserved: the two files were identical afterwards, which is
    exactly why nothing noticed until a test that reads `solvated.pdb` failed.
    """
    import shutil

    source = tmp_path / "md" / "row"
    source.mkdir(parents=True)
    (source / "solvated.pdb").write_text(_pdb([("C1", "STI", 1, "C"),
                                               ("H1", "STI", 2, "H")]))
    staged = tmp_path / "work" / "row"
    staged.mkdir(parents=True)
    shutil.copy2(source / "solvated.pdb", staged / "topology.pdb")

    assert ((staged / "topology.pdb").stat().st_ino
            != (source / "solvated.pdb").stat().st_ino), "a copy, not a link"

    module = _packer()
    module._merge_ligand_residues(staged / "topology.pdb", ["STI"])
    merged = {int(line[22:26]) for line in
              (staged / "topology.pdb").read_text().splitlines()
              if line.startswith("HETATM")}
    original = {int(line[22:26]) for line in
                (source / "solvated.pdb").read_text().splitlines()
                if line.startswith("HETATM")}
    assert merged == {1}, "the staged copy is merged"
    assert original == {1, 2}, "the simulation's own output is untouched"


def test_the_packer_copies_the_topology_rather_than_linking_it():
    """The guard for the mistake above, in the code that makes the decision."""
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "pack_run_fxa.py").read_text(encoding="utf-8")
    staging = source[source.index("def _stage("):]
    staging = staging[: staging.index("\ndef ")]
    assert 'shutil.copy2(source / "solvated.pdb", run_dir / "topology.pdb")' in staging
    assert "hardlink_to(src)" not in staging, (
        "the topology is rewritten in place, so it cannot share an inode with "
        "the simulation's output")


def test_a_lipid_near_the_far_face_is_imaged_whole_rather_than_split():
    """A molecule is moved by one shift, not atom by atom.

    Imaging per atom rounds the atoms of one lipid in opposite directions when
    the lipid sits about half a box from the receptor, and the two halves land
    a box apart: on screen it is a lipid stretched across the whole picture.
    It does not show at the first frame, where every lipid in the patch is
    within 18 angstrom of the receptor, and it shows more the longer the run
    goes and the further the lipids have diffused.
    """
    import numpy as np

    packer = _packer()
    box = np.array([60.0, 60.0, 60.0])
    centre = np.array([30.0, 30.0, 30.0])
    # Four atoms of one molecule, whole, straddling the far face. Their
    # offsets from the receptor are 29 to 32 angstrom, so half of them round
    # to zero boxes and half to one.
    molecule = np.array([[59.0, 30.0, 30.0], [60.0, 30.0, 30.0],
                         [61.0, 30.0, 30.0], [62.0, 30.0, 30.0]])

    split = centre + (molecule - centre) - box * np.round((molecule - centre) / box)
    assert np.ptp(split[:, 0]) > 50, "per atom, this molecule is torn in two"

    whole = packer._image_whole(molecule, [np.arange(4)], centre, box)
    assert np.ptp(whole[:, 0]) == pytest.approx(3.0), "it stays one molecule"
    # And on the near side of the box, which is the point of imaging it.
    assert abs(whole[:, 0].mean() - centre[0]) < box[0] / 2


def test_the_molecules_are_grouped_by_residue():
    """One shift per residue means one group per residue, in order."""
    import types

    import numpy as np

    packer = _packer()
    # Atom index 7 and 9 are one lipid, 8 is another: the selection is not
    # contiguous and not sorted by residue, which is what a patch chosen by
    # distance looks like.
    residues = {7: 3, 8: 4, 9: 3}

    class _Topology:
        def atom(self, index):
            return types.SimpleNamespace(
                residue=types.SimpleNamespace(index=residues[index]))

    groups = packer._residue_groups(_Topology(), np.array([7, 8, 9]))
    assert [sorted(group.tolist()) for group in groups] == [[0, 2], [1]]


def test_both_bilayer_paths_image_whole_molecules():
    """The guard for the fault above, at both places that image lipids.

    The first frame and every later frame have to be imaged the same way, or
    the bilayer changes shape the moment a reader presses play.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "pack_run_fxa.py").read_text(encoding="utf-8")
    frames = source[source.index("def _stage_lipid_frames("):]
    frames = frames[: frames.index("\ndef ")]
    assert "_image_whole(" in frames
    assert "np.round(offset / box)" not in frames, "per atom, this tears lipids"

    members = source[source.index("def _viewer_members("):]
    first = members[members.index("def imaged("):]
    first = first[: first.index("\n        centre =")]
    assert "_image_whole(" in first
    assert "np.round(offset / box)" not in first


def test_no_panel_prints_its_own_heading_twice():
    """The card's <h3> names the panel, so the figure must not name it again.

    The title is taken from the figure, which keeps each builder able to
    describe itself when a figure is used on its own, and then blanked, which
    is what `plots.add` does for the panels above these. Before this, five
    charts carried the same words as a heading and as a chart title, one above
    the other, and the duplicate took vertical space from the plot.
    """
    panels_built = panels.build(_payload())
    charts = [panel for panel in panels_built if panel["kind"] == "figure"]
    assert charts, "the fixture should build at least one chart"
    for panel in charts:
        assert panel["title"], f"{panel['id']} lost its heading"
        assert panel["figure"]["layout"]["title"] == {"text": ""}, (
            f"{panel['id']} still carries an in-figure title")
        # And the space the title took goes back to the plot.
        assert panel["figure"]["layout"]["margin"]["t"] == 12


def test_the_packer_sends_the_ligand_chemistry_into_the_bundle():
    """Without it there is no binding energy, and nothing said so.

    `bundle.build` takes the SDFs from its caller. The web app has called
    `collect_ligands` before `build` since the ligand path was written; this
    packer is the second caller and did not, so it built bundles with an empty
    `ligands/` directory. The analysis pass then declined MM/GBSA with "holds
    no .sdf, so the ligand cannot be parameterised", on stdout that was thrown
    away unless the pass failed, and both ligand runs went to the site reading
    "no binding energy in this file".
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "pack_run_fxa.py").read_text(encoding="utf-8")
    assert "bundle.collect_ligands(cfg)" in source
    assert "ligands=ligands" in source
    collect = source.index("bundle.collect_ligands(cfg)")
    assert collect < source.index("bundle.build("), (
        "the chemistry has to be fetched before the bundle is built")


def test_the_packer_keeps_the_analysis_output_even_when_it_succeeds():
    """A step that declines to run is what the log is for.

    The output was printed only when the pass exited non-zero, so a pass that
    succeeded while skipping the energy said nothing at all.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "pack_run_fxa.py").read_text(encoding="utf-8")
    block = source[source.index("def _analyse("):]
    block = block[: block.index("\ndef ")]
    assert "-analyse.log" in block
    assert "analysis warning" in block, (
        "a warning from the pass should reach the packer's own log")
