"""Results-file reading, plotting, and the ways an upload can be hostile.

The .fxa arrives from an unauthenticated upload form, so every test here that
looks paranoid is testing a real attack shape: zip traversal, zip bombs, JSON
that is valid but not a mapping, and a manifest that lies about its contents.

The happy-path fixture is a real .fxa produced by an actual OpenMM run of
lysozyme, not a synthetic one -- see tests/fixtures/lysozyme.fxa.
"""

from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pytest

from flexappeal import fxa, plots
from flexappeal.fxa import FxaError

from conftest import FIXTURES

FIXTURE = FIXTURES / "lysozyme.fxa"


@pytest.fixture(scope="module")
def raw():
    if not FIXTURE.is_file():
        pytest.skip("tests/fixtures/lysozyme.fxa is missing")
    return FIXTURE.read_bytes()


@pytest.fixture(scope="module")
def results(raw):
    return fxa.load(raw)


def _make_fxa(members: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)
    return buffer.getvalue()


def _minimal(**manifest_overrides) -> bytes:
    manifest = {"fxa_version": 1, "job_name": "test", "metrics_computed": ["rmsd"]}
    manifest.update(manifest_overrides)
    return _make_fxa({
        "manifest.json": json.dumps(manifest).encode(),
        "metrics.json": json.dumps({"time_ns": [0, 1], "rmsd_nm": [0.0, 0.1]}).encode(),
    })


# ---------------------------------------------------------------------------
#  The real file
# ---------------------------------------------------------------------------


def test_loads_a_real_results_file(results):
    assert results.job_name
    assert results.metrics["time_ns"]
    assert results.topology_pdb, "the viewer needs a topology"
    assert results.trajectory_xtc, "the viewer needs coordinates"


def test_checksums_verify_on_an_untouched_file(results):
    assert not results.warnings, f"an intact file should not warn: {results.warnings}"


def test_summary_pulls_the_headline_facts(results):
    summary = fxa.summarise(results)
    assert summary["duration_ns"] > 0
    assert summary["force_field"]
    assert summary["platform"]
    assert summary["openmm_version"]


def test_arrays_load_as_numpy(results):
    if "contact_map" in results.arrays:
        matrix = results.array("contact_map")
        assert isinstance(matrix, np.ndarray)
        assert matrix.ndim == 2
        assert matrix.shape[0] == matrix.shape[1], "a contact map must be square"


def test_a_modified_member_is_flagged(raw):
    """The checksums exist so a corrupted transfer is visible rather than silent."""
    original = fxa.load(raw)
    members = dict(original.members)
    members["metrics.json"] = json.dumps(original.metrics).encode() + b"\n"
    tampered = _make_fxa(members)

    loaded = fxa.load(tampered)
    assert any("checksum" in w for w in loaded.warnings)


# ---------------------------------------------------------------------------
#  Malformed input
# ---------------------------------------------------------------------------


def test_empty_file():
    with pytest.raises(FxaError, match="empty"):
        fxa.load(b"")


def test_not_a_zip():
    with pytest.raises(FxaError, match="not a valid"):
        fxa.load(b"this is a trajectory, honest" * 100)


def test_missing_required_members():
    with pytest.raises(FxaError, match="missing"):
        fxa.load(_make_fxa({"manifest.json": b"{}"}))


def test_corrupt_json():
    with pytest.raises(FxaError, match="corrupt"):
        fxa.load(_make_fxa({"manifest.json": b"{not json", "metrics.json": b"{}"}))


@pytest.mark.parametrize("body", [b"3", b"[1,2,3]", b'"a string"', b"null", b"true"])
def test_metadata_that_is_valid_json_but_not_an_object(body):
    """json.loads("3") is an int; .get on it is an AttributeError, not a 400."""
    with pytest.raises(FxaError, match="expected form"):
        fxa.load(_make_fxa({"manifest.json": body, "metrics.json": b"{}"}))


def test_a_newer_format_version_is_refused_with_advice():
    with pytest.raises(FxaError, match="newer FlexAppeal"):
        fxa.load(_minimal(fxa_version=fxa.FXA_VERSION + 1))


def test_a_missing_version_warns_rather_than_failing():
    loaded = fxa.load(_make_fxa({
        "manifest.json": json.dumps({"job_name": "x"}).encode(),
        "metrics.json": b"{}",
    }))
    assert any("version" in w for w in loaded.warnings)


def test_an_invalid_version_type_is_refused():
    with pytest.raises(FxaError, match="invalid format version"):
        fxa.load(_minimal(fxa_version="one"))


# ---------------------------------------------------------------------------
#  Hostile archives
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "../../../etc/passwd",
    "/etc/passwd",
    "..\\..\\windows\\system32\\config",
    "C:/windows/system32",
])
def test_path_traversal_is_refused(name):
    """Nothing is written to disk here, but the check must hold regardless.

    /reanalyse will extract to a scratch directory, and a validator that only
    works because the current caller happens not to write is not a validator.
    """
    with pytest.raises(FxaError, match="unsafe path"):
        fxa.load(_make_fxa({name: b"x", "manifest.json": b"{}", "metrics.json": b"{}"}))


def test_too_many_entries_is_refused():
    members = {f"file{i}.txt": b"x" for i in range(fxa.MAX_ENTRIES + 5)}
    members["manifest.json"] = b"{}"
    members["metrics.json"] = b"{}"
    with pytest.raises(FxaError, match="entries"):
        fxa.load(_make_fxa(members))


def test_a_zip_bomb_is_refused():
    """A member that expands enormously from almost nothing."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        archive.writestr("manifest.json", b"{}")
        archive.writestr("metrics.json", b"{}")
        archive.writestr("bomb.bin", b"\0" * (60 * 1024 * 1024))
    with pytest.raises(FxaError, match="compression ratio"):
        fxa.load(buffer.getvalue())


def test_an_oversized_member_is_refused(monkeypatch):
    """Size alone must trip the guard, independently of compression ratio.

    The cap is lowered rather than a 512 MB member built, so the test stays
    fast; zipfile recomputes file_size on write, so a forged header is not an
    option.
    """
    monkeypatch.setattr(fxa, "MAX_SINGLE_MEMBER", 1024)
    payload = bytes(range(256)) * 8  # 2 KB of incompressible-ish data
    with pytest.raises(FxaError, match="larger than"):
        fxa.load(_make_fxa({
            "manifest.json": b"{}", "metrics.json": b"{}", "big.bin": payload,
        }))


def test_a_pickled_object_array_is_refused_not_executed():
    """np.load with allow_pickle=True would make an .npz arbitrary code execution.

    numpy refuses an object array outright when allow_pickle is off, which is
    the property being pinned: the load fails loudly and the page still renders
    rather than the pickle being deserialised.
    """
    payload = io.BytesIO()
    np.savez(payload, evil=np.array([{"anything": 1}], dtype=object))

    loaded = fxa.load(_make_fxa({
        "manifest.json": json.dumps({"fxa_version": 1}).encode(),
        "metrics.json": b"{}",
        "arrays.npz": payload.getvalue(),
    }))
    assert loaded.arrays == {}, "nothing from a pickled npz should reach the page"
    assert any("arrays" in w for w in loaded.warnings)


def test_a_non_npy_member_is_not_deserialised():
    """A member that is not a .npy comes back as opaque bytes, never an object."""
    import pickle

    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr("evil.npy", pickle.dumps({"anything": 1}))

    loaded = fxa.load(_make_fxa({
        "manifest.json": json.dumps({"fxa_version": 1}).encode(),
        "metrics.json": b"{}",
        "arrays.npz": payload.getvalue(),
    }))
    assert not isinstance(loaded.arrays.get("evil"), dict)


# ---------------------------------------------------------------------------
#  Plots
# ---------------------------------------------------------------------------


def test_every_panel_is_json_serialisable(results):
    json.dumps(plots.build_all(results))


def test_the_real_file_produces_the_expected_panels(results):
    ids = {p["id"] for p in plots.build_all(results)["panels"]}
    assert {"rmsd", "rgyr", "rmsf"} <= ids


def test_panels_have_no_duplicate_title(results):
    """The card heading names the panel; an in-figure title would repeat it."""
    for panel in plots.build_all(results)["panels"]:
        assert panel["figure"]["layout"]["title"]["text"] == ""
        assert panel["title"], "the card still needs a heading"


def test_missing_metrics_produce_no_panel():
    empty = fxa.Results(manifest={"fxa_version": 1}, metrics={})
    assert plots.build_all(empty)["panels"] == []


def test_single_series_line_panels_carry_no_legend(results):
    """A legend box for one line is noise -- the heading already names it.

    Heatmap panels are excluded: their categories live in the z values rather
    than in separate traces, so they need a legend even though they draw as one
    trace plus zero-width legend proxies.
    """
    for panel in plots.build_all(results)["panels"]:
        traces = panel["figure"]["data"]
        if any(t.get("type") == "heatmap" for t in traces):
            continue
        drawn = [t for t in traces if t.get("x") and t["x"] != [None]]
        if len(drawn) == 1:
            assert not panel["figure"]["layout"].get("showlegend"), \
                f"{panel['id']} shows a legend for a single series"


def test_multi_category_panels_do_carry_a_legend(results):
    """Identity must never be colour-alone."""
    panels = {p["id"]: p for p in plots.build_all(results)["panels"]}
    for key in ("dssp", "ss_fractions"):
        if key in panels:
            assert panels[key]["figure"]["layout"].get("showlegend"), \
                f"{key} encodes category by colour and needs a legend"


def test_no_panel_uses_a_second_y_axis(results):
    """A dual-axis chart invents a correlation that is not in the data."""
    for panel in plots.build_all(results)["panels"]:
        assert "yaxis2" not in panel["figure"]["layout"], \
            f"{panel['id']} has a second y-axis"
        for trace in panel["figure"]["data"]:
            assert trace.get("yaxis") in (None, "y")


def test_convergence_is_small_multiples_not_one_crowded_plot(results):
    columns = plots.parse_state_data(results.members["state_data.csv"])
    figures = plots.convergence(columns)
    assert len(figures) >= 3, "expected one figure per measure"
    for figure in figures:
        assert len(figure["figure"]["data"]) == 1, \
            "each measure gets its own plot; scales differ by orders of magnitude"


def test_state_data_parses_openmm_header(results):
    columns = plots.parse_state_data(results.members["state_data.csv"])
    # OpenMM writes '#"Step","Time (ps)",...' -- the leading # and the quotes
    # both have to come off or every column name is wrong.
    assert any(name.startswith("Time") for name in columns), list(columns)
    assert all(not name.startswith(("#", '"')) for name in columns)


def test_parse_state_data_survives_a_truncated_csv():
    """A run killed mid-write leaves a partial final line."""
    truncated = b'#"Step","Time (ps)","Temperature (K)"\n1,0.004,300.0\n2,0.008\n'
    columns = plots.parse_state_data(truncated)
    assert len(columns["Time (ps)"]) == 1


def test_parse_state_data_on_rubbish():
    assert plots.parse_state_data(b"") == {}


def test_tiles_degrade_when_metrics_are_absent():
    tiles = plots.tiles({"duration_ns": 5.0}, {})
    assert any(t["label"] == "Simulated" for t in tiles)
    assert all(t["value"] is not None for t in tiles)


def test_nan_becomes_a_gap_not_a_spike():
    """An exploded frame writes nan; plotting it as 0 would be a lie."""
    metrics = {"time_ns": [0, 1, 2], "rmsd_nm": [0.1, float("nan"), 0.3]}
    figure = plots.rmsd(metrics)
    assert figure["data"][0]["y"][1] is None


# ---------------------------------------------------------------------------
#  A real protein-ligand run
# ---------------------------------------------------------------------------

LIGAND_FIXTURE = FIXTURES / "trypsin_ben.fxa"


@pytest.fixture(scope="module")
def trypsin():
    if not LIGAND_FIXTURE.is_file():
        pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")
    return fxa.load(LIGAND_FIXTURE.read_bytes())


def test_ligand_run_produces_ligand_panels(trypsin):
    ids = {p["id"] for p in plots.build_all(trypsin)["panels"]}
    assert {"ligand_rmsd", "ligand_contacts"} <= ids


def test_the_ligand_is_one_residue_not_two(trypsin):
    """Chem.AddHs leaves new hydrogens without residue metadata.

    Left unstamped they form a separate UNK residue, and the symptom in the
    results is a contact list computed against 9 of the ligand's 17 atoms.
    """
    topology = trypsin.topology_pdb.decode()
    assert "UNK" not in topology, "ligand hydrogens ended up in their own residue"
    assert topology.count(" BEN ") >= 17, "the ligand should carry its hydrogens"


def test_benzamidine_finds_the_trypsin_s1_pocket(trypsin):
    """The chemistry test: does the pipeline reproduce known binding?

    Benzamidine in trypsin is textbook -- its amidinium forms a salt bridge with
    Asp189 at the bottom of the S1 specificity pocket. If parameterisation,
    placement or contact analysis were wrong, this is what would break, and no
    amount of "it ran without error" would tell us.
    """
    contacts = {c["residue"]: c["occupancy"] for c in trypsin.metrics["ligand_contacts"]}

    assert "ASP189" in contacts, f"the defining S1 contact is absent: {list(contacts)[:10]}"
    assert contacts["ASP189"] > 0.9, "Asp189 should be in contact essentially always"

    # The rest of the canonical S1 pocket.
    residue_numbers = {c["residue"][3:] for c in trypsin.metrics["ligand_contacts"]}
    assert {"189", "190", "215", "216", "219"} <= residue_numbers


def test_the_ligand_stays_in_the_pocket(trypsin):
    """A tight binder should not drift over a short run."""
    rmsd_angstrom = [v * 10 for v in trypsin.metrics["ligand_rmsd_nm"]]
    assert max(rmsd_angstrom) < 3.0, "the ligand left the site, or alignment is wrong"


def test_ligand_contacts_are_sorted_by_occupancy(trypsin):
    occupancies = [c["occupancy"] for c in trypsin.metrics["ligand_contacts"]]
    assert occupancies == sorted(occupancies, reverse=True)


def test_ligand_contact_chart_uses_one_colour_for_every_bar(trypsin):
    """One series: colouring by value would re-encode what bar length shows."""
    panels = {p["id"]: p for p in plots.build_all(trypsin)["panels"]}
    marker = panels["ligand_contacts"]["figure"]["data"][0]["marker"]
    assert isinstance(marker["color"], str), "a per-bar colour list is a value ramp"


# ---------------------------------------------------------------------------
#  A real membrane run
# ---------------------------------------------------------------------------

MEMBRANE_FIXTURE = FIXTURES / "glycophorin.fxa"


@pytest.fixture(scope="module")
def membrane():
    if not MEMBRANE_FIXTURE.is_file():
        pytest.skip("tests/fixtures/glycophorin.fxa is missing")
    return fxa.load(MEMBRANE_FIXTURE.read_bytes())


def test_membrane_run_produces_membrane_panels(membrane):
    ids = {p["id"] for p in plots.build_all(membrane)["panels"]}
    assert {"membrane_apl", "membrane_thickness", "membrane_scd"} <= ids


def test_the_bilayer_actually_got_built(membrane):
    assert membrane.has_membrane
    assert membrane.metrics["lipid_count"] > 50, "a bilayer needs lipids in both leaflets"


def test_lipid_selection_survives_pdb_name_truncation(membrane):
    """A PDB residue name is three characters, so POPC is stored as POP.

    Selecting on the configured name alone matches nothing and every membrane
    metric silently returns empty -- which is exactly what happened before this
    was fixed, with no error anywhere.
    """
    assert membrane.metrics.get("area_per_lipid_nm2"), \
        "the lipid selection matched nothing"
    assert membrane.metrics.get("lipid_order_parameters"), \
        "the lipid selection matched nothing"


def test_bilayer_thickness_matches_a_popc_bilayer(membrane):
    """POPC phosphate-to-phosphate is about 3.7-4.0 nm."""
    thickness = membrane.metrics["bilayer_thickness_nm"]
    mean = sum(thickness) / len(thickness)
    assert 3.0 < mean < 4.5, f"{mean:.2f} nm is not a POPC bilayer"


def test_lipid_order_parameters_have_the_fluid_bilayer_shape(membrane):
    """A fluid bilayer plateaus near the headgroup and falls toward the tail.

    This is the shape that says the lipids are behaving like a membrane rather
    than a frozen slab or a disordered mess, and it is the check that would
    catch order parameters computed against the wrong atoms.
    """
    profile = membrane.metrics["lipid_order_parameters"]
    assert len(profile) >= 8, "too few carbons resolved to judge the profile"

    values = {p["carbon"]: p["scd"] for p in profile}
    plateau = [values[c] for c in sorted(values) if c <= 8]
    tail = [values[c] for c in sorted(values) if c >= 14]

    assert all(0.05 < v < 0.45 for v in plateau), \
        f"plateau out of range for a fluid POPC bilayer: {plateau}"
    assert sum(tail) / len(tail) < sum(plateau) / len(plateau), \
        "order should fall toward the tail end; a flat or rising profile is wrong"


def test_area_per_lipid_is_labelled_as_a_diagnostic(membrane):
    """The absolute value includes the protein cross-section, so it must not be
    presented as a measurement to compare against a pure-bilayer value."""
    assert "protein" in membrane.metrics.get("area_per_lipid_note", "")
    panels = {p["id"]: p for p in plots.build_all(membrane)["panels"]}
    assert "cross-section" in panels["membrane_apl"]["blurb"]


def test_a_soluble_run_has_no_membrane_panels(results):
    ids = {p["id"] for p in plots.build_all(results)["panels"]}
    assert not any(i.startswith("membrane") for i in ids)
