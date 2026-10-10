"""Server-side re-analysis: request validation, budgets, locking, and physics.

This is the only code that does science on the droplet, on input from an
unauthenticated form, so most of these tests are about refusing things. The
physics tests exist because the one bug that mattered here was silent: periodic
wrapping made a 4 A contact measure 65 A while every other panel looked fine.
"""

from __future__ import annotations

import json
import shutil

import pytest

from flexappeal import analysis, plots
from flexappeal.analysis import Busy, ReanalysisError

from conftest import FIXTURES


@pytest.fixture
def session(tmp_path):
    """A scratch session holding a real results file, as the web app lays it out."""
    root = tmp_path / "scratch"
    (root / "sess").mkdir(parents=True)
    shutil.copy(FIXTURES / "trypsin_ben.fxa", root / "sess" / "results.fxa")
    return root


# ---------------------------------------------------------------------------
#  Request validation
# ---------------------------------------------------------------------------


def test_a_valid_request_parses():
    request = analysis.parse_request(
        {"metrics": ["rmsd", "rgyr"], "selection": "protein and resSeq 40 to 60"})
    assert request.metrics == ["rmsd", "rgyr"]


def test_at_least_one_metric_is_required():
    with pytest.raises(ReanalysisError, match="at least one"):
        analysis.parse_request({"metrics": [], "selection": "protein"})


def test_unknown_metrics_are_named_in_the_error():
    with pytest.raises(ReanalysisError, match="nonsense"):
        analysis.parse_request({"metrics": ["nonsense"], "selection": "protein"})


def test_a_distance_needs_two_selections():
    with pytest.raises(ReanalysisError, match="second selection"):
        analysis.parse_request({"metrics": ["distance"], "selection": "protein"})


def test_a_request_that_is_not_an_object_is_refused():
    for raw in (3, [1, 2], "protein", None):
        with pytest.raises(ReanalysisError):
            analysis.parse_request(raw)


@pytest.mark.parametrize("hostile", [
    "__import__('os').system('id')",
    "protein or __import__('os')",
    "(1).__class__.__bases__",
    "eval('1+1')",
    "protein\nimport os",
    "x" * 300,
    "",
])
def test_hostile_selections_are_refused(hostile):
    """MDTraj's parser rejects these too; this is the layer in front of it.

    The parser is a pyparsing grammar and does not evaluate Python, verified
    directly -- but it does compile what it parses, so a cheap gate in front
    costs nothing and does not depend on that staying true.
    """
    with pytest.raises(ReanalysisError):
        analysis.parse_request({"metrics": ["rmsd"], "selection": hostile})


def test_a_decimal_in_a_selection_is_still_allowed():
    """The dunder/attribute gate must not break `mass > 12.5`."""
    request = analysis.parse_request(
        {"metrics": ["rmsd"], "selection": "protein and mass > 12.5"})
    assert "12.5" in request.selection


# ---------------------------------------------------------------------------
#  Budgets
# ---------------------------------------------------------------------------


def test_an_oversized_job_is_refused_before_any_array_is_allocated():
    with pytest.raises(ReanalysisError, match="atom-frames"):
        analysis.check_budget(50_000, 2_000)


def test_a_reasonable_job_passes():
    analysis.check_budget(5_000, 500, 300, ["rmsd"])


def test_a_huge_contact_map_is_refused_separately():
    """Pairwise work is quadratic in residues, so the atom budget does not cover it."""
    with pytest.raises(ReanalysisError, match="contact map"):
        analysis.check_budget(1_000, 100, 2_000, ["contacts"])


# ---------------------------------------------------------------------------
#  Concurrency
# ---------------------------------------------------------------------------


def test_only_one_job_runs_at_a_time(tmp_path):
    analysis.acquire_lock(tmp_path)
    with pytest.raises(Busy):
        analysis.acquire_lock(tmp_path)
    analysis.release_lock(tmp_path)
    analysis.acquire_lock(tmp_path)          # free again
    analysis.release_lock(tmp_path)


def test_a_stale_lock_expires(tmp_path):
    """A worker killed mid-job must not wedge the slot forever."""
    import os
    import time

    lock = analysis.acquire_lock(tmp_path)
    old = time.time() - 10_000
    os.utime(lock, (old, old))
    analysis.acquire_lock(tmp_path, ttl=60)   # takes it over rather than raising
    analysis.release_lock(tmp_path)


def test_releasing_a_lock_that_is_not_held_is_harmless(tmp_path):
    analysis.release_lock(tmp_path)


# ---------------------------------------------------------------------------
#  The work itself
# ---------------------------------------------------------------------------


def test_computes_a_metric_the_run_did_not(session):
    """The whole point: SASA was never requested at build time."""
    request = analysis.parse_request({"metrics": ["sasa"], "selection": "protein"})
    metrics = analysis.run(session / "sess" / "results.fxa", request)
    assert len(metrics["sasa_total_nm2"]) == metrics["n_frames"]
    assert all(v > 0 for v in metrics["sasa_total_nm2"])


def test_restricts_a_metric_to_part_of_the_molecule(session):
    request = analysis.parse_request(
        {"metrics": ["rmsd"], "selection": "protein and resSeq 180 to 230"})
    metrics = analysis.run(session / "sess" / "results.fxa", request)
    assert metrics["n_atoms"] < 1000, "the selection should be a fraction of the protein"
    assert metrics["rmsd_nm"][0] == pytest.approx(0.0, abs=1e-6)


def test_a_selection_that_matches_nothing_says_so(session):
    request = analysis.parse_request(
        {"metrics": ["rmsd"], "selection": "resname ZZZ"})
    with pytest.raises(ReanalysisError, match="matched no atoms"):
        analysis.run(session / "sess" / "results.fxa", request)


def test_an_invalid_selection_gets_a_usable_message(session):
    request = analysis.Request(metrics=["rmsd"], selection="protein and and")
    with pytest.raises(ReanalysisError, match="not a valid MDTraj selection"):
        analysis.run(session / "sess" / "results.fxa", request)


def test_periodic_wrapping_is_undone_before_measuring(session):
    """The bug this test exists for: a bound ligand a whole box length away.

    Trajectories are written with enforcePeriodicBox, which wraps each molecule
    into the primary cell independently. Benzamidine sits in trypsin's S1
    pocket about 7 A from Asp189's alpha carbon; wrapped, it measures ~68 A --
    one box length. compute_contacts applies the minimum image convention
    itself, so the contact map looked perfectly correct while the distance
    beside it was nonsense.
    """
    request = analysis.parse_request({
        "metrics": ["distance"],
        "selection": "resname BEN",
        "selection_b": "resSeq 189 and name CA",
    })
    metrics = analysis.run(session / "sess" / "results.fxa", request)
    mean = sum(metrics["distance_nm"]) / len(metrics["distance_nm"]) * 10

    assert mean < 15, (
        f"benzamidine measures {mean:.1f} A from Asp189 CA -- it is bound in the "
        f"S1 pocket, so this is periodic wrapping, not chemistry"
    )
    assert mean > 3, "closer than a covalent bond; something else is wrong"


def test_reimaging_names_its_anchor(session):
    """MDTraj's own anchor heuristic cannot work for protein-plus-ligand.

    Left to choose, it looks for a molecule with *more* atoms than the largest,
    which is unsatisfiable, and raises rather than re-imaging.
    """
    import mdtraj as md

    from flexappeal import fxa

    results = fxa.load((FIXTURES / "trypsin_ben.fxa").read_bytes(),
                       verify_checksums=False)
    work = session / "sess"
    (work / "t.pdb").write_bytes(results.topology_pdb)
    (work / "t.xtc").write_bytes(results.trajectory_xtc)
    traj = md.load(str(work / "t.xtc"), top=str(work / "t.pdb"))

    with pytest.raises(ValueError, match="anchor molecules"):
        traj.image_molecules(inplace=False)

    assert analysis._reimage(traj) is True


def test_distance_between_two_protein_selections(session):
    request = analysis.parse_request({
        "metrics": ["distance"],
        "selection": "resSeq 189 and name CA",
        "selection_b": "resSeq 195 and name CA",
    })
    metrics = analysis.run(session / "sess" / "results.fxa", request)
    mean = sum(metrics["distance_nm"]) / len(metrics["distance_nm"]) * 10
    # Asp189 and Ser195 both line the active site; a few nanometres apart at most.
    assert 3 < mean < 25, f"{mean:.1f} A is not two residues in one active site"


def test_results_are_json_serialisable(session):
    request = analysis.parse_request(
        {"metrics": ["rmsd", "rgyr", "rmsf", "pca"], "selection": "protein"})
    json.dumps(analysis.run(session / "sess" / "results.fxa", request))


def test_reanalysis_panels_reuse_the_main_builders(session):
    """A re-analysed RMSD should look identical to one from the bundle."""
    request = analysis.parse_request(
        {"metrics": ["rmsd", "rgyr"], "selection": "protein"})
    metrics = analysis.run(session / "sess" / "results.fxa", request)
    panels = plots.build_reanalysis(metrics)
    assert {p["id"] for p in panels} == {"re-rmsd", "re-rgyr"}
    for panel in panels:
        assert panel["figure"]["layout"]["title"]["text"] == ""


# ---------------------------------------------------------------------------
#  The worker entry point
# ---------------------------------------------------------------------------


def test_the_worker_writes_its_result_and_frees_the_lock(session):
    request_path = session / "sess" / "request.json"
    output_path = session / "sess" / "result.json"
    request_path.write_text(json.dumps({"metrics": ["rgyr"], "selection": "protein"}))

    analysis.acquire_lock(session)
    assert analysis.run_to_file(session / "sess" / "results.fxa",
                                request_path, output_path) == 0

    payload = json.loads(output_path.read_text())
    assert payload["status"] == "ready"
    assert payload["metrics"]["rgyr_nm"]
    assert not (session / ".reanalyse.lock").exists(), "the slot was not released"


def test_the_worker_reports_failure_instead_of_dying(session):
    request_path = session / "sess" / "request.json"
    output_path = session / "sess" / "result.json"
    request_path.write_text(json.dumps({"metrics": ["rmsd"], "selection": "resname ZZZ"}))

    analysis.acquire_lock(session)
    assert analysis.run_to_file(session / "sess" / "results.fxa",
                                request_path, output_path) == 1

    payload = json.loads(output_path.read_text())
    assert payload["status"] == "error"
    assert "matched no atoms" in payload["message"]
    assert not (session / ".reanalyse.lock").exists(), \
        "a failed job must still free the slot"


def test_the_worker_survives_a_corrupt_request(session):
    request_path = session / "sess" / "request.json"
    output_path = session / "sess" / "result.json"
    request_path.write_text("{not json")

    analysis.acquire_lock(session)
    assert analysis.run_to_file(session / "sess" / "results.fxa",
                                request_path, output_path) == 1
    assert json.loads(output_path.read_text())["status"] == "error"
    assert not (session / ".reanalyse.lock").exists()


# ---------------------------------------------------------------------------
#  Binding-site panels, from a version 1 upload with a ligand in it
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def measured():
    """The committed trypsin and benzamidine run, measured as a complex."""
    from flexappeal import fxa, package

    fixture = FIXTURES / "trypsin_ben.fxa"
    if not fixture.is_file():
        pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")
    pytest.importorskip("mdtraj")
    result, _ = package.measure_results(fxa.load(fixture.read_bytes()), "trypsin_ben")
    if result is None:
        pytest.skip("the fixture could not be measured")
    return result


def test_an_apo_run_has_no_binding_panels_and_says_why():
    """A run with no ligand has no binding site, and that is not an error."""
    from flexappeal import fxa, package

    fixture = FIXTURES / "lysozyme.fxa"
    if not fixture.is_file():
        pytest.skip("tests/fixtures/lysozyme.fxa is missing")
    result, reasons = package.measure_results(fxa.load(fixture.read_bytes()))
    assert result is None
    assert any("no ligand" in reason for reason in reasons)


def test_every_binding_panel_builds_for_a_real_complex(measured):
    panels = plots.build_binding(measured)
    ids = {panel["id"] for panel in panels}
    # The six sub-score inputs, the matrices and the ligand's own chemistry.
    assert {"binding-pose", "binding-centroid", "binding-occupancy",
            "binding-interactions", "binding-barcode", "binding-clusters",
            "binding-strain", "binding-torsions", "binding-ligand_rmsf",
            "binding-pocket_rmsf", "binding-rmsd_2d", "binding-dccm"} <= ids


def test_every_binding_panel_is_json_serialisable(measured):
    """The template serialises each figure into the page.

    Every series here came out of numpy, and a stray float32 reaches the
    template as an unserialisable object rather than as a number: the panel then
    fails to render with a traceback about the whole page.
    """
    json.dumps(plots.build_binding(measured))


def test_every_binding_panel_is_titled_and_explained(measured):
    for panel in plots.build_binding(measured):
        assert panel["title"]
        assert len(panel["blurb"]) > 40, panel["id"]
        # The card's own heading names the panel, so an in-figure title would
        # say the same thing twice.
        assert panel["figure"]["layout"]["title"] == {"text": ""}


def test_the_page_weight_stays_reasonable(measured):
    """A 500 by 500 matrix is a quarter of a million numbers in the HTML.

    The file keeps them, which is right for a stored payload; the page cannot
    show a quarter of them, so both matrices are block-reduced on the way out.
    """
    blob = json.dumps(plots.build_binding(measured))
    assert len(blob) < 400_000, f"{len(blob) // 1024} KB of panel JSON"


def test_a_long_run_does_not_put_its_every_frame_on_the_page():
    """A page is a few hundred pixels wide; a 100 ns run has ten thousand frames.

    Undecimated, the 100 ns example came to 24 MB on one page, 18.7 of it the
    DSSP heat map alone, which carries a cell and a state name per residue per
    frame. The file still holds every frame -- this is about the wire, and
    `/reanalyse` works from the full arrays.
    """
    import json

    import numpy as np

    from flexappeal import plots

    frames = 10_000
    rng = np.random.default_rng(0)
    metrics = {
        "time_ns": [i * 0.01 for i in range(frames)],
        "rmsd_nm": list(rng.normal(0.2, 0.01, frames)),
        "rgyr_nm": list(rng.normal(1.4, 0.01, frames)),
        "sasa_total_nm2": list(rng.normal(70.0, 1.0, frames)),
        "dssp_helix_fraction": list(rng.normal(0.4, 0.01, frames)),
        "dssp_sheet_fraction": list(rng.normal(0.1, 0.01, frames)),
        "pca_projection": [[float(a), float(b)] for a, b in rng.normal(0, 1, (frames, 2))],
        "pca_variance_ratio": [0.4, 0.2],
    }
    arrays = {"dssp": rng.integers(0, 4, (frames, 162))}

    builders = {
        "rmsd": plots.rmsd(metrics),
        "rgyr": plots.radius_of_gyration(metrics),
        "sasa": plots.sasa(metrics),
        "ss_fractions": plots.secondary_structure_fractions(metrics),
        "pca": plots.pca(metrics),
        "dssp": plots.dssp(metrics, arrays),
    }
    for name, figure in builders.items():
        assert figure is not None, name
        size = len(json.dumps(figure))
        assert size < 1_000_000, f"{name} is {size // 1024} KB of JSON"

    for name in ("rmsd", "rgyr", "sasa", "ss_fractions", "pca"):
        trace = builders[name]["data"][0]
        assert len(trace["x"]) <= plots.WIRE_SERIES_POINTS, (name, len(trace["x"]))

    heatmap = builders["dssp"]["data"][0]
    assert len(heatmap["x"]) <= plots.WIRE_HEATMAP_COLUMNS, len(heatmap["x"])
    # The residue axis is the one a reader reads; it is never decimated.
    assert len(heatmap["z"]) == 162


def test_a_decimated_series_keeps_its_axis_the_same_length():
    """A strided series against an unstrided axis draws the wrong curve."""
    import numpy as np

    from flexappeal import plots

    frames = 5_000
    metrics = {
        "time_ns": [i * 0.02 for i in range(frames)],
        "rmsd_nm": list(np.linspace(0.1, 0.3, frames)),
        "dssp_helix_fraction": list(np.linspace(0.3, 0.5, frames)),
        "dssp_sheet_fraction": list(np.linspace(0.2, 0.1, frames)),
    }
    for figure in (plots.rmsd(metrics),
                   plots.secondary_structure_fractions(metrics)):
        for trace in figure["data"]:
            assert len(trace["x"]) == len(trace["y"]), trace.get("name")

    # The decimation keeps the ends, so the axis still spans the whole run.
    axis = plots.rmsd(metrics)["data"][0]["x"]
    assert axis[0] == metrics["time_ns"][0]
    assert axis[-1] >= metrics["time_ns"][-1] - 0.2


# A real difference block from the benzene run, which is what the panel draws.
_GB_TERMS = {
    "BOND": [-0.0, 0.0, 0.0], "ANGLE": [-0.0, 0.0, 0.0], "DIHED": [0.0, 0.0, 0.0],
    "VDWAALS": [-11.2868, 7.6385, 0.7639], "EEL": [-1.1099, 0.8477, 0.0848],
    "1-4 VDW": [0.0, 0.0001, 0.0], "1-4 EEL": [-0.0, 0.0001, 0.0],
    "EGB": [4.3853, 2.9658, 0.2966], "ESURF": [-1.7661, 1.1847, 0.1185],
    "DELTA G gas": [-12.3967, 8.3824, 0.8382],
    "DELTA G solv": [2.6192, 1.7972, 0.1797],
    "DELTA TOTAL": [-9.7775, 6.643, 0.6643],
}


def test_the_energy_terms_panel_is_absent_without_an_energy(measured):
    """Most runs have no MM/GBSA, and a panel of nothing is worse than no panel."""
    from flexappeal import plots

    measured.metrics.pop("mmpbsa_gb_terms", None)
    assert plots.energy_terms(measured) is None
    ids = {p["id"] for p in plots.build_binding(measured)}
    assert "binding-energy_terms" not in ids


def test_the_energy_terms_panel_draws_the_difference_block(measured):
    """The last panel DeepFlex had and this did not."""
    from flexappeal import plots

    measured.metrics["mmpbsa_gb_terms"] = dict(_GB_TERMS)
    trace = plots.energy_terms(measured)["data"][0]

    assert trace["y"] == list(plots.ENERGY_TERM_ORDER)
    assert trace["x"][0] == pytest.approx(-11.2868)
    # The spread, not the standard error: a standard error shrinks with frame
    # count and would hide how variable a term is.
    assert trace["error_x"]["array"][0] == pytest.approx(7.6385)
    assert trace["orientation"] == "h"

    ids = {p["id"] for p in plots.build_binding(measured)}
    assert "binding-energy_terms" in ids


def test_the_bonded_terms_are_not_drawn(measured):
    """They cancel to zero by construction in a single trajectory, so five bars
    at zero would say nothing. That they cancel is a check on the calculation,
    which `energy._check_cancellation` makes, not a finding to plot."""
    from flexappeal import plots

    measured.metrics["mmpbsa_gb_terms"] = dict(_GB_TERMS)
    drawn = plots.energy_terms(measured)["data"][0]["y"]
    for term in ("BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"):
        assert term not in drawn, term
    # And the total, which is printed beside the chart rather than competing
    # with its own parts for width. DeepFlex's panel filters it the same way.
    assert "DELTA TOTAL" not in drawn


def test_a_favourable_term_is_coloured_differently_from_a_costly_one(measured):
    """The sign is the first thing anybody asks of an energy term."""
    from flexappeal import plots

    measured.metrics["mmpbsa_gb_terms"] = dict(_GB_TERMS)
    trace = plots.energy_terms(measured)["data"][0]
    colours = dict(zip(trace["y"], trace["marker"]["color"]))
    # VDWAALS is favourable here, EGB is the desolvation penalty.
    assert colours["VDWAALS"] != colours["EGB"]


def test_the_matrices_are_reduced_for_the_wire(measured):
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    for key in ("binding-rmsd_2d", "binding-dccm"):
        z = panels[key]["figure"]["data"][0]["z"]
        assert len(z) <= plots.WIRE_MATRIX_SIDE, (key, len(z))
        assert len(z[0]) == len(z), "a reduced matrix must stay square"


def test_the_matrix_axes_stay_the_same_length_as_the_matrix(measured):
    """A decimated matrix with an undecimated axis silently mislabels every cell."""
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    for key in ("binding-rmsd_2d", "binding-dccm"):
        trace = panels[key]["figure"]["data"][0]
        assert len(trace["x"]) == len(trace["z"]), key
        assert len(trace["y"]) == len(trace["z"]), key


def test_the_pose_panel_shows_both_fits(measured):
    """Fitted on the pocket and fitted on the ligand: the pair is the point."""
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    names = {trace["name"] for trace in panels["binding-pose"]["figure"]["data"]}
    assert names == {"pocket fit", "ligand fit"}


def test_the_discarded_head_is_shaded_on_the_scored_panels(measured):
    """A chart that includes the equilibration beside a number that does not."""
    if measured.metrics["equil_start_frame"] == 0:
        pytest.skip("nothing was discarded in this run")
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    shapes = panels["binding-pose"]["figure"]["layout"]["shapes"]
    assert any(s["type"] == "rect" for s in shapes)


def test_the_occupancy_panel_draws_the_threshold(measured):
    """The line is what `in_pocket_fraction` counts, so it has to be visible."""
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    shapes = panels["binding-occupancy"]["figure"]["layout"]["shapes"]
    assert any(s["type"] == "line" for s in shapes)


def test_the_interaction_panel_is_ordered_by_occupancy(measured):
    """The question is which contacts carried the pose, not where they sit."""
    panels = {p["id"]: p for p in plots.build_binding(measured)}
    values = panels["binding-interactions"]["figure"]["data"][0]["x"]
    assert values == sorted(values)
    assert all(0.0 <= v <= 1.0 for v in values)


def test_the_panels_survive_a_server_without_scipy_or_rdkit():
    """The droplet installs a deliberately thin set of packages.

    `requirements.txt` carries SciPy for the pose clusters but not RDKit, which
    is an order of magnitude larger and only two panels need. Either missing has
    to cost those panels and not the page: an unguarded import here would 500
    the Analysis tab for every ligand run on the deployed site.
    """
    import sys

    from flexappeal import fxa, package

    fixture = FIXTURES / "trypsin_ben.fxa"
    if not fixture.is_file():
        pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")

    blocked = ["scipy", "scipy.cluster", "scipy.cluster.hierarchy",
               "scipy.spatial", "scipy.spatial.distance",
               "rdkit", "rdkit.Chem", "rdkit.Chem.AllChem", "rdkit.Geometry",
               "rdkit.Chem.rdMolTransforms", "rdkit.Chem.rdDetermineBonds",
               "rdkit.Chem.Draw", "rdkit.RDLogger"]
    saved = {name: sys.modules.get(name) for name in blocked}
    try:
        # A None in sys.modules makes `import name` raise ImportError, which is
        # exactly what a server without the package does.
        for name in blocked:
            sys.modules[name] = None
        measured, _ = package.measure_results(fxa.load(fixture.read_bytes()))
        assert measured is not None and measured.status == "ok"
        ids = {panel["id"] for panel in plots.build_binding(measured)}
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module

    # The nine that need neither still draw.
    assert {"binding-pose", "binding-centroid", "binding-occupancy",
            "binding-interactions", "binding-barcode", "binding-ligand_rmsf",
            "binding-pocket_rmsf", "binding-rmsd_2d", "binding-dccm"} <= ids
    # The three that do are absent rather than broken.
    assert not ids & {"binding-clusters", "binding-strain", "binding-torsions"}


def test_a_measurement_that_raises_costs_the_panels_not_the_page(monkeypatch):
    """A traceback here would take down a page that otherwise works."""
    from flexappeal import package, webapp

    def explode(*args, **kwargs):
        raise RuntimeError("measurement blew up")

    monkeypatch.setattr(package, "measure_results", explode)
    webapp._BINDING_CACHE.clear()

    application = webapp.create_app()
    application.config.update(TESTING=True)
    with application.app_context():
        from flexappeal import fxa

        fixture = FIXTURES / "trypsin_ben.fxa"
        if not fixture.is_file():
            pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")
        content = fixture.read_bytes()
        record = webapp._binding_panels(fxa.load(content), content, "x")
    assert record["panels"] == []
    assert record["depiction"] is None
    assert any("could not be built" in note for note in record["notes"])
    webapp._BINDING_CACHE.clear()
