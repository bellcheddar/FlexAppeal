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
