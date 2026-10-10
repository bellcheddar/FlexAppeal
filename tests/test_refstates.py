"""Tests for the reference-states analysis.

Every fixture here is built in the test. Nothing is downloaded and nothing is
committed, which is the rule the rest of this suite follows: a test that needs
the network is a test that fails when the network does, and a committed
structure is a structure somebody has to check before it is published.

The cases chosen are the ones where this build actually got the answer wrong.
Each bug in FINDINGS section 3 produced plausible numbers rather than an error,
so what is being protected is not "does it run" but "does it still refuse the
specific wrong answer it used to give".
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

from flexappeal.refstates import ensembles, membrane, sse, state, store, traj


# ===========================================================================
#  The state score
# ===========================================================================


def test_the_same_word_means_opposite_states_in_different_metrics():
    """`broken` is inactive for the R-spine and active for the ionic lock.

    A single table keyed on the word alone labelled every broken lock inactive
    and then agreed with the geometry for the wrong reason.
    """
    assert state.rule_state("r_spine", "broken") == "inactive"
    assert state.rule_state("ionic_lock", "broken") == "active"


def test_an_unknown_metric_or_label_does_not_guess():
    assert state.rule_state("not_a_metric", "broken") is None
    assert state.rule_state("r_spine", "sideways") is None
    assert state.rule_state("r_spine", "unresolved") is None


def test_geometry_and_rule_disagreeing_reports_conflict():
    motif = state.score_motif(
        "xDFG", active={"2F4J": 1.0}, inactive={"3PYY": 4.0}, n_residues=4,
        rule_label="broken", rule_metric="r_spine")
    assert motif.geometric_label == "active"
    assert motif.label == "conflict"


def test_geometry_and_rule_agreeing_reports_the_state():
    motif = state.score_motif(
        "TM6", active={"3SN6": 1.0}, inactive={"2RH1": 4.0}, n_residues=38,
        rule_label="broken", rule_metric="ionic_lock")
    assert motif.label == "active"


def test_a_motif_with_no_rule_falls_back_to_geometry_and_says_so():
    motif = state.score_motif("hinge", active={"a": 1.0}, inactive={"b": 4.0},
                              n_residues=3)
    assert motif.label == "active"
    assert "rule_label" not in motif.as_row()


def test_an_unmeasured_motif_reports_nothing_rather_than_zero():
    motif = state.score_motif("absent", active={}, inactive={}, n_residues=0)
    assert motif.s is None
    row = motif.as_row()
    for key in ("s", "d_act_a", "d_inact_a", "label"):
        assert key not in row, f"{key} was written for an unmeasured motif"


def test_the_whole_score_is_weighted_by_residue_count():
    """A twenty-residue loop is more of the protein than one tryptophan."""
    frame = state.FrameState(time_ns=1.0, motifs=[
        state.score_motif("big", active={"a": 1.0}, inactive={"b": 3.0},
                          n_residues=20),
        state.score_motif("small", active={"a": 3.0}, inactive={"b": 1.0},
                          n_residues=1),
    ])
    # Unweighted these cancel exactly; weighted the large motif decides.
    assert frame.whole_s is not None
    assert frame.whole_s > 0.3
    assert frame.whole_label == "active"


def test_smoothing_leaves_gaps_as_gaps():
    """An interpolated point on a state trace reads as a measurement."""
    values = [0.1, 0.2, None, 0.4, 0.5]
    out = state.smooth(values, frame_ps=100.0, window_ns=0.2)
    assert out[2] is None
    assert all(v is not None for i, v in enumerate(out) if i != 2)


def test_transitions_ignore_flickers_but_keep_real_changes():
    labels = (["active"] * 20 + ["inactive"] + ["active"] * 20
              + ["inactive"] * 20)
    times = [i * 0.1 for i in range(len(labels))]
    events = state.transitions(labels, times, min_dwell_ns=1.0)
    assert len(events) == 1, f"expected one real transition, got {events}"
    assert events[0]["to"] == "inactive"


# ===========================================================================
#  Superposition
# ===========================================================================


def test_superposition_recovers_a_known_rotation_and_shift():
    rng = np.random.RandomState(0)
    original = rng.rand(30, 3) * 10.0
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = original @ rotation.T + np.array([5.0, 1.0, 2.0])
    transform = traj.superpose(original, moved)
    assert traj.rmsd(traj.apply_transform(original, *transform), moved) < 1e-8


def test_superposition_refuses_to_fit_a_mirror_image():
    """An improper rotation would report a small RMSD for the wrong handedness."""
    rng = np.random.RandomState(1)
    original = rng.rand(30, 3) * 10.0
    mirrored = original * np.array([1.0, 1.0, -1.0])
    transform = traj.superpose(original, mirrored)
    fitted = traj.apply_transform(original, *transform)
    assert traj.rmsd(fitted, mirrored) > 1.0
    assert np.linalg.det(transform[0]) > 0


def test_references_in_unrelated_frames_are_brought_into_one():
    """The bug that read a confident 0.877 "active" on a nonsense comparison."""
    rng = np.random.RandomState(2)
    core = rng.rand(20, 3) * 10.0
    motif = rng.rand(5, 3) * 10.0

    def rotated(seed):
        angle = 0.7 * seed
        spin = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle), np.cos(angle), 0.0],
                         [0.0, 0.0, 1.0]])
        shift = np.array([seed * 30.0, 0.0, 0.0])
        return core @ spin.T + shift, motif @ spin.T + shift

    built = {}
    for which, seeds in (("active", (0, 1)), ("inactive", (2, 3))):
        found = ensembles.MotifEnsemble(state=which)
        for seed in seeds:
            core_i, motif_i = rotated(seed)
            found.by_code[f"R{seed}"] = {ensembles.CORE_KEY: core_i,
                                         "m": motif_i}
        built[which] = found

    # Before: the same motif sits tens of angstrom from itself across frames.
    spread = [traj.rmsd(built[w].by_code[c]["m"], built["active"].by_code["R0"]["m"])
              for w in built for c in built[w].by_code]
    assert max(spread) > 10.0

    report = ensembles.align_to_common_frame(built)
    assert report["n_aligned"] == 4
    after = [traj.rmsd(built[w].by_code[c]["m"], built["active"].by_code["R0"]["m"])
             for w in built for c in built[w].by_code]
    assert max(after) < 1e-6, "identical motifs should coincide once aligned"


# ===========================================================================
#  Secondary structure
# ===========================================================================


def test_an_unassigned_residue_does_not_shift_the_motif():
    """mdtraj writes "NA", two characters, and a joined string misaligns."""
    row = ["H", "H", sse.UNASSIGNED, "H", "E", "E", " "]
    motif = sse.motif_sse(row, "test", list(range(7)), list(range(100, 107)))
    assert len(motif.string) == 7
    assert motif.labels[2] == sse.UNASSIGNED


def test_unassigned_residues_are_not_counted_as_coil():
    helix, strand = sse.fractions(["H", "H", "H", sse.UNASSIGNED])
    assert helix == pytest.approx(1.0)
    assert strand == pytest.approx(0.0)


def test_sse_agreement_omits_residues_with_no_data():
    frame = ["H", "E", sse.UNASSIGNED]
    active = sse.reference_profile({"a": ["H", "H", "H"], "b": ["H", "E", "H"]})
    inactive = sse.reference_profile({"c": ["E", "E", "E"]})
    rows = sse.agreement(frame, active, inactive)
    assert [r["index"] for r in rows] == [0, 1]
    assert rows[0]["active_fraction"] == pytest.approx(1.0)
    assert rows[0]["inactive_fraction"] == pytest.approx(0.0)


# ===========================================================================
#  Membrane
# ===========================================================================


def test_leaflets_split_at_their_own_midplane_not_at_zero():
    """A bilayer drifts in z, and a split at the origin halves the area."""
    drifted = [40.0, 41.0, 42.0, 60.0, 61.0, 62.0]
    upper, lower = membrane.split_leaflets(drifted)
    assert len(upper) == 3 and len(lower) == 3


def test_area_per_lipid_verdicts_follow_the_spec_windows():
    assert membrane.MembraneChunk(0, 10, apl_a2=63.0).apl_verdict == "expected"
    assert membrane.MembraneChunk(0, 10, apl_a2=57.0).apl_verdict == "outside expected"
    assert membrane.MembraneChunk(0, 10, apl_a2=40.0).apl_verdict == "G6"


def test_gate_six_tolerates_one_bad_chunk_but_not_a_bad_run():
    good = [membrane.MembraneChunk(i, i + 10, apl_a2=64.0) for i in range(0, 40, 10)]
    assert membrane.gate_six(good + [membrane.MembraneChunk(40, 50, apl_a2=40.0)])["open"] is False
    assert membrane.gate_six([membrane.MembraneChunk(0, 10, apl_a2=40.0)])["open"] is True


def test_receptor_tilt_is_folded_into_the_first_quadrant():
    along_z = np.array([[0.0, 0.0, float(i)] for i in range(10)])
    assert membrane.receptor_tilt(along_z) == pytest.approx(0.0, abs=1e-6)
    flipped = np.array([[0.0, 0.0, -float(i)] for i in range(10)])
    assert membrane.receptor_tilt(flipped) == pytest.approx(0.0, abs=1e-6)


# ===========================================================================
#  The store
# ===========================================================================


def _analysis_payload() -> dict:
    return {
        "name": "fixture", "kind": "kinase",
        "core": {"n_residues": 43, "source": "fixture core"},
        "motif_residues": {"xDFG": [380, 381, 382, 383], "hinge": [316, 317, 318]},
        "ensembles": {"active": {"codes": ["2F4J"]},
                      "inactive": {"codes": ["3PYY"]}},
        "frames": [
            {"time_ns": 0.0, "n_conflicts": 0, "motifs": [
                {"motif": "xDFG", "d_act_a": 1.0, "d_inact_a": 4.0, "s": 0.6,
                 "label": "active"},
                {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0, "s": 0.0,
                 "label": "intermediate"}]},
            {"time_ns": 1.0, "n_conflicts": 0, "motifs": [
                {"motif": "xDFG", "d_act_a": 0.5, "d_inact_a": 4.5, "s": 0.8,
                 "label": "active"},
                {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0, "s": 0.0,
                 "label": "intermediate"}]},
        ],
        "notes": [],
    }


def test_the_schema_executes(tmp_path: Path):
    """Section 7 names a table `references`, which is a SQL keyword."""
    store.create(tmp_path / "t.db")
    with sqlite3.connect(tmp_path / "t.db") as connection:
        names = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'")}
    assert "reference_structures" in names
    assert "references" not in names


def test_ingesting_twice_does_not_double_the_rows(tmp_path: Path):
    path = tmp_path / "t.db"
    payload = _analysis_payload()
    first = store.ingest_analysis(path, payload, system_id="s", run_id="r",
                                  family="kinase", environment="water")
    store.ingest_analysis(path, payload, system_id="s", run_id="r",
                          family="kinase", environment="water")
    with sqlite3.connect(path) as connection:
        rows = connection.execute("SELECT COUNT(*) FROM motif_frames").fetchone()[0]
    assert rows == first["frames"]


def test_an_in_flight_run_keeps_its_own_status(tmp_path: Path):
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water",
                          run_fields={"status": "running"})
    with sqlite3.connect(path) as connection:
        status = connection.execute(
            "SELECT status FROM runs WHERE run_id = 'r'").fetchone()[0]
    assert status == "running"


def test_membrane_rows_require_their_run(tmp_path: Path):
    """The foreign keys are a constraint here, not decoration.

    SQLite leaves foreign keys off by default, which would have made the
    REFERENCES clauses in section 7's schema purely documentary. With them on, a
    membrane row naming a run that does not exist is refused.
    """
    path = tmp_path / "t.db"
    store.create(path)
    with pytest.raises(sqlite3.IntegrityError):
        with store.connect(path) as connection:
            store.put_membrane(connection, "no_such_run",
                               [{"start_ns": 0.0, "apl_a2": 63.0}])


def test_membrane_rows_are_replaced_rather_than_appended(tmp_path: Path):
    path = tmp_path / "t.db"
    # Through the real path, which creates the run before its membrane rows.
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    chunks = [{"start_ns": 0.0, "apl_a2": 63.0, "thickness_a": 38.0}]
    for _ in range(3):
        with store.connect(path) as connection:
            store.put_membrane(connection, "r", chunks)
    with sqlite3.connect(path) as connection:
        assert connection.execute(
            "SELECT COUNT(*) FROM membrane").fetchone()[0] == 1


def test_summarise_reports_the_final_tenth(tmp_path: Path):
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    summary = store.summarise(path, "r")
    assert summary["span_ns"] == [0.0, 1.0]
    xdfg = summary["motifs"]["s:xDFG"]
    assert xdfg["start_s"] == pytest.approx(0.6)
    assert xdfg["final_s"] == pytest.approx(0.8)


# ===========================================================================
#  The pages
# ===========================================================================


def test_the_pages_render_without_a_database(tmp_path: Path):
    """The feature is additive: a deployment without the file must not 500."""
    from flexappeal.webapp import create_app

    app = create_app()
    app.config["REFSTATES_DB"] = str(tmp_path / "absent.db")
    client = app.test_client()
    assert client.get("/refstates/").status_code == 200
    assert client.get("/refstates/run/nothing").status_code == 404


def test_the_finding_is_generated_from_the_data(tmp_path: Path):
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    with app.test_request_context():
        from flexappeal.refstates.web import finding

        text = finding("r")
    assert "xDFG" in text
    assert "towards active" in text


def test_the_run_page_shows_every_motif(tmp_path: Path):
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    client = app.test_client()
    page = client.get("/refstates/run/r").get_data(as_text=True)
    assert "xDFG" in page and "hinge" in page
    series = client.get("/refstates/run/r/frames.json").get_json()
    assert set(series["motifs"]) == {"s:xDFG", "s:hinge"}
    assert series["motifs"]["s:xDFG"]["s"] == [0.6, 0.8]
