"""The binding grade.

`flexappeal/grade.py` is a transcription of DeepFlex's `config/grade.toml` and
`pipeline/grade.py`, and a transcription is only worth having while it agrees
with the original. So most of what follows pins the numbers themselves against
the values quoted in `DEEPFLEX_FXA2_SPEC.md` section 7, rather than against
whatever the module happens to contain: a drift then fails here, loudly, instead
of showing a reader one grade on this site and another on DeepFlex.
"""

from __future__ import annotations

import numpy as np
import pytest

from flexappeal import grade


# ===========================================================================
#  The constants, against the spec
# ===========================================================================


def test_the_weights_are_the_published_ones():
    """Spec section 7: pose 0.30, occupancy 0.20, interactions 0.20, energy 0.15,
    convergence 0.10, conformation 0.05."""
    assert grade.WEIGHTS == {
        "pose": 0.30, "occupancy": 0.20, "interactions": 0.20,
        "energy": 0.15, "convergence": 0.10, "conformation": 0.05,
    }
    assert sum(grade.WEIGHTS.values()) == pytest.approx(1.0)


def test_the_bands_are_the_published_ones():
    """A at 75 and above, B from 50, C below, sub-bands within 5 points."""
    assert (grade.A_MIN, grade.B_MIN, grade.EDGE_MARGIN) == (75.0, 50.0, 5.0)


def test_the_ramps_are_the_published_ones():
    assert grade.RAMPS["pose"] == (1.5, 4.0)
    assert grade.RAMPS["occupancy"] == (0.95, 0.50)
    assert grade.RAMPS["interactions"] == (0.80, 0.30)
    assert grade.RAMPS["convergence"] == (0.05, 0.30)
    assert grade.RAMPS["strain"] == (0.0, 20.0)
    assert grade.RAMPS["flips"] == (1.0, 4.0)


def test_the_gates_are_the_published_ones():
    assert grade.GATE_DRIFT_A == 5.0
    assert grade.GATE_DRIFT_FRACTION == 0.20
    assert grade.DRIFTED_CAPS_AT == "C"
    assert grade.ALTERNATIVE_POSE_CAPS_AT == "B"


# ===========================================================================
#  ramp
# ===========================================================================


@pytest.mark.parametrize(("value", "expected"), [
    (1.5, 100.0),      # full marks
    (4.0, 0.0),        # zero
    (2.75, 50.0),      # halfway
    (0.5, 100.0),      # better than full marks still scores 100, never 120
    (9.9, 0.0),        # worse than zero never scores negative
])
def test_a_lower_is_better_ramp(value, expected):
    assert grade.ramp(value, 1.5, 4.0) == pytest.approx(expected)


@pytest.mark.parametrize(("value", "expected"), [
    (0.95, 100.0),
    (0.50, 0.0),
    (0.725, 50.0),
    (1.0, 100.0),
])
def test_a_higher_is_better_ramp(value, expected):
    assert grade.ramp(value, 0.95, 0.50) == pytest.approx(expected)


def test_a_missing_measurement_propagates_rather_than_scoring_zero():
    """The difference between "bad" and "not measured", kept all the way to the
    weighting. A None that became a 0.0 here would be a silent bad grade."""
    assert grade.ramp(None, 1.5, 4.0) is None


def test_a_degenerate_ramp_does_not_divide_by_zero():
    assert grade.ramp(1.0, 2.0, 2.0) == 50.0


# ===========================================================================
#  letter and cap_to
# ===========================================================================


@pytest.mark.parametrize(("score", "expected"), [
    (100.0, "A"), (80.0, "A"), (79.9, "A-"), (75.0, "A-"),
    (74.9, "B+"), (70.0, "B+"), (69.9, "B"), (55.0, "B"),
    (54.9, "B-"), (50.0, "B-"), (49.9, "C+"), (45.0, "C+"),
    (44.9, "C"), (0.0, "C"),
])
def test_the_letter_for_a_score(score, expected):
    assert grade.letter(score) == expected


def test_a_cap_holds_a_grade_down_and_never_lifts_one():
    assert grade.cap_to("A", "C") == "C"
    assert grade.cap_to("A", "B") == "B"
    # Already at or below the ceiling: left exactly as it was.
    assert grade.cap_to("C", "B") == "C"
    assert grade.cap_to("B-", "B") == "B-"
    assert grade.cap_to("B", "B") == "B"


# ===========================================================================
#  score_one
# ===========================================================================


class _Measured:
    """The surface of measure.Measured that grading actually reads."""

    def __init__(self, metrics, interactions=None, series=None):
        self.metrics = metrics
        self.interactions = interactions
        self.series = series

    @property
    def status(self):
        return str(self.metrics.get("status") or "ok")


def _good(**overrides):
    metrics = {
        "status": "ok",
        "equil_start_frame": 0,
        "ligand_rmsd_median": 1.5,          # full marks
        "in_pocket_fraction": 0.95,         # full marks
        "conv_js_mean": 0.05,               # full marks
        "strain_kcal": 0.0,                 # full marks
        "flips_per_ns_per_torsion": 1.0,    # full marks
        "changed_pose_cluster": 0,
    }
    metrics.update(overrides)
    inter = {"status": "ok", "n_t0_interactions": 10,
             "key_interaction_occupancy_mean": 0.80,
             "tanimoto_t0_mean": 0.80}
    series = {"com_drift_a": np.zeros(100, dtype=np.float32)}
    return _Measured(metrics, inter, series)


def test_a_perfect_complex_scores_one_hundred():
    g = grade.score_one(_good())
    assert g.score == pytest.approx(100.0)
    assert g.grade == "A"
    assert g.gate_flag is None


def test_the_energy_weight_is_renormalised_away_not_filled():
    """The whole reason a single complex can be graded at all.

    A neutral 50 in place of the percentile would put 15% of every grade on a
    value nobody measured. DeepFlex renormalises instead, and so does this.
    """
    g = grade.score_one(_good())
    assert g.subscores["energy"] is None
    assert g.pending == ["energy"]
    assert "energy" not in g.weights_used
    assert sum(g.weights_used.values()) == pytest.approx(1.0)
    # Each surviving weight is its own share of the 0.85 that remains.
    assert g.weights_used["pose"] == pytest.approx(0.30 / 0.85)
    assert g.weights_used["conformation"] == pytest.approx(0.05 / 0.85)


def test_a_renormalised_grade_is_not_dragged_down_by_the_missing_sub_score():
    """A pending sub-score must not behave like a zero.

    With energy scored 0 the composite would be 85, not 100, and a flawless run
    would come out A- rather than A.
    """
    g = grade.score_one(_good())
    assert g.score == pytest.approx(100.0)


def test_the_composite_is_the_weighted_mean_of_what_was_measured():
    g = grade.score_one(_good(ligand_rmsd_median=2.75))     # S_pose 50
    expected = (50 * 0.30 + 100 * 0.20 + 100 * 0.20
                + 100 * 0.10 + 100 * 0.05) / 0.85
    assert g.score == pytest.approx(expected)


def test_the_interaction_sub_score_averages_both_readings():
    """Occupancy 0.80 and Tanimoto 0.30 average to 0.55, which is halfway."""
    m = _good()
    m.interactions["tanimoto_t0_mean"] = 0.30
    g = grade.score_one(m)
    assert g.subscores["interactions"] == pytest.approx(50.0)


def test_a_sparse_fingerprint_falls_back_only_when_that_is_harsher():
    """The protection must not become a reward.

    A thin frame-0 fingerprint with one well-held key interaction scored 100 on
    the fallback while the full measurement said much less. Applied
    unconditionally that promoted complexes it was written to protect against.
    """
    m = _good()
    m.interactions.update({"n_t0_interactions": 2,
                           "key_interaction_occupancy_mean": 0.80,  # alone: 100
                           "tanimoto_t0_mean": 0.30})               # mean: 50
    g = grade.score_one(m)
    assert g.subscores["interactions"] == pytest.approx(50.0), "a caution raised the score"
    assert any("not applied" in n for n in g.notes)


def test_a_sparse_fingerprint_does_pull_a_score_down():
    m = _good()
    m.interactions.update({"n_t0_interactions": 2,
                           "key_interaction_occupancy_mean": 0.30,  # alone: 0
                           "tanimoto_t0_mean": 0.80})               # mean: 60
    g = grade.score_one(m)
    assert g.subscores["interactions"] == pytest.approx(0.0)
    assert any("lower reading" in n for n in g.notes)


def test_conformation_rests_on_strain_alone_without_a_rotatable_bond():
    """A rigid ligand has no torsion to flip, so the flip term is simply absent
    and must not average in as a zero."""
    m = _good()
    del m.metrics["flips_per_ns_per_torsion"]
    g = grade.score_one(m)
    assert g.subscores["conformation"] == pytest.approx(100.0)


def test_a_ligand_that_relaxes_takes_full_marks_for_strain():
    """Negative excess strain is a good outcome, not a neutral one."""
    g = grade.score_one(_good(strain_kcal=-4.0))
    assert g.subscores["conformation"] == pytest.approx(100.0)


# --- the gates -------------------------------------------------------------


def test_sustained_centroid_drift_caps_the_grade_at_c():
    drift = np.zeros(100, dtype=np.float32)
    drift[:30] = 7.0                      # 30% of frames beyond 5 A
    m = _good()
    m.series = {"com_drift_a": drift}
    g = grade.score_one(m)
    assert g.gate_flag == "drifted"
    assert g.grade == "C"
    assert g.score == pytest.approx(100.0), "the gate caps the letter, not the score"
    assert g.gate_applied == "A to C"


def test_brief_drift_does_not_trip_the_gate():
    """The threshold is a fifth of the frames, so a tenth must pass."""
    drift = np.zeros(100, dtype=np.float32)
    drift[:10] = 7.0
    m = _good()
    m.series = {"com_drift_a": drift}
    g = grade.score_one(m)
    assert g.gate_flag is None
    assert g.grade == "A"


def test_the_drift_gate_ignores_the_discarded_equilibration():
    """Every scored number starts at `equil_start_frame`, and the gate is one of
    them. Counting the head would fail a run on frames nothing else scored."""
    drift = np.zeros(100, dtype=np.float32)
    drift[:40] = 7.0
    m = _good(equil_start_frame=40)
    m.series = {"com_drift_a": drift}
    g = grade.score_one(m)
    assert g.gate_flag is None


def test_a_changed_dominant_pose_caps_the_grade_at_b():
    g = grade.score_one(_good(changed_pose_cluster=1))
    assert g.gate_flag == "alternative pose"
    assert g.grade == "B"
    assert g.gate_applied == "A to B"


def test_a_gate_does_not_lift_a_poor_grade():
    """Capping at B must leave a C where it is."""
    m = _good(ligand_rmsd_median=4.0, in_pocket_fraction=0.50,
              conv_js_mean=0.30, strain_kcal=20.0,
              flips_per_ns_per_torsion=4.0, changed_pose_cluster=1)
    m.interactions.update({"key_interaction_occupancy_mean": 0.30,
                           "tanimoto_t0_mean": 0.30})
    g = grade.score_one(m)
    assert g.score == pytest.approx(0.0)
    assert g.grade == "C"
    assert g.gate_applied is None


# --- refusals --------------------------------------------------------------


def test_a_failed_complex_is_not_graded():
    """A missing result and a poor result are different findings."""
    m = _Measured({"status": "failed:no_ligand"})
    assert grade.score_one(m) is None


def test_a_complex_with_nothing_measurable_is_not_graded():
    m = _Measured({"status": "ok"}, None, None)
    assert grade.score_one(m) is None


# ===========================================================================
#  meter
# ===========================================================================


def test_a_measured_energy_is_scored_and_shown():
    """A measured energy scores on the absolute ramp and nothing stays pending.

    The row carries both: the score, which is what the grade uses, and the
    kcal/mol it was read off, which is the number a reader came for. The
    sub-score having a value is what gives the row a bar.
    """
    m = _good(mmpbsa_gb_dg=-9.7775, mmpbsa_gb_sd=6.643)
    g = grade.score_one(m)

    assert g.subscores["energy"] == pytest.approx(16.29, abs=0.01)
    assert g.pending == []
    assert g.weights_used["energy"] == pytest.approx(0.15), "no renormalisation"
    assert g.energy_absolute is True

    row = next(r for r in grade.meter(g)["rows"] if r["key"] == "energy")
    assert row["pending"] is False
    assert row["score"] == pytest.approx(16.3, abs=0.1)
    assert row["measured"] == "-9.8 kcal/mol ± 6.6"


def test_every_scored_row_has_a_number_and_can_draw_a_bar():
    """Including a zero. A sub-score of 0.0 is a measured nought, and a row that
    looked empty instead read as a missing panel."""
    g = grade.score_one(_good(conv_js_mean=0.30,                 # S_conv 0
                              mmpbsa_gb_dg=-9.7775))
    rows = {r["key"]: r for r in grade.meter(g)["rows"]}
    assert rows["convergence"]["score"] == pytest.approx(0.0)
    assert rows["convergence"]["pending"] is False
    for key, row in rows.items():
        assert row["pending"] is False, key
        assert row["score"] is not None, key


def test_the_energy_ramp_is_the_one_deepflex_falls_back_to():
    """A transcription, not a second opinion.

    -60 kcal/mol to 0, from DeepFlex's `[energy]` block. Zero is physics -- no
    binding free energy, no binding. -60 is the strong end of the 27 complexes
    it has measured.
    """
    assert grade.ENERGY_SINGLE_RAMP == (-60.0, 0.0)
    assert grade.ramp(-60.0, *grade.ENERGY_SINGLE_RAMP) == pytest.approx(100.0)
    assert grade.ramp(0.0, *grade.ENERGY_SINGLE_RAMP) == pytest.approx(0.0)
    # Stronger than full marks is still full marks, never more.
    assert grade.ramp(-200.0, *grade.ENERGY_SINGLE_RAMP) == pytest.approx(100.0)
    # A positive energy is no binding at all, and floors rather than going below.
    assert grade.ramp(5.0, *grade.ENERGY_SINGLE_RAMP) == pytest.approx(0.0)


def test_an_absent_energy_is_still_pending_and_renormalised():
    """Scoring a measured energy must not invent one where none was measured."""
    g = grade.score_one(_good())
    assert g.subscores["energy"] is None
    assert g.pending == ["energy"]
    assert g.weights_used["pose"] == pytest.approx(0.30 / 0.85)


def test_an_absent_energy_says_so_differently():
    row = next(r for r in grade.meter(grade.score_one(_good()))["rows"]
               if r["key"] == "energy")
    assert row["measured"] is None
    assert row["reason"] == "no binding energy in this file"


def test_the_energy_score_says_which_instrument_produced_it():
    """An absolute ramp and a percentile are different claims.

    Q007 recommended against this ramp and was overruled, so what it cannot do
    has to travel with the number: within one pocket these energies do not
    separate the reference agonists from the rest, which is why nothing anchors
    -60 beyond the spread of the data.
    """
    g = grade.score_one(_good(mmpbsa_gb_dg=-30.0))
    assert g.energy_absolute is True
    note = " ".join(g.notes)
    assert "absolute thresholds" in note
    assert "no cohort" in note or "cohort to rank in" in note
    assert "better evidence" in note


def test_the_meter_names_every_sub_score_including_the_pending_one():
    m = grade.meter(grade.score_one(_good()))
    assert [row["key"] for row in m["rows"]] == list(grade.WEIGHTS)
    energy = next(r for r in m["rows"] if r["key"] == "energy")
    assert energy["pending"] is True
    assert energy["score"] is None
    assert energy["weight"] == 0.15
    assert energy["applied"] == 0.0


def test_the_meter_carries_the_band_edges_it_was_scored_against():
    """The scale is drawn from these, so a threshold change moves the ticks
    rather than leaving them where a stylesheet put them."""
    m = grade.meter(grade.score_one(_good()))
    assert (m["a_min"], m["b_min"]) == (grade.A_MIN, grade.B_MIN)


def test_the_meter_is_none_for_an_ungraded_complex():
    assert grade.meter(None) is None


def test_the_applied_weights_on_the_meter_sum_to_a_hundred():
    m = grade.meter(grade.score_one(_good()))
    assert sum(r["applied"] for r in m["rows"]) == pytest.approx(100.0, abs=0.2)


# ===========================================================================
#  Against the committed example
# ===========================================================================


def test_the_committed_example_grades_from_its_own_measurements():
    """End to end on a real 100 ns run, so the ramps are exercised on numbers
    nobody chose to make the arithmetic come out."""
    pytest.importorskip("mdtraj")
    from flexappeal import fxa, package
    from flexappeal.webapp import EXAMPLE_FXA

    if not EXAMPLE_FXA.is_file():
        pytest.skip("the example run has not been built")

    results = fxa.load(EXAMPLE_FXA.read_bytes(), verify_checksums=False)
    measured, _ = package.measure_results(results, "example")
    if measured is None:
        pytest.skip("the example has no ligand to measure")

    g = grade.score_one(measured)
    assert g is not None
    assert 0.0 <= g.score <= 100.0
    assert g.grade in ("A", "A-", "B+", "B", "B-", "C+", "C")
    # The committed example computes its own MM/GBSA, so nothing is pending.
    assert g.pending == []
    assert g.energy_dg is not None and g.energy_dg < 0.0
    # Benzene in a buried cavity: it cannot leave, so occupancy is full marks
    # and the pose term is what the grade turns on.
    assert g.subscores["occupancy"] == pytest.approx(100.0)
    assert g.subscores["pose"] is not None
