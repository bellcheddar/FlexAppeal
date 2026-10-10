"""The DeepFlex binding grade, computed for one complex.

This is a *copy* of a judgement that lives somewhere else, which is the thing to
be careful about. The weights, the ramps, the bands and the two hard gates below
are DeepFlex's `config/grade.toml` and `pipeline/grade.py`, transcribed, and the
transcription is only worth having while it agrees with them. `tests/test_grade.py`
pins every number against the values quoted in `DEEPFLEX_FXA2_SPEC.md` section 7
so a drift fails here rather than showing a reader two different grades for one
run.

A single complex is always graded, which section 0 of the spec now requires
rather than merely permits. The arithmetic makes it possible: five of the six
sub-scores are ramps between two absolute thresholds and compare a complex with
nothing but those thresholds, so 85% of the weight and the letter bands
themselves need no cohort at all.

Only `S_E` does, being a percentile among the other ligands run against the same
receptor. For a cohort of one there is no percentile to take, so a measured
binding energy is scored against the absolute thresholds in `ENERGY_SINGLE_RAMP`
instead: DeepFlex's own fallback, not a second opinion invented here. A complex
scored that way says so on its row, because it is the weaker of the two
instruments. Where the file carries no binding energy at all, the sub-score is
recorded pending and its 15% renormalised across the other five, which is what
DeepFlex does with any sub-score a file does not carry.

What it is never scored as is zero. That was a real defect on the DeepFlex side,
found on 2026-10-10 and fixed, where a lone complex was handed the 0th percentile
and quietly lost up to 15 points for having nobody to be compared with.

What a single complex still does not get is a league table, a rank in its
receptor, or a place in a grade distribution. Those are properties of a package,
and a package of one has none of them.

Nothing in this module may be written into a `.fxa`. Spec section 7 reserves
`grade`, `score`, `s_*`, `gate_flag`, `gate_applied` and `pending_subscores` to
DeepFlex on ingest, precisely so a second copy cannot drift from the first; this
one exists to be rendered and then thrown away.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# ===========================================================================
#  DeepFlex's config/grade.toml, transcribed
# ===========================================================================

WEIGHTS = {
    "pose": 0.30,
    "occupancy": 0.20,
    "interactions": 0.20,
    "energy": 0.15,
    "convergence": 0.10,
    "conformation": 0.05,
}

# Composite score to letter, with sub-bands within `EDGE_MARGIN` of a boundary.
A_MIN = 75.0
B_MIN = 50.0
EDGE_MARGIN = 5.0

# Each sub-score is a linear ramp between full marks and zero. The direction is
# carried by which threshold is larger, so a lower-is-better quantity needs no
# special case.
RAMPS = {
    # Median ligand RMSD against frame 0, pocket-Calpha fit, in angstrom.
    "pose": (1.5, 4.0),
    # Fraction of scored frames holding at least 25% of the frame-0 contacts.
    "occupancy": (0.95, 0.50),
    # Mean of key-interaction occupancy and the Tanimoto to frame 0.
    "interactions": (0.80, 0.30),
    # Mean pairwise Jensen-Shannon divergence across four blocks of the window.
    "convergence": (0.05, 0.30),
    # MMFF94 strain in excess of the frame-0 pose, kcal/mol. Full marks at or
    # below zero: a ligand that relaxes out of a strained starting pose has done
    # something good, not something neutral.
    "strain": (0.0, 20.0),
    # Rotamer changes per nanosecond per rotatable torsion.
    "flips": (1.0, 4.0),
}

# Scoring the energy where there is no cohort to take a percentile in, from
# DeepFlex's `[energy]` block. Used only there: wherever a real cohort exists the
# percentile is the better instrument and wins.
#
# Zero is physics -- no binding free energy, no binding -- so that end needs no
# calibration. -60 kcal/mol is empirical: it is the strong end of DeepFlex's own
# reference distribution, taken off that distribution rather than off its single
# strongest value so one outlier cannot define the scale. The distribution itself
# is not reproduced here, and this repository is public, which is why.
#
# What this scale cannot do is recorded in DeepFlex's Q007 and is not small.
# There is no external anchor behind -60 beyond the spread of that data, and the
# number is pocket-dependent: a mean binding energy is a property of a site as
# much as of a ligand, so the same compound scores differently in two pockets. A
# complex scored this way says so on its row. A package is better evidence
# wherever there is one.
ENERGY_SINGLE_RAMP = (-60.0, 0.0)

# Below this many typed interactions at frame 0, the fingerprint is too coarse to
# trust on its own. The fallback to key-interaction occupancy alone is a caution
# and must never raise a score, so it applies only where it is the lower reading.
MIN_T0_INTERACTIONS = 4

# The hard gates, which override the composite rather than feeding into it.
GATE_DRIFT_A = 5.0
GATE_DRIFT_FRACTION = 0.20
DRIFTED_CAPS_AT = "C"
ALTERNATIVE_POSE_CAPS_AT = "B"

_BAND_ORDER = ["C", "C+", "B-", "B", "B+", "A-", "A"]

# What each sub-score is called on the page, and the measurement it reads.
SUBSCORE_LABELS = {
    "pose": "Pose",
    "occupancy": "Occupancy",
    "interactions": "Interactions",
    "energy": "Energy",
    "convergence": "Convergence",
    "conformation": "Conformation",
}


def ramp(value: float | None, full_marks_at: float, zero_at: float) -> float | None:
    """Linear 0..100 between two thresholds, in either direction.

    None propagates, because an absent measurement has to stay absent all the
    way to the weighting rather than becoming a zero on the way.
    """
    if value is None:
        return None
    if full_marks_at == zero_at:
        return 50.0
    fraction = (float(value) - zero_at) / (full_marks_at - zero_at)
    return float(min(max(fraction, 0.0), 1.0) * 100.0)


def letter(score: float) -> str:
    """The band, with a sub-band where the score sits near a boundary."""
    if score >= A_MIN:
        return "A" if score >= A_MIN + EDGE_MARGIN else "A-"
    if score >= B_MIN:
        if score >= A_MIN - EDGE_MARGIN:
            return "B+"
        return "B" if score >= B_MIN + EDGE_MARGIN else "B-"
    return "C+" if score >= B_MIN - EDGE_MARGIN else "C"


def cap_to(grade: str, ceiling: str) -> str:
    """Hold a grade at a ceiling, leaving anything already below it alone."""
    if _BAND_ORDER.index(grade) <= _BAND_ORDER.index(ceiling):
        return grade
    return ceiling


@dataclass
class Grade:
    """One complex's grade, and enough of its working to be auditable."""

    score: float | None = None
    grade: str | None = None
    gate_flag: str | None = None
    gate_applied: str | None = None

    subscores: dict[str, float | None] = field(default_factory=dict)
    weights_used: dict[str, float] = field(default_factory=dict)
    pending: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    # The measured binding energy, where the run computed one. Carried separately
    # from the sub-scores because it is a measurement that exists and a score
    # that does not: see the note in `score_one`.
    energy_dg: float | None = None
    energy_sd: float | None = None
    # Whether the energy sub-score came from the absolute thresholds rather than
    # a percentile. Said on the row, so a reader knows which instrument it was.
    energy_absolute: bool = False
    # The rotatable-bond count, where the run reported one. Zero is a real
    # measurement and is what explains a pending conformation sub-score.
    n_rotatable_bonds: int | None = None

    # The inputs behind each sub-score, so the meter can show what it read.
    inputs: dict[str, Any] = field(default_factory=dict)

    @property
    def band(self) -> str:
        """The letter without its sub-band, for styling."""
        return (self.grade or "")[:1]


def _drift_fraction(measured) -> float | None:
    """Share of scored frames with the ligand centroid over the gate distance.

    The gate is a fraction of frames rather than a mean, so it needs the series
    and not the summary: a ligand that sits 1 A out for most of a run and 20 A
    out for a tenth of it has a mean under the threshold and is exactly what the
    gate exists to catch.
    """
    series = getattr(measured, "series", None) or {}
    drift = series.get("com_drift_a")
    if drift is None or len(drift) == 0:
        return None
    start = int(measured.metrics.get("equil_start_frame") or 0)
    scored = drift[start:] if start < len(drift) else drift
    if len(scored) == 0:
        return None
    return float((scored > GATE_DRIFT_A).mean())


def score_one(measured) -> Grade | None:
    """Grade a single complex, or None if it was never measured.

    Every sub-score the run carries inputs for is computed, and any whose inputs
    are missing is left out and its weight renormalised across the rest, which is
    what DeepFlex does on ingest. For one complex that is usually the energy (no
    cohort, so the absolute ramp) and, for a ligand with no rotatable bond, the
    conformation (no torsion, so no strain to measure).
    """
    if measured is None or measured.status != "ok":
        return None

    metrics = measured.metrics or {}
    inter = measured.interactions or {}
    result = Grade()

    # --- S_pose, S_occ, S_conv: one measurement each --------------------------
    pose_input = metrics.get("ligand_rmsd_median")
    occ_input = metrics.get("in_pocket_fraction")
    conv_input = metrics.get("conv_js_mean")

    result.subscores["pose"] = ramp(pose_input, *RAMPS["pose"])
    result.subscores["occupancy"] = ramp(occ_input, *RAMPS["occupancy"])
    result.subscores["convergence"] = ramp(conv_input, *RAMPS["convergence"])

    # --- S_int: the mean of two readings of the same fingerprint --------------
    occupancy = inter.get("key_interaction_occupancy_mean")
    tanimoto = inter.get("tanimoto_t0_mean")
    n_t0 = inter.get("n_t0_interactions")
    parts = [v for v in (occupancy, tanimoto) if v is not None]
    s_int = ramp(sum(parts) / len(parts), *RAMPS["interactions"]) if parts else None

    if (n_t0 or 0) < MIN_T0_INTERACTIONS and occupancy is not None:
        alone = ramp(occupancy, *RAMPS["interactions"])
        if s_int is None or alone < s_int:
            result.notes.append(
                f"only {n_t0} typed interaction(s) at the first frame, below the "
                f"floor of {MIN_T0_INTERACTIONS}: the interaction sub-score falls "
                f"back to key-interaction occupancy alone, which is the lower "
                f"reading here")
            s_int = alone
        else:
            result.notes.append(
                f"only {n_t0} typed interaction(s) at the first frame, below the "
                f"floor of {MIN_T0_INTERACTIONS}. The sparse-fingerprint fallback "
                f"would score {alone:.0f} against {s_int:.0f} on the full "
                f"measurement, so it is not applied: a caution must not raise a "
                f"grade")
    result.subscores["interactions"] = s_int

    # --- S_conf: strain and torsion stability, averaged ----------------------
    strain = ramp(metrics.get("strain_kcal"), *RAMPS["strain"])
    flips = ramp(metrics.get("flips_per_ns_per_torsion"), *RAMPS["flips"])
    conf_parts = [s for s in (strain, flips) if s is not None]
    result.subscores["conformation"] = (
        sum(conf_parts) / len(conf_parts) if conf_parts else None)

    # --- S_E: measured here, ranked only in a cohort -------------------------
    # The better instrument for this one is a percentile among the other ligands
    # run against the same receptor, and a lone complex has nothing to be a
    # percentile of. It is never scored as zero for that: a 0th percentile costs
    # a flawless lone complex up to 15 points, which is the defect spec section 0
    # records against DeepFlex and had fixed.
    #
    # The *measurement* is a different matter and is carried through where the
    # run made one, because "we did not measure this" and "this cannot be ranked
    # against one ligand" are different statements and the page used to make the
    # first when the second was true.
    #
    # So it is scored against `ENERGY_SINGLE_RAMP` instead, which is DeepFlex's
    # own fallback and not a second opinion invented here. Q007 recommended
    # leaving it pending and was overruled, on the grounds that a visible 15%
    # hole in every single-run report is worse than an imperfect number that says
    # it is imperfect. The limits are real and are stated on the row: nothing
    # anchors -60 beyond the spread of the data it came off, and the scale is
    # pocket-dependent, so the same compound scores differently in two sites.
    rotatable = metrics.get("n_rotatable_bonds")
    if isinstance(rotatable, int):
        result.n_rotatable_bonds = rotatable
    result.energy_dg = metrics.get("mmpbsa_gb_dg")
    result.energy_sd = metrics.get("mmpbsa_gb_sd")
    if result.energy_dg is not None:
        result.subscores["energy"] = ramp(result.energy_dg, *ENERGY_SINGLE_RAMP)
        result.energy_absolute = True
        result.notes.append(
            f"the energy sub-score is against absolute thresholds "
            f"({ENERGY_SINGLE_RAMP[0]:.0f} to {ENERGY_SINGLE_RAMP[1]:.0f} "
            f"kcal/mol) rather than a percentile, because one complex has no "
            f"cohort to rank in. Within a single pocket these energies do not "
            f"separate good compounds from poor ones, so a study package is "
            f"better evidence where there is one")
    else:
        result.subscores["energy"] = None

    result.inputs = {
        "pose": pose_input,
        "occupancy": occ_input,
        "interactions": (sum(parts) / len(parts)) if parts else None,
        "energy": None,
        "convergence": conv_input,
        "conformation": metrics.get("strain_kcal"),
    }

    # --- the composite -------------------------------------------------------
    used = {k: w for k, w in WEIGHTS.items() if result.subscores.get(k) is not None}
    result.pending = [k for k in WEIGHTS if result.subscores.get(k) is None]
    if not used:
        return None

    total = sum(used.values())
    result.weights_used = {k: w / total for k, w in used.items()}
    result.score = sum(result.subscores[k] * w for k, w in used.items()) / total

    # --- the gates, which override rather than contribute --------------------
    fraction = _drift_fraction(measured)
    if fraction is not None and fraction > GATE_DRIFT_FRACTION:
        result.gate_flag = "drifted"
        result.notes.append(
            f"the ligand centroid is over {GATE_DRIFT_A:.0f} A from its starting "
            f"position for {fraction * 100:.0f}% of scored frames, so the grade is "
            f"capped at {DRIFTED_CAPS_AT}")
    elif metrics.get("changed_pose_cluster"):
        result.gate_flag = "alternative pose"
        result.notes.append(
            "the dominant pose cluster of the second half differs from the first "
            f"half and the new pose holds, so the grade is capped at "
            f"{ALTERNATIVE_POSE_CAPS_AT} and flagged as an alternative pose")

    base = letter(result.score)
    final = base
    if result.gate_flag == "drifted":
        final = cap_to(base, DRIFTED_CAPS_AT)
    elif result.gate_flag == "alternative pose":
        final = cap_to(base, ALTERNATIVE_POSE_CAPS_AT)
    if final != base:
        result.gate_applied = f"{base} to {final}"
    result.grade = final
    return result


def meter(graded: Grade | None) -> dict[str, Any] | None:
    """The grade shaped for the template: a dial and six rows under it."""
    if graded is None or graded.score is None:
        return None

    rows = []
    metrics_rotatable = graded.n_rotatable_bonds
    for key, weight in WEIGHTS.items():
        value = graded.subscores.get(key)
        row = {
            "key": key,
            "label": SUBSCORE_LABELS[key],
            "score": None if value is None else round(value, 1),
            "weight": weight,
            "applied": round(graded.weights_used.get(key, 0.0) * 100, 1),
            "pending": value is None,
            # The band this sub-score would fall in on its own, from the same
            # `letter` the composite uses, so a coloured bar and the dial agree
            # about where 75 and 50 are. A sub-score is not itself graded, but
            # it is a 0-100 reading on the same scale, so the thresholds
            # transfer and nothing new has to be invented to colour it.
            "band": None if value is None else letter(value)[:1],
            # What to put in the cell of a row with no score. Default says the
            # measurement is missing; the energy row overrides it when there is a
            # number but no cohort to rank it in.
            "reason": "not measured here",
            "measured": None,
        }
        if key == "conformation" and value is None:
            # Say which kind of gap it is. A rigid ligand has no torsional
            # strain to measure and no torsion to flip, which is a fact about
            # the molecule; reporting it as zero strain would be full marks.
            row["reason"] = (
                "the ligand has no rotatable bond"
                if metrics_rotatable == 0 else
                "no strain measurement in this file")
        if key == "energy" and graded.energy_dg is not None:
            # Measured but not scored, which is a different statement from not
            # measured and the one the page used to get wrong.
            row["measured"] = (
                f"{graded.energy_dg:+.1f} kcal/mol"
                + (f" ± {graded.energy_sd:.1f}" if graded.energy_sd else ""))
            row["reason"] = "measured, but ranked only within a cohort"
        elif key == "energy":
            row["reason"] = "no binding energy in this file"
        rows.append(row)

    return {
        "score": round(graded.score, 1),
        "grade": graded.grade,
        "band": graded.band,
        "rows": rows,
        "pending": [SUBSCORE_LABELS[k] for k in graded.pending],
        "gate_flag": graded.gate_flag,
        "gate_applied": graded.gate_applied,
        "notes": graded.notes,
        # Where the bands sit, so the scale is drawn from the same numbers the
        # letter came from rather than from hardcoded CSS percentages.
        "a_min": A_MIN,
        "b_min": B_MIN,
    }
