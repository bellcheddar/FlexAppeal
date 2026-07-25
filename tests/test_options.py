"""Integrity tests for the option registry.

These are the tests that make a single-source-of-truth registry safe to rely on.
Every one of them catches a class of mistake that would otherwise surface as a
broken form, a rejected valid config, or -- worst -- a generated script that
imports a symbol OpenMM no longer has, eight hours into someone's run.
"""

from __future__ import annotations

import pytest

from flexappeal import options as opts
from flexappeal import schema
from flexappeal.schema import PredicateError, evaluate_predicate


# ---------------------------------------------------------------------------
#  Structural integrity
# ---------------------------------------------------------------------------


def test_option_ids_are_unique():
    ids = [o.id for o in opts.OPTIONS]
    assert len(ids) == len(set(ids)), \
        f"duplicate option ids: {sorted({i for i in ids if ids.count(i) > 1})}"


def test_every_option_belongs_to_a_declared_group():
    known = {g.id for g in opts.GROUPS}
    for o in opts.OPTIONS:
        assert o.group in known, f"{o.id} is in unknown group {o.group!r}"


def test_every_group_has_at_least_one_option():
    for g in opts.GROUPS:
        assert opts.BY_GROUP[g.id], f"group {g.id!r} has no options"


def test_widgets_are_known():
    known = {"text", "number", "int", "select", "multiselect",
             "checkbox", "textarea", "file", "pdbid"}
    for o in opts.OPTIONS:
        assert o.widget in known, f"{o.id} uses unknown widget {o.widget!r}"


def test_select_options_have_choices_and_others_do_not():
    for o in opts.OPTIONS:
        if o.widget in ("select", "multiselect"):
            # `dynamic` options are populated from the parsed structure at
            # request time (chains, heteroatoms), so the registry cannot know
            # their choices -- but nothing else is allowed to be choice-less.
            assert o.choices or o.dynamic, f"{o.id} is a {o.widget} with no choices"
        elif o.widget != "checkbox":
            assert not o.choices, f"{o.id} is a {o.widget} but declares choices"


def test_only_multiselects_are_dynamic():
    for o in opts.OPTIONS:
        if o.dynamic:
            assert o.widget in ("select", "multiselect"), \
                f"{o.id} is dynamic but is a {o.widget}"


def test_choice_values_are_unique_within_an_option():
    for o in opts.OPTIONS:
        values = [c.value for c in o.choices]
        assert len(values) == len(set(values)), f"{o.id} has duplicate choice values"


def test_help_text_is_a_real_sentence():
    """Help must explain, not restate the label -- that is the whole convention."""
    for o in opts.OPTIONS:
        assert len(o.help) >= 40, f"{o.id} has help text that is too short to explain anything"
        # A leading digit is legitimate -- plenty of these open with the default
        # value itself ("310 K is physiological", "1 bar is atmospheric").
        assert o.help[0].isupper() or o.help[0].isdigit() or o.help.startswith(("ff", "pH")), \
            f"{o.id} help does not start with a capital or a digit"
        assert o.help.rstrip().endswith("."), f"{o.id} help does not end in a full stop"


def test_numeric_options_declare_bounds():
    for o in opts.OPTIONS:
        if o.widget in ("number", "int") and o.id != "random_seed":
            assert o.minimum is not None, f"{o.id} has no minimum"


def test_defaults_are_within_bounds():
    for o in opts.OPTIONS:
        if o.widget in ("number", "int"):
            if o.minimum is not None:
                assert o.default >= o.minimum, f"{o.id} default is below its own minimum"
            if o.maximum is not None:
                assert o.default <= o.maximum, f"{o.id} default is above its own maximum"


def test_select_defaults_are_offered_choices():
    for o in opts.OPTIONS:
        if o.widget == "select":
            assert o.default in o.choice_values(), \
                f"{o.id} default {o.default!r} is not one of its choices"
        elif o.widget == "multiselect" and o.choices:
            for v in o.default:
                assert v in o.choice_values() or v == "*", \
                    f"{o.id} default contains {v!r}, which is not one of its choices"


# ---------------------------------------------------------------------------
#  Predicates
# ---------------------------------------------------------------------------


def test_every_requires_predicate_parses_and_references_real_options():
    for o in opts.OPTIONS:
        if o.requires:
            evaluate_predicate(o.requires, opts.defaults())
        for c in o.choices:
            if c.requires:
                evaluate_predicate(c.requires, opts.defaults())


def test_predicate_rejects_unknown_option():
    with pytest.raises(PredicateError, match="unknown option"):
        evaluate_predicate("no_such_option == 'x'", opts.defaults())


def test_predicate_rejects_arbitrary_code():
    """Predicates are authored, not user input -- but the whitelist is the proof."""
    for hostile in ("__import__('os').system('true')", "open('/etc/passwd')", "x.y"):
        with pytest.raises(PredicateError):
            evaluate_predicate(hostile, opts.defaults())


@pytest.mark.parametrize(
    "expr,overrides,expected",
    [
        ("solvent_mode == 'explicit'", {}, True),
        ("solvent_mode == 'explicit'", {"solvent_mode": "implicit"}, False),
        ("use_membrane", {}, False),
        ("not use_membrane", {}, True),
        ("input_source in ('rcsb', 'opm')", {"input_source": "opm"}, True),
        ("input_source in ('rcsb', 'opm')", {"input_source": "upload"}, False),
        ("heat_duration > 0", {}, True),
        ("heat_duration > 0", {"heat_duration": 0.0}, False),
        ("replicates > 1", {"replicates": 3}, True),
        ("solvent_mode == 'explicit' and not use_membrane", {"use_membrane": True}, False),
    ],
)
def test_predicate_evaluation(expr, overrides, expected):
    cfg = opts.defaults() | overrides
    assert evaluate_predicate(expr, cfg) is expected


def test_default_config_is_self_consistent():
    """Every default must be reachable: no option whose default choice is hidden."""
    cfg = opts.defaults()
    for o in opts.OPTIONS:
        if not schema.is_active(o, cfg) or o.widget != "select":
            continue
        allowed = {c.value for c in schema.active_choices(o, cfg)}
        assert o.default in allowed, (
            f"{o.id} defaults to {o.default!r}, but that choice is hidden by its own "
            f"requires predicate in the default config"
        )


# ---------------------------------------------------------------------------
#  Validation
# ---------------------------------------------------------------------------


def test_defaults_validate_cleanly():
    result = schema.validate(opts.defaults())
    assert result.ok, f"the default config does not validate: {result.errors}"
    assert not result.warnings, f"the default config warns: {result.warnings}"


def test_unknown_key_warns_but_does_not_block():
    result = schema.validate(opts.defaults() | {"nonsense_option": 1})
    assert result.ok
    assert any("nonsense_option" in w.message for w in result.warnings)


def test_unknown_key_can_be_made_fatal():
    result = schema.validate(opts.defaults() | {"nonsense_option": 1}, strict_unknown=True)
    assert not result.ok


def test_out_of_range_value_is_an_error():
    result = schema.validate(opts.defaults() | {"ph": 99.0})
    assert not result.ok
    assert any(i.option_id == "ph" for i in result.errors)


def test_invalid_choice_is_an_error():
    result = schema.validate(opts.defaults() | {"protein_ff": "made_up.xml"})
    assert not result.ok


def test_hidden_options_are_not_validated():
    """A membrane field left at a silly value must not block a non-membrane run."""
    result = schema.validate(opts.defaults() | {"use_membrane": False, "lipid_type": "POPC"})
    assert result.ok


def test_form_post_treats_absent_checkbox_as_false():
    raw = {"job_name": "x", "use_hmr": "on"}
    result = schema.validate(raw, form_post=True)
    assert result.config["use_hmr"] is True
    assert result.config["minimize"] is False  # absent from the POST


def test_partial_config_keeps_checkbox_defaults():
    result = schema.validate({"job_name": "x"}, form_post=False)
    assert result.config["minimize"] is True


# ---------------------------------------------------------------------------
#  Physics rules
# ---------------------------------------------------------------------------


def test_long_timestep_without_constraints_is_rejected():
    result = schema.validate(opts.defaults() | {"constraints": "None", "use_hmr": False, "timestep": 4.0})
    assert not result.ok
    assert any(i.option_id == "timestep" for i in result.errors)


def test_hmr_without_constraints_is_rejected():
    result = schema.validate(opts.defaults() | {"constraints": "None", "timestep": 0.5})
    assert any(i.option_id == "use_hmr" for i in result.errors)


def test_padding_smaller_than_cutoff_is_rejected():
    result = schema.validate(opts.defaults() | {"padding": 0.8, "nonbonded_cutoff": 1.0})
    assert not result.ok
    assert any(i.option_id == "padding" for i in result.errors)


def test_switch_distance_beyond_cutoff_is_rejected():
    result = schema.validate(
        opts.defaults() | {"use_switching": True, "switch_distance": 1.5, "nonbonded_cutoff": 1.0}
    )
    assert not result.ok


def test_skipping_minimisation_in_explicit_solvent_is_rejected():
    result = schema.validate(opts.defaults() | {"minimize": False})
    assert not result.ok


def test_bad_job_name_is_rejected():
    for bad in ("", "../escape", "has space", "-leading"):
        result = schema.validate(opts.defaults() | {"job_name": bad})
        assert not result.ok, f"{bad!r} should not be an acceptable job name"


def test_bad_pdb_id_is_rejected():
    result = schema.validate(opts.defaults() | {"input_source": "rcsb", "pdb_id": "NOPE"})
    assert not result.ok


def test_good_pdb_id_is_accepted():
    result = schema.validate(opts.defaults() | {"input_source": "rcsb", "pdb_id": "1AKI"})
    assert result.ok, result.errors


# ---------------------------------------------------------------------------
#  Normalisation
# ---------------------------------------------------------------------------


def test_ff19sb_forces_opc_water():
    cfg = schema.normalise(opts.defaults() | {"protein_ff": "amber/protein.ff19SB.xml"})
    assert cfg["water_model"] == "opc"


def test_charmm_forces_its_own_water():
    cfg = schema.normalise(opts.defaults() | {"protein_ff": "charmm36.xml"})
    assert cfg["water_model"] == "charmm_tip3p"


def test_membrane_forces_the_membrane_barostat():
    cfg = schema.normalise(opts.defaults() | {"use_membrane": True})
    assert cfg["barostat"] == "MonteCarloMembrane"


def test_implicit_solvent_drops_the_barostat_and_pme():
    cfg = schema.normalise(opts.defaults() | {"solvent_mode": "implicit"})
    assert cfg["barostat"] == "none"
    assert cfg["nonbonded_method"] == "CutoffNonPeriodic"


def test_keeping_a_heteroatom_enables_ligand_parameters():
    cfg = schema.normalise(opts.defaults() | {"keep_heteroatoms": ["HEM"]})
    assert cfg["has_ligands"] is True


# ---------------------------------------------------------------------------
#  Derived quantities
# ---------------------------------------------------------------------------


def test_step_counts_follow_from_the_timestep():
    cfg = opts.defaults() | {"timestep": 4.0, "production_duration": 100.0}
    d = schema.derive(cfg)
    assert d["production_steps"] == 25_000_000  # 100 ns / 4 fs


def test_frame_count_follows_from_the_interval():
    cfg = opts.defaults() | {"production_duration": 100.0, "traj_interval": 10.0}
    assert schema.derive(cfg)["traj_frames"] == 10_000


def test_replicates_multiply_the_totals():
    one = schema.derive(opts.defaults() | {"replicates": 1})
    three = schema.derive(opts.defaults() | {"replicates": 3})
    assert three["total_steps"] == one["total_steps"] * 3
    assert three["traj_frames"] == one["traj_frames"] * 3


def test_size_estimate_scales_with_saved_atoms():
    base = opts.defaults() | {"_estimated_atoms": 50000, "_solute_atoms": 5000}
    everything = schema.derive(base | {"traj_selection": "all"})
    ca_only = schema.derive(base | {"traj_selection": "ca"})
    assert everything["traj_bytes"] > ca_only["traj_bytes"] * 50


def test_wall_time_is_honest_about_having_no_structure():
    est = schema.estimate_wall_time(opts.defaults())
    assert est["basis"] == "unknown"


def test_wall_time_uses_a_benchmark_when_given_one():
    cfg = opts.defaults() | {"_estimated_atoms": 25000}
    est = schema.estimate_wall_time(cfg, ns_per_day=100.0)
    assert est["basis"] == "benchmarked"
    assert est["hours"] > 0
