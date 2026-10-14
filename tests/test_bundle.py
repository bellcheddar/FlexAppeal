"""Bundle generation tests.

The generated scripts run on someone else's machine, hours after we last saw
them, with no way to reach back. So the bar here is higher than "it rendered":
every configuration shape must produce Python that compiles, and the physics
choices that differ between shapes must actually appear in the output.

`test_every_shape_compiles` is the load-bearing one. A Jinja conditional that
emits a dangling comma or an unclosed bracket is invisible until someone runs it.
"""

from __future__ import annotations

import ast
import json
import py_compile

import pytest

from flexappeal import bundle, options as opts
from flexappeal.bundle import BundleError

from conftest import FIXTURES


@pytest.fixture(scope="module")
def structure():
    return (FIXTURES / "1aki.pdb").read_bytes()


def _cfg(**overrides):
    return opts.defaults() | {"job_name": "test_job"} | overrides


def _build(structure, **overrides):
    return bundle.build(_cfg(**overrides), structure, "1aki.pdb")


# ---------------------------------------------------------------------------
#  Shape
# ---------------------------------------------------------------------------


def test_bundle_contains_everything_needed_to_run(structure):
    result = _build(structure)
    files = bundle.unpack(result.content)
    assert set(files) == {
        "README.md", "analyse.py", "config.json", "input.pdb", "pixi.toml", "run.py",
    }
    assert files["input.pdb"] == structure, "the structure must survive the round trip byte for byte"


def test_bundle_filename_follows_the_job_name(structure):
    assert _build(structure, job_name="my_run").filename == "flexappeal_my_run.command"


def test_bundle_is_a_runnable_shell_script(structure):
    content = _build(structure).content
    assert content.startswith(b"#!/usr/bin/env bash")
    assert b"\n__FLEXAPPEAL_PAYLOAD__\n" in content
    assert b"chmod +x" in content, "the instructions must mention the executable bit"


def test_payload_lines_are_wrapped(structure):
    """An unwrapped multi-megabyte line breaks editors, mail and some awks."""
    content = _build(structure).content
    payload = content.split(b"\n__FLEXAPPEAL_PAYLOAD__\n", 1)[1]
    assert payload, "no payload"
    assert all(len(line) <= 76 for line in payload.split(b"\n"))


def test_bundles_are_deterministic(structure, monkeypatch):
    """The same configuration must produce the same bytes.

    The clock is frozen for the comparison. config.json records generated_at to
    the minute, and it is inside the payload, so two builds either side of a
    minute boundary differ legitimately -- this test used to fail roughly once
    an hour for that reason alone, which is a flaky test rather than a finding.
    """
    import datetime as _dt

    class Frozen(_dt.datetime):
        @classmethod
        def now(cls, tz=None):
            return _dt.datetime(2026, 1, 1, 12, 0, tzinfo=tz)

    monkeypatch.setattr(bundle, "datetime", Frozen)
    first = _build(structure).content
    second = _build(structure).content
    assert first == second, "the same configuration produced different bytes"


def test_config_json_carries_the_full_configuration(structure):
    files = bundle.unpack(_build(structure).content)
    payload = json.loads(files["config.json"])
    assert payload["config"]["job_name"] == "test_job"
    assert len(payload["config"]) == len(opts.OPTIONS)
    assert payload["structure"]["sha256"], "the structure checksum is the provenance link"
    assert not any(k.startswith("_") for k in payload["config"])


# ---------------------------------------------------------------------------
#  Unpacking
# ---------------------------------------------------------------------------


def test_unpack_rejects_something_that_is_not_a_bundle():
    with pytest.raises(BundleError, match="payload marker"):
        bundle.unpack(b"#!/bin/bash\necho hello\n")


@pytest.mark.parametrize("fraction", [0.1, 0.25, 0.5, 0.75, 0.9, 0.99])
def test_unpack_rejects_a_truncated_payload(structure, fraction):
    """Every cut point, not one.

    This used to test a single offset and passed on where that happened to
    land. Growing the run script moved it, and the new offset produced gzip's
    EOFError -- which is not an OSError, was not caught, and reached the user
    as a traceback.
    """
    content = _build(structure).content
    with pytest.raises(BundleError, match="corrupt or truncated"):
        bundle.unpack(content[: int(len(content) * fraction)])


# ---------------------------------------------------------------------------
#  Every configuration shape must produce valid Python
# ---------------------------------------------------------------------------

SHAPES = {
    "apo_default": {},
    "implicit": {"solvent_mode": "implicit", "nonbonded_method": "CutoffNonPeriodic"},
    "vacuum": {"solvent_mode": "vacuum", "nonbonded_method": "NoCutoff", "minimize": True},
    "membrane": {"use_membrane": True, "protein_ff": "charmm36.xml",
                 "input_source": "opm", "pdb_id": "1AKI"},
    "ligands": {"keep_heteroatoms": ["HEM"], "has_ligands": True},
    "ff19sb": {"protein_ff": "amber/protein.ff19SB.xml"},
    "charmm": {"protein_ff": "charmm36.xml", "nonbonded_cutoff": 1.2,
               "use_switching": True, "switch_distance": 1.0,
               "dispersion_correction": False},
    "no_hmr": {"use_hmr": False, "timestep": 2.0},
    "no_constraints": {"constraints": "None", "use_hmr": False, "timestep": 0.5},
    "verlet_thermostat": {"integrator": "Verlet", "use_thermostat": True},
    "nose_hoover": {"integrator": "NoseHoover"},
    "variable_timestep": {"integrator": "VariableLangevin"},
    "nvt": {"barostat": "none"},
    "anisotropic": {"barostat": "MonteCarloAnisotropic",
                    "anisotropic_pressure": "1.0, 1.0, 1.5"},
    "flexible_barostat": {"barostat": "MonteCarloFlexible"},
    "no_restraints": {"use_positional_restraints": False},
    "dcd": {"traj_format": "dcd"},
    "hdf5": {"traj_format": "hdf5"},
    "custom_selection": {"traj_selection": "custom",
                         "traj_custom_selection": "protein or resname LIG"},
    "replicates": {"replicates": 3},
    "explicit_box_vectors": {"box_sizing": "vectors", "box_vectors": "6.0, 6.0, 6.0"},
    "fixed_water_count": {"box_sizing": "num_waters", "num_waters": 5000},
    "mutations": {"mutations": "A:ALA-57-GLY\nA:VAL-2-ILE"},
    "no_repair": {"fix_missing_residues": False, "replace_nonstandard": False,
                  "add_missing_atoms": False, "strip_hydrogens": False},
    "named_platform": {"platform": "CPU", "cpu_threads": 8},
    "minimal_output": {"save_system_xml": False, "save_solvated_pdb": False,
                       "enforce_periodic_box": False},
    "no_memory_guard": {"memory_guard": False},
    "seeded": {"random_seed": 42},
}


@pytest.mark.parametrize("name,overrides", sorted(SHAPES.items()))
def test_every_shape_compiles(structure, tmp_path, name, overrides):
    result = bundle.build(_cfg(**overrides), structure, "1aki.pdb")
    files = bundle.unpack(result.content)
    for script in ("run.py", "analyse.py"):
        path = tmp_path / f"{name}_{script}"
        path.write_bytes(files[script])
        try:
            py_compile.compile(str(path), doraise=True)
        except py_compile.PyCompileError as exc:
            pytest.fail(f"{name} produced invalid {script}: {exc}")


# ---------------------------------------------------------------------------
#  The physics choices must actually reach the script
# ---------------------------------------------------------------------------


def _run_py(structure, **overrides):
    return bundle.unpack(bundle.build(_cfg(**overrides), structure, "1aki.pdb").content)["run.py"].decode()


def _pixi(structure, **overrides):
    """The pixi.toml with comment lines removed.

    The comments legitimately mention packages the file does not depend on
    (explaining why 3.11 is pinned, for instance), so a substring check against
    the raw text tests the prose rather than the dependency list.
    """
    text = bundle.unpack(
        bundle.build(_cfg(**overrides), structure, "1aki.pdb").content
    )["pixi.toml"].decode()
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))


def _code(structure, **overrides):
    """The generated run.py as executable code only: no comments, no docstrings.

    Assertions about what the script *does* must not be satisfiable by a comment
    that happens to mention the same word. Six tests in this file originally
    passed or failed on their own prose before this helper existed.

    ast.unparse normalises string literals to single quotes, so assertions below
    are written in that style.
    """
    tree = ast.parse(_run_py(structure, **overrides))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:] or [ast.Pass()]
    return ast.unparse(tree)


def _function(structure, name, **overrides):
    """One function's code, for assertions about ordering within it."""
    tree = ast.parse(_code(structure, **overrides))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(f"the generated script has no {name}()")


def test_water_model_emits_geometry_not_the_model_name(structure):
    """OPC is placed with tip4pew geometry and parameterised from opc.xml.

    Emitting model='opc' would not raise -- addSolvent simply does not know that
    name -- so this asserts the distinction survives into the script.
    """
    code = _code(structure, water_model="opc")
    assert "model='tip4pew'" in code
    assert "amber14/opc.xml" in code
    assert "model='opc'" not in code


def test_tip3p_uses_its_own_geometry(structure):
    code = _code(structure, water_model="tip3p")
    assert "model='tip3p'" in code
    assert "amber14/tip3p.xml" in code


def test_membrane_emits_addmembrane_and_the_membrane_barostat(structure):
    code = _code(structure, use_membrane=True, protein_ff="charmm36.xml")
    assert "modeller.addMembrane(" in code
    assert "MonteCarloMembraneBarostat(" in code
    assert "modeller.addSolvent(" not in code


def test_soluble_system_never_emits_a_membrane_barostat(structure):
    code = _code(structure)
    assert "MonteCarloMembraneBarostat(" not in code
    assert "modeller.addSolvent(" in code


def test_implicit_solvent_skips_solvation_and_pressure(structure):
    code = _code(structure, solvent_mode="implicit")
    assert "modeller.addSolvent(" not in code
    assert "Barostat(" not in code
    assert "implicitSolvent=app.GBn2" in code


def test_hmr_appears_only_when_enabled(structure):
    assert "hydrogenMass" in _code(structure, use_hmr=True)
    assert "hydrogenMass" not in _code(structure, use_hmr=False, timestep=2.0)


def test_restraints_appear_only_when_enabled(structure):
    assert "CustomExternalForce" in _code(structure)
    assert "CustomExternalForce" not in _code(structure, use_positional_restraints=False)


def test_values_are_inlined_not_read_from_config(structure):
    """The generated script must be readable as a record of what was run."""
    code = _code(structure, temperature=300.0, timestep=3.0)
    assert "300.0 * unit.kelvin" in code
    assert "3.0 * unit.femtoseconds" in code
    assert "cfg['temperature']" not in code


def test_switching_function_reaches_the_system(structure):
    code = _code(structure, protein_ff="charmm36.xml", nonbonded_cutoff=1.2,
                 use_switching=True, switch_distance=1.0)
    assert "switchDistance=1.0 * unit.nanometer" in code


def test_named_platform_skips_benchmarking(structure):
    assert "benchmarking" not in _code(structure, platform="CPU")
    assert "benchmarking" in _code(structure, platform="auto")


def test_relaxation_precedes_benchmarking(structure):
    """Benchmarking unminimised coordinates sends every platform to NaN.

    A real bug, caught by running the script rather than reading it, so the
    ordering gets pinned. Checked inside main() because both names also appear
    as function definitions earlier in the file.
    """
    main = _function(structure, "main")
    assert main.index("relax(") < main.index("select_platform(")


def test_restraints_reference_relaxed_coordinates(structure):
    main = _function(structure, "main")
    assert main.index("relax(") < main.index("add_restraints(")


# ---------------------------------------------------------------------------
#  Environment
# ---------------------------------------------------------------------------


def test_apo_bundle_does_not_pull_ambertools(structure):
    """A protein-only run should not drag in a gigabyte it will never use."""
    assert "openmmforcefields" not in _pixi(structure)
    assert "openff-toolkit" not in _pixi(structure)


def test_ligands_pull_the_parameterisation_stack(structure):
    text = _pixi(structure, keep_heteroatoms=["HEM"], has_ligands=True)
    assert "openmmforcefields" in text
    assert "openff-toolkit" in text


def test_ff19sb_pulls_openmmforcefields_even_without_ligands(structure):
    """ff19SB ships with openmmforcefields rather than OpenMM itself."""
    assert "openmmforcefields" in _pixi(structure, protein_ff="amber/protein.ff19SB.xml")


def test_the_reason_for_an_extra_dependency_is_written_down(structure):
    """Someone reading a generated pixi.toml should learn why it pulls a gigabyte."""
    raw = bundle.unpack(
        bundle.build(_cfg(protein_ff="amber/protein.ff19SB.xml"), structure, "1aki.pdb").content
    )["pixi.toml"].decode()
    assert "ff19SB ships with openmmforcefields" in raw


def test_openmm_is_pinned_to_the_tested_series(structure):
    assert f'openmm = "{bundle.OPENMM_PIN}"' in _pixi(structure)


def test_hdf5_output_pulls_mdtraj_reporter(structure):
    assert "from mdtraj.reporters import HDF5Reporter" in _run_py(structure, traj_format="hdf5")
    assert "HDF5Reporter" not in _run_py(structure, traj_format="xtc")


# ---------------------------------------------------------------------------
#  Free-text parsing
# ---------------------------------------------------------------------------


def test_mutations_are_grouped_by_chain(structure):
    script = _run_py(structure, mutations="A:ALA-57-GLY\nB:VAL-2-ILE\nA:LYS-9-ARG")
    assert 'applyMutations(["ALA-57-GLY", "LYS-9-ARG"], "A")' in script
    assert 'applyMutations(["VAL-2-ILE"], "B")' in script


def test_a_malformed_mutation_is_rejected_with_a_useful_message(structure):
    with pytest.raises(BundleError, match="CHAIN:WT-RESID-MUT"):
        _build(structure, mutations="this is not a mutation")


def test_a_malformed_restraint_schedule_is_rejected(structure):
    with pytest.raises(BundleError, match="not a"):
        _build(structure, restraint_schedule="1000, banana, 0")


def test_box_vectors_need_three_numbers(structure):
    with pytest.raises(BundleError, match="three numbers"):
        _build(structure, box_sizing="vectors", box_vectors="6.0, 6.0")


def test_an_invalid_config_is_refused_before_rendering(structure):
    with pytest.raises(BundleError, match="does not validate"):
        _build(structure, padding=0.5, nonbonded_cutoff=1.0)


# ---------------------------------------------------------------------------
#  Documentation
# ---------------------------------------------------------------------------


def test_readme_describes_the_actual_run(structure):
    files = bundle.unpack(_build(structure, production_duration=50.0).content)
    readme = files["README.md"].decode()
    assert "50.0 ns" in readme
    assert "chmod +x" in readme
    assert ".fxa" in readme
    assert "flexappeal_test_job.command" in readme


def test_readme_lists_replicates_only_when_there_are_several(structure):
    assert "Replicates" in bundle.unpack(_build(structure, replicates=3).content)["README.md"].decode()
    assert "Replicates" not in bundle.unpack(_build(structure).content)["README.md"].decode()


# ---------------------------------------------------------------------------
#  Ligands
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def trypsin():
    return (FIXTURES / "3ptb.pdb").read_bytes()


def _report(raw, name):
    from flexappeal import structure as structure_mod
    return structure_mod.analyse(raw, name).to_dict()


def _fake_ccd(name):
    """A stand-in for the RCSB fetch, so the suite stays offline."""
    from flexappeal import sources
    return sources.FetchResult(b"fake sdf\n", f"{name}.sdf", "ccd", "", f"CCD {name}")


def test_a_metal_cofactor_is_refused_with_the_real_reason(structure):
    """HEM cannot go through any organic force field, and saying so early matters.

    Left to run, this fails minutes in with a message about a missing residue
    template, which does not tell the user that the problem is the iron.
    """
    report = _report(structure, "1aki.pdb")
    report["heteroatoms"] = [{"name": "HEM", "category": "cofactor",
                              "elements": ["C", "FE", "N", "O"], "count": 1}]
    with pytest.raises(BundleError, match="transition-metal"):
        bundle.collect_ligands({"keep_heteroatoms": ["HEM"]}, report, fetch=_fake_ccd)


def test_an_organic_ligand_is_accepted(trypsin):
    files, warnings = bundle.collect_ligands(
        {"keep_heteroatoms": ["BEN"], "ph": 7.4}, _report(trypsin, "3ptb.pdb"),
        fetch=_fake_ccd)
    assert "BEN.sdf" in files
    assert any("protonation" in w for w in warnings), \
        "the CCD's fixed protonation state is a real caveat and must be surfaced"


def test_monatomic_ions_need_no_chemical_definition(structure):
    """The protein force field already has ion parameters."""
    report = _report(structure, "1aki.pdb")
    report["heteroatoms"] = [{"name": "NA", "category": "ion",
                              "elements": ["NA"], "count": 1}]
    files, _ = bundle.collect_ligands({"keep_heteroatoms": ["NA"]}, report,
                                      fetch=_fake_ccd)
    assert files == {}


def test_a_ligand_bundle_carries_its_chemistry(trypsin):
    files, _ = bundle.collect_ligands({"keep_heteroatoms": ["BEN"]},
                                      _report(trypsin, "3ptb.pdb"), fetch=_fake_ccd)
    result = bundle.build(_cfg(keep_heteroatoms=["BEN"], has_ligands=True),
                          trypsin, "3ptb.pdb", ligands=files)
    unpacked = bundle.unpack(result.content)
    assert "ligands/BEN.sdf" in unpacked


def test_ligand_run_script_stamps_residue_metadata_on_new_hydrogens(structure):
    """Chem.AddHs leaves the new hydrogens without residue info.

    They then land in a separate UNK residue: `resname BEN` selects 9 atoms
    instead of 17, and a residue with zero heavy atoms breaks MDTraj's
    closest-heavy contact analysis. Caught by running it, so it gets a test.
    """
    code = _code(structure, keep_heteroatoms=["BEN"], has_ligands=True)
    assert "AtomPDBResidueInfo()" in code
    assert "SetResidueName" in code


def test_am1bcc_is_checked_for_before_it_is_needed(structure):
    """AmberTools is found on PATH, not by import, so an unactivated environment
    fails minutes later with an opaque toolkit-registry error."""
    code = _code(structure, keep_heteroatoms=["BEN"], has_ligands=True,
                 ligand_charge_method="am1bcc")
    assert "shutil.which('sqm')" in code or 'shutil.which("sqm")' in code


def test_am1bcc_declares_ambertools_explicitly(structure):
    """It arrives transitively today; a hard requirement should not rely on that."""
    text = _pixi(structure, keep_heteroatoms=["BEN"], has_ligands=True,
                 ligand_charge_method="am1bcc")
    assert "ambertools" in text


def test_gasteiger_needs_no_ambertools(structure):
    text = _pixi(structure, keep_heteroatoms=["BEN"], has_ligands=True,
                 ligand_charge_method="gasteiger")
    assert "ambertools" not in text


@pytest.mark.parametrize("ff,generator", [
    ("openff-2.2.1", "SMIRNOFFTemplateGenerator"),
    ("gaff-2.11", "GAFFTemplateGenerator"),
    ("espaloma-0.3.2", "EspalomaTemplateGenerator"),
])
def test_each_ligand_force_field_selects_its_generator(structure, ff, generator):
    code = _code(structure, keep_heteroatoms=["BEN"], has_ligands=True, ligand_ff=ff)
    assert generator in code


def test_apo_run_has_no_ligand_machinery(structure):
    code = _code(structure)
    assert "TemplateGenerator" not in code
    assert "AssignBondOrdersFromTemplate" not in code


# ---------------------------------------------------------------------------
#  Membrane
# ---------------------------------------------------------------------------


def _membrane_cfg(**overrides):
    return _cfg(use_membrane=True, protein_ff="charmm36.xml",
                input_source="opm", pdb_id="1AFO", **overrides)


def test_membrane_uses_the_membrane_barostat_not_the_isotropic_one(structure):
    """An isotropic barostat squeezes the bilayer plane and the normal together,
    which is not physical for a membrane."""
    code = _code(structure, **{k: v for k, v in _membrane_cfg().items()
                               if k not in ("job_name",)})
    assert "MonteCarloMembraneBarostat(" in code
    assert "openmm.MonteCarloBarostat(" not in code


def test_membrane_normalisation_cannot_be_overridden_by_accident(structure):
    """Choosing a membrane and an isotropic barostat is a contradiction; the
    normaliser resolves it rather than letting it through."""
    from flexappeal import schema
    cfg = schema.normalise(opts.defaults() | {"use_membrane": True,
                                              "barostat": "MonteCarlo"})
    assert cfg["barostat"] == "MonteCarloMembrane"


def test_an_amber_membrane_loads_lipid_parameters_where_they_are_missing(structure):
    """ff19SB is a protein force field on its own and needs lipid17 beside it."""
    code = _code(structure, use_membrane=True,
                 protein_ff="amber/protein.ff19SB.xml")
    assert "amber14/lipid17.xml" in code


def test_a_force_field_that_has_the_lipids_is_not_given_them_twice(structure):
    """This test asserted the fault rather than the fix.

    It read "an AMBER force field needs lipid17" and checked that
    `amber14-all.xml` loaded it, which is exactly what OpenMM refuses:
    `amber14-all.xml` already carries the lipid17 residue templates, so the
    second registration raises "Residue template DLPC with the same override
    level 0 already exists". It is the default protein force field, so every
    membrane bundle failed on its force field minutes into the run.

    CHARMM36 carries its own lipids and AMOEBA carries everything, and the real
    check is `test_every_offered_force_field_builds_a_membrane`, which
    constructs each one rather than reading the list.
    """
    for ff in ("amber14-all.xml", "charmm36.xml", "amoeba2018.xml"):
        assert "lipid17" not in _code(structure, use_membrane=True,
                                      protein_ff=ff), ff


def test_lipid_residue_names_cover_the_pdb_truncation(structure):
    """A PDB residue name field is three characters wide.

    Matching only 'POPC' means a lipid read back from PDB (as 'POP') never
    matches, so selections meant to exclude the bilayer silently include it.
    """
    code = _code(structure, use_membrane=True, protein_ff="charmm36.xml")
    for full, truncated in (("POPC", "POP"), ("DPPC", "DPP"), ("DOPC", "DOP")):
        assert f"'{full}'" in code, f"{full} missing from the lipid set"
        assert f"'{truncated}'" in code, f"{truncated} missing from the lipid set"


def test_membrane_analysis_handles_both_lipid_spellings(structure):
    analyse = bundle.unpack(
        bundle.build(_cfg(use_membrane=True, protein_ff="charmm36.xml",
                          analysis_metrics=["membrane_apl", "membrane_scd"]),
                     structure, "1aki.pdb").content
    )["analyse.py"].decode()
    assert "def lipid_selection(" in analyse
    assert "name[:3]" in analyse


def test_a_soluble_run_emits_no_membrane_code(structure):
    code = _code(structure)
    assert "addMembrane" not in code
    assert "lipid_order" not in code


# ---------------------------------------------------------------------------
#  Progress reporting
# ---------------------------------------------------------------------------


def test_every_slow_stage_reports_progress(structure):
    """A stage that prints nothing for minutes is indistinguishable from a hang."""
    code = _code(structure)
    for machinery in ("class Progress", "class Spinner",
                      "class MinimisationProgress", "def run_steps", "def run_chunked"):
        assert machinery in code, f"{machinery} missing"

    # The stages that used to be silent.
    assert "Spinner('solvating" in code or 'Spinner("solvating' in code
    assert "Spinner('building the system" in code or 'Spinner("building the system' in code
    assert "MinimisationProgress('minimising'" in code or 'MinimisationProgress("minimising"' in code
    assert "run_steps(simulation" in code or "run_chunked(simulation" in code


def test_minimisation_reporter_takes_four_arguments(structure):
    """OpenMM calls report(iteration, x, grad, args).

    Python's introspection of the base class reports three parameters, and a
    three-argument override fails with a TypeError the moment minimisation
    starts -- after solvation, so several minutes into a run.
    """
    import ast as ast_mod

    tree = ast_mod.parse(_run_py(structure))
    for node in ast_mod.walk(tree):
        if isinstance(node, ast_mod.FunctionDef) and node.name == "report":
            names = [a.arg for a in node.args.args]
            assert names == ["self", "iteration", "x", "grad", "args"], names
            return
    raise AssertionError("no report() override found")


def test_the_minimisation_energy_is_read_by_subscript(structure):
    """`args` is OpenMM's mapstringdouble: subscriptable, but with no .get().

    Guarding with hasattr(args, "get") silently skipped the energy on every
    call and the progress line came out blank.
    """
    code = _code(structure)
    assert "args['system energy']" in code or 'args["system energy"]' in code
    assert "args.get(" not in code


def test_minimisation_progress_does_not_claim_a_percentage(structure):
    """OpenMM's L-BFGS restarts, so `iteration` is not monotonic.

    Measured on a real minimisation: 1530 calls, three resets, a final value of
    29 after peaking at 499. A bar driven by it runs backwards.
    """
    code = _code(structure)
    # Bounded at the next top-level definition: MinimisationProgress is the last
    # class in the file, so splitting only on "class " runs to the end and picks
    # up run_steps' own use of Progress.
    import re as _re

    tail = code.split("class MinimisationProgress", 1)[1]
    cut = _re.search(r"\n(?:class |def )", tail)
    cls = tail[: cut.start()] if cut else tail
    assert "self.iterations += 1" in cls, "the count must be its own, not the argument"

    # State the property rather than banning a spelling. This used to assert
    # 'Progress(' was absent, which caught the old hand-rolled bar but began
    # matching `RichProgress(` -- the spinner -- the moment rich arrived, even
    # though the spinner is exactly what it was asking for. What must not appear
    # is a *determinate* bar: a BarColumn, or a task given a finite total.
    assert "BarColumn" not in cls, "a percentage bar cannot track a resetting counter"
    assert "total=None" in cls, "the minimisation task must be indeterminate"


def test_production_has_no_scrolling_state_reporter(structure):
    """A StateDataReporter to stdout prints a row per interval, for hours."""
    code = _code(structure)
    assert "sys.stdout" not in code.split("def simulate", 1)[1].split("StateDataReporter")[0] \
        or "StateDataReporter(sys.stdout" not in code
    # The CSV one must survive: it is the data behind the convergence panels.
    assert "state_data.csv" in code


def test_progress_degrades_when_output_is_not_a_terminal(structure):
    """Redrawing with \\r into a log file makes one unreadable multi-megabyte line."""
    code = _code(structure)
    for cls in ("class Progress", "class Spinner", "class MinimisationProgress"):
        body = code.split(cls, 1)[1].split("\nclass ", 1)[0]
        assert "_TTY" in body, f"{cls} does not check whether stdout is a terminal"


# ---------------------------------------------------------------------------
#  Instability
# ---------------------------------------------------------------------------


def test_a_nan_is_diagnosed_rather_than_raised(structure):
    """OpenMM's "Particle coordinate is NaN" says what, never why.

    It also arrives as a traceback through OpenMM's internals, twenty minutes
    into a run. The script catches it and names the stage and the likely cause.
    """
    code = _code(structure)
    assert "def explain_instability" in code
    assert "def check_stable" in code
    assert "except openmm.OpenMMException" in code
    assert "'NaN' not in str(exc)" in code or '"NaN" not in str(exc)' in code


def test_stability_is_checked_at_every_stage_boundary(structure):
    """So the message names the stage, instead of a frame inside integrator.step."""
    code = _code(structure)
    for stage in ("minimisation", "heating", "equilibration"):
        assert f"check_stable(simulation, '{stage}')" in code \
            or f'check_stable(simulation, "{stage}")' in code, f"{stage} unchecked"


def test_the_diagnosis_names_a_rebuilt_loop(structure):
    """A long rebuilt loop is the commonest cause and the hardest to guess at."""
    code = _code(structure)
    body = code.split("def explain_instability", 1)[1].split("\ndef ", 1)[0]
    assert "REBUILT['longest']" in body or 'REBUILT["longest"]' in body
    assert "rebuilt" in body


def test_a_long_rebuilt_loop_is_flagged_at_preparation(structure):
    """Said once when it happens, not only after it has ruined the run."""
    code = _code(structure)
    assert "an invention, not a measurement" in code


def test_positive_energy_after_minimisation_warns(structure):
    """A solvated system's potential energy is strongly negative."""
    code = _code(structure)
    assert "still positive after minimisation" in code


def test_an_implicit_solvent_run_does_not_use_the_solvated_energy_check(structure):
    """The threshold is only meaningful with explicit water in the box."""
    assert "still positive after minimisation" not in _code(
        structure, solvent_mode="implicit")


def _template_theme(text):
    """The THEME literal out of a generated script.

    Those two are plain string literals, so they parse. console.py's is built
    from a BRAND dict and is not a literal, which is why it is imported as a
    live object rather than parsed.
    """
    import ast
    import re

    m = re.search(r"THEME = Theme\((\{.*?\})\)", text, re.S)
    assert m, "no THEME literal found"
    return ast.literal_eval(m.group(1))


def test_the_three_theme_copies_have_not_drifted():
    """The CLI, run.py and analyse.py each carry the palette.

    They have to be separate copies: the generated scripts run inside the
    bundle's own pixi environment, which has rich but no FlexAppeal package to
    import a shared theme from. Copies drift silently, so this compares them.
    The generated pair may define fewer keys than the CLI -- analyse.py has no
    remaining-time column -- but every key they share must agree.
    """
    import pathlib

    from rich.style import Style

    from flexappeal.console import THEME as cli_theme

    root = pathlib.Path(__file__).resolve().parent.parent
    for name in ("run.py.j2", "analyse.py.j2"):
        theme = _template_theme((root / "flexappeal" / "runtime" / name).read_text())
        assert theme, f"{name} has no theme"
        for key, value in theme.items():
            assert key in cli_theme.styles, f"{name} defines {key!r}, console.py does not"
            assert Style.parse(value) == cli_theme.styles[key], (
                f"{name} styles {key!r} as {value!r}, "
                f"console.py as {cli_theme.styles[key]!r}")


def test_the_bundle_declares_the_console_dependencies(structure):
    """rich and psutil are imported by both generated scripts.

    They arrive transitively in the dev environment, which is exactly why the
    bundle has to name them: the bundle solves its own environment from this
    file alone, and an unlisted import is an ImportError on the user's machine.
    """
    files = bundle.unpack(_build(structure).content)
    pixi = files["pixi.toml"].decode()
    run = files["run.py"].decode()
    analyse = files["analyse.py"].decode()
    for dep in ("rich", "psutil"):
        assert f"{dep} = " in pixi, f"the bundle environment does not declare {dep}"
    assert "import psutil" in run and "import psutil" in analyse
    assert "from rich" in run and "from rich" in analyse


def test_the_run_projects_memory_growth_rather_than_assuming_a_rate(structure):
    """The growth warning must measure, not hard-code.

    OpenMM's OpenCL platform leaks about 3 kB per integration step on an M1 Max
    -- measured with no reporters attached, so it is the platform and not this
    script. The rate depends on the system, the driver and the OpenMM build, so
    a constant baked in today would be wrong on someone else's machine. The
    script samples its own resident size and extrapolates instead.
    """
    code = _code(structure)
    assert "def watch_growth" in code
    tail = code.split("def watch_growth", 1)[1]
    body = tail[: tail.index("\ndef ")]

    assert "memory_info()" in body, "the projection must sample real memory"
    assert "_GROWTH" in body, "it needs a baseline to extrapolate from"
    # No magic per-step constant: the numbers present should be fractions and
    # thresholds, not a bytes-per-step rate.
    import re

    literals = [float(n) for n in re.findall(r"\b\d+\.?\d*\b", body)]
    assert not any(n > 100 for n in literals), \
        f"a hard-coded rate looks to have crept in: {[n for n in literals if n > 100]}"


def test_the_growth_warning_names_the_way_out(structure):
    """A warning with no remedy just tells someone their run is doomed.

    Restarting works because the bundle resumes from its checkpoint, which is
    not obvious, so the message has to say it.
    """
    run = bundle.unpack(_build(structure).content)["run.py"].decode()
    assert "resumes from the last checkpoint" in run
    assert "cannot be reclaimed while the process lives" in run


# ---------------------------------------------------------------------------
#  The command line, which has to default the same way the form does
# ---------------------------------------------------------------------------


def test_the_cli_keeps_a_ligand_and_ships_its_chemistry(trypsin, tmp_path,
                                                        monkeypatch, capsys):
    """`./FlexAppeal.py bundle 3ptb.pdb` must produce the complex, not the apo run.

    This path set no heteroatoms at all, so a structure with a ligand in it
    silently became an apo simulation, and the bundle carried no chemical
    definition for the run script to transfer bond orders from.
    """
    import FlexAppeal
    from flexappeal import sources

    monkeypatch.setattr(sources, "fetch_ccd_ideal", _fake_ccd)
    source = tmp_path / "3ptb.pdb"
    source.write_bytes(trypsin)
    out = tmp_path / "trypsin.command"

    code = FlexAppeal.main(["bundle", str(source), "-o", str(out), "-n", "trypsin"])
    assert code == 0

    files = bundle.unpack(out.read_bytes())
    assert "ligands/BEN.sdf" in files
    cfg = json.loads(files["config.json"])["config"]
    assert cfg["keep_heteroatoms"] == ["BEN"]
    assert cfg["has_ligands"] is True
    # The ion was left out, and the reason was printed rather than swallowed.
    assert "CA" not in cfg["keep_heteroatoms"]
    assert "is an ion" in capsys.readouterr().out


def test_an_explicit_config_still_wins_on_the_cli(trypsin, tmp_path, monkeypatch):
    """Auto-selection is a default, not an override.

    A config from a previous bundle already says what to keep, including when
    the answer is nothing.
    """
    import FlexAppeal
    from flexappeal import sources

    monkeypatch.setattr(sources, "fetch_ccd_ideal", _fake_ccd)
    source = tmp_path / "3ptb.pdb"
    source.write_bytes(trypsin)
    config = tmp_path / "apo.json"
    config.write_text(json.dumps({"config": {"keep_heteroatoms": []}}))
    out = tmp_path / "apo.command"

    assert FlexAppeal.main(["bundle", str(source), "-c", str(config),
                            "-o", str(out), "-n", "apo"]) == 0
    files = bundle.unpack(out.read_bytes())
    assert not [name for name in files if name.startswith("ligands/")]
    assert json.loads(files["config.json"])["config"]["keep_heteroatoms"] == []


# ===========================================================================
#  The measurement library a ligand bundle carries
# ===========================================================================
#
# Spec section 0: a run with a ligand emits version 2, and a one-complex package
# is a first-class graded run. The bundle runs on someone else's machine with no
# FlexAppeal checkout on it, so the measurement code has to travel with it.


def test_a_ligand_bundle_carries_the_measurement_library(trypsin):
    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)

    for name in bundle.VENDORED_MODULES:
        assert f"{bundle.VENDOR_PACKAGE}/{name}" in files, name
    # The modules import each other relatively, so the package marker has to be
    # there or `from . import fxa2` fails on the user's machine.
    assert f"{bundle.VENDOR_PACKAGE}/__init__.py" in files


def test_the_vendored_modules_are_byte_identical_to_the_installed_ones(trypsin):
    """A copy that can drift is worse than no copy.

    These are 220 KB of measurement code that must behave the same whether it
    runs on the server or in a bundle. Templating them, or keeping an edited
    copy, would mean two versions of the grade within one release.
    """
    from pathlib import Path

    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)

    source_dir = Path(bundle.__file__).resolve().parent
    for name in bundle.VENDORED_MODULES:
        assert files[f"{bundle.VENDOR_PACKAGE}/{name}"] == (
            source_dir / name).read_bytes(), f"{name} has drifted from the package"


def test_an_apo_bundle_carries_none_of_it(structure):
    """An apo run writes a version 1 .fxa and that is the complete answer for it
    (spec section 0), so there is nothing to measure a pocket for and no reason
    to grow the bundle by 220 KB."""
    cfg = _cfg(job_name="apo", keep_heteroatoms=[], has_ligands=False)
    files = bundle.unpack(bundle.build(cfg, structure, "1aki.pdb").content)
    assert not [n for n in files if n.startswith(f"{bundle.VENDOR_PACKAGE}/")]


def test_the_ligand_analysis_script_builds_the_version_two_package(trypsin):
    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)
    script = files["analyse.py"].decode()

    assert "def package_v2(" in script
    assert "from fxalib import package as packlib" in script
    # Measured from the run directory, never from the version 1 file: that
    # payload is decimated for viewing and a grade must not rest on it.
    assert "packlib.build(" in script
    compile(script, "analyse.py", "exec")


def test_the_apo_analysis_script_has_no_package_step(structure):
    cfg = _cfg(job_name="apo", keep_heteroatoms=[], has_ligands=False)
    files = bundle.unpack(bundle.build(cfg, structure, "1aki.pdb").content)
    script = files["analyse.py"].decode()
    assert "package_v2" not in script
    compile(script, "analyse.py", "exec")


def test_a_failed_package_does_not_cost_the_version_one_file(trypsin):
    """The version 1 file is written first and the package step is wrapped.

    A pocket measurement that raises on someone's machine must cost them the
    graded package and not the eight hours of simulation behind it.
    """
    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)
    script = files["analyse.py"].decode()

    body = script[script.index("def package_v2("):]
    assert "except Exception" in body
    assert "could not be built" in body
    assert "is unaffected" in body
    # The package step is called after the version 1 files are packed. The
    # energy, which both files carry, is computed before them: the version 1
    # metrics have to contain it by the time pack() runs or the Analysis tab
    # never sees a binding energy at all.
    assert script.index("binding_energy(run_dir)") < script.index("written.append(pack(")
    assert script.index("written.append(pack(") < script.rindex("package_v2(run_dirs,")


def test_version_two_is_the_default_for_a_ligand_run_not_an_option(trypsin):
    """No flag, no config key, no opt-in: spec section 0 makes it the default.

    A run with a ligand emits version 2 because version 1 has nowhere to put a
    ligand measurement, so the grade would be lost rather than deferred. Putting
    that behind a switch would mean most runs never produced one, which is the
    situation section 0 was written to end.
    """
    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)
    script = files["analyse.py"].decode()

    # Called unconditionally from main(), not guarded by an option.
    assert "package_v2(run_dirs," in script
    for switch in ("if PACKAGE", "--package", "args.package", "WRITE_PACKAGE",
                   "if cfg.get(\"package\")"):
        assert switch not in script, f"version 2 is behind {switch}"

    # And the config that produced it carries no key that could turn it off.
    config = json.loads(files["config.json"])["config"]
    assert not [k for k in config if "fxa_version" in k or k == "package"]


def test_the_energy_reaches_both_files_from_one_calculation(trypsin):
    """Computed once, before either file is written.

    The Analysis tab reads the version 1 payload and nothing else, so an energy
    that exists only in the package is an energy its reader never sees. And the
    frames are subsampled, so computing it twice would be two samples of one
    trajectory disagreeing in the last decimal.
    """
    cfg = _cfg(job_name="lig", keep_heteroatoms=["BEN"], has_ligands=True)
    files = bundle.unpack(bundle.build(cfg, trypsin, "3ptb.pdb").content)
    script = files["analyse.py"].decode()

    # Once as the def, once as the call, and no more: a second call site would
    # be a second subsample of the same trajectory.
    assert script.count("binding_energy(run_dir)") == 2, "computed more than once"
    assert script.count("def binding_energy(run_dir)") == 1
    assert "metrics.update(energy_metrics(" in script
    # Handed to the package rather than recomputed there.
    assert "energies=supplied or None" in script
    # Comments are stripped first: the template explains in prose why it passes
    # `energies` rather than `mmgbsa=True`, and matching that would be matching
    # the explanation instead of the code.
    code = "\n".join(line for line in script.splitlines()
                     if not line.lstrip().startswith("#"))
    assert "mmgbsa=True" not in code, "the package would compute it a second time"


# ---------------------------------------------------------------------------
#  Reference states, as a request rather than a step
# ---------------------------------------------------------------------------


def test_the_reference_state_options_are_off_and_silent_by_default(structure):
    """A run that did not ask for them should not hear about them."""
    from flexappeal import schema

    cfg = schema.validate({"job_name": "plain"}).config
    assert cfg["use_refstates"] is False

    files = bundle.unpack(bundle.build(cfg, structure).content)
    assert "reference states were requested" not in files["analyse.py"].decode()
    assert "## Reference states" not in files["README.md"].decode()


def test_the_request_travels_and_is_scored_by_the_bundle(structure):
    """The annotation is part of the run now, not a command printed after it.

    Resolving which residues a motif is needs KLIFS or GPCRdb, their deposited
    structures and gemmi; scoring a frame against the answer is numpy. So the
    first half runs where the bundle is built and its answer travels, which is
    what the ligand chemistry beside it already does.
    """
    import ast

    from flexappeal import schema

    cfg = schema.validate({
        "job_name": "abl1_demo", "use_refstates": "on",
        "refstates_family": "kinase", "refstates_protein": "ABL1",
        "refstates_reference_pdb": "2HYY", "refstates_max_references": 25,
    }).config
    files = bundle.unpack(bundle.build(
        cfg, structure, structure_filename="2hyy.pdb",
        refstates=b'{"version":1,"kind":"kinase"}').content)

    # The request is in the config the results file will carry.
    assert json.loads(files["config.json"])["config"]["refstates_protein"] == "ABL1"
    # The resolved ensemble, and the code that scores against it.
    assert files[bundle.REFSTATES_MEMBER] == b'{"version":1,"kind":"kinase"}'
    assert "fxalib/refstates/traj.py" in files

    script = files["analyse.py"].decode()
    ast.parse(script)
    assert "def reference_states(" in script
    assert "refstates annotate" not in script, (
        "the bundle scores this itself; a command to run afterwards would be "
        "telling the user to do what has already been done")

    readme = files["README.md"].decode()
    assert "## Reference states" in readme
    assert "deposited number" in readme, (
        "the renumbering is the step a reader has to understand")


def test_an_apo_run_can_be_scored_too(structure):
    """Two of the four library rows are apo, and both have reference states.

    The scoring step was written beside the binding energy at first, which put
    it inside the ligand-only section of the template: a receptor simulated
    without a ligand would have had the option and no scoring.
    """
    import ast

    from flexappeal import schema

    cfg = schema.validate({"job_name": "apo", "use_refstates": "on",
                           "refstates_family": "gpcr"}).config
    assert cfg["has_ligands"] is False
    script = bundle.unpack(bundle.build(
        cfg, structure, refstates=b"{}").content)["analyse.py"].decode()
    ast.parse(script)
    assert "def reference_states(" in script


def test_a_membrane_build_does_not_hand_its_relaxation_to_the_gpu(structure):
    """This repository already knew, and the generated script did not.

    `addMembrane` builds its own Context to inflate the protein back into the
    lipid patch, and with no platform given OpenMM takes the fastest available.
    On a Mac that is OpenCL, which Apple implements in single precision only,
    and the relaxation dies at `integrator.step(20)` with "Particle coordinate
    is NaN". `tools/build_membrane.py` carries both that finding and its fix in
    a comment; `run.py.j2` carried neither, so no bundle could build a membrane
    on a Mac at all.

    The second half is the unit cell: a placeholder CRYST1 of 1 1 1 makes
    `addMembrane` size the box from the cell rather than from the protein, the
    receptor wraps onto itself through the periodic boundary, and the same
    relaxation produces the same NaN.
    """
    from flexappeal import schema

    cfg = schema.validate({"job_name": "membrane", "use_membrane": "on",
                           "lipid_type": "POPC",
                           "barostat": "MonteCarloMembrane"}).config
    script = bundle.unpack(bundle.build(cfg, structure).content)["run.py"].decode()

    block = script[script.index("modeller.addMembrane("):]
    block = block[: block.index(")\n") + 2]
    assert 'platform=openmm.Platform.getPlatformByName("CPU")' in block

    before = script[:script.index("modeller.addMembrane(")]
    assert "setUnitCellDimensions(None)" in before, (
        "the placeholder cell has to go before the membrane is built")

    # And a run with no membrane neither needs nor mentions either.
    plain = schema.validate({"job_name": "water"}).config
    plain_script = bundle.unpack(
        bundle.build(plain, structure).content)["run.py"].decode()
    # The call, not the word: a comment about which OpenMM calls get a spinner
    # names `addMembrane` in every script.
    assert "modeller.addMembrane(" not in plain_script
    assert "setUnitCellDimensions(None)" not in plain_script


def test_a_membrane_measurement_without_the_membrane_in_the_trajectory_warns():
    """The panels are otherwise absent with nothing to say why.

    The default trajectory selection is protein and ligands, which excludes
    lipids by name, and the membrane metrics are measured on the saved
    trajectory rather than on the system. The end-to-end receptor run asked for
    all three, ran all three steps, found no lipids and omitted all three keys:
    correct behaviour for a value that cannot be produced, and useless to a
    reader who has just waited for a run.
    """
    from flexappeal import schema

    asked = {"job_name": "x", "use_membrane": "on", "lipid_type": "POPC",
             "barostat": "MonteCarloMembrane",
             "analysis_metrics": ["rmsd", "membrane_apl"]}

    warned = schema.validate(asked)
    assert warned.ok, "a selection that loses the lipids is a warning, not a refusal"
    assert [i for i in warned.issues if i.option_id == "traj_selection"]

    kept = schema.validate({**asked, "traj_selection": "all"})
    assert not [i for i in kept.issues if i.option_id == "traj_selection"]

    # And a run that asked for no membrane metric is not told about them.
    quiet = schema.validate({**asked, "analysis_metrics": ["rmsd"]})
    assert not [i for i in quiet.issues if i.option_id == "traj_selection"]


def test_no_resname_selection_is_built_unquoted(structure):
    """One unquoted interpolation is enough to lose a ligand silently.

    There were six, across the package and the generated script, and a ligand
    whose code begins with a digit fell out of every one of them: the packed
    trajectory, the ligand metrics, the binding-site measurement and the
    energy. The guard is on the source rather than on a behaviour, because the
    failure is silence.
    """
    import re

    from flexappeal import schema

    import pathlib

    root = pathlib.Path(__file__).resolve().parent.parent
    roots = [root / "flexappeal" / name for name in
             ("measure.py", "energy.py", "package.py")]
    cfg = schema.validate({"job_name": "x", "has_ligands": True,
                           "keep_heteroatoms": ["4VO"]}).config
    generated = bundle.unpack(bundle.build(cfg, structure).content)
    sources = {str(p): p.read_text() for p in roots}
    sources["analyse.py"] = generated["analyse.py"].decode()

    for where, text in sources.items():
        bad = re.findall(r"resname \{[^}']*\}", text)
        assert not bad, f"{where} interpolates a residue name unquoted: {bad[:3]}"
