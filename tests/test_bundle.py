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


def test_bundles_are_deterministic(structure):
    """The same configuration must produce the same bytes.

    Only the generation timestamp may differ, so compare the payloads rather
    than the whole file.
    """
    first = _build(structure).content.split(b"__FLEXAPPEAL_PAYLOAD__", 1)[1]
    second = _build(structure).content.split(b"__FLEXAPPEAL_PAYLOAD__", 1)[1]
    assert first == second


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


def test_unpack_rejects_a_truncated_payload(structure):
    content = _build(structure).content
    with pytest.raises(BundleError, match="corrupt or truncated"):
        bundle.unpack(content[: len(content) // 2])


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


def test_amber_membrane_loads_lipid_parameters(structure):
    """CHARMM36 carries its own lipids; an AMBER force field needs lipid17."""
    code = _code(structure, use_membrane=True, protein_ff="amber14-all.xml")
    assert "amber14/lipid17.xml" in code


def test_charmm_membrane_does_not_double_load_lipids(structure):
    code = _code(structure, use_membrane=True, protein_ff="charmm36.xml")
    assert "lipid17" not in code


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
