"""Resolve every `openmm=` reference in the option registry against the real
installed packages.

This is the test that stops the registry drifting away from the API it claims to
describe. Without it, an upstream rename surfaces as a NameError inside a
generated script on someone else's machine, hours into their run and far away
from anything we can see.

Skipped entirely when the MD stack is not installed, so the fast suite still
runs on the droplet and in a bare checkout. Run `./install.sh` to get the real
coverage.
"""

from __future__ import annotations

import importlib

import pytest

from flexappeal import options as opts

from conftest import FIXTURES

pytest.importorskip("openmm", reason="OpenMM is not installed; run ./install.sh")


def _resolve(dotted: str):
    """Resolve 'openmm.app.PME' by importing the longest importable prefix."""
    parts = dotted.split(".")
    module = None
    consumed = 0
    for i in range(len(parts), 0, -1):
        try:
            module = importlib.import_module(".".join(parts[:i]))
            consumed = i
            break
        except ImportError:
            continue
    if module is None:
        raise ImportError(f"no importable module prefix in {dotted!r}")

    obj = module
    for attr in parts[consumed:]:
        obj = getattr(obj, attr)
    return obj


SYMBOLS = sorted({o.openmm: o.id for o in opts.OPTIONS if o.openmm}.items())


@pytest.mark.parametrize("dotted,option_id", SYMBOLS)
def test_declared_symbol_resolves(dotted, option_id):
    if dotted.startswith(("pdbfixer.", "openmmforcefields.")):
        root = dotted.split(".")[0]
        pytest.importorskip(root, reason=f"{root} is not installed")
    assert _resolve(dotted) is not None, \
        f"option {option_id!r} claims to map to {dotted!r}, which does not resolve"


def test_force_field_files_exist():
    """Every force-field XML we offer must actually ship with OpenMM."""
    from openmm import app

    ff_option = opts.BY_ID["protein_ff"]
    for choice in ff_option.choices:
        try:
            app.ForceField(choice.value)
        except Exception as exc:  # noqa: BLE001 -- we want the file name in the message
            pytest.fail(
                f"protein force field {choice.value!r} ({choice.label}) does not load: {exc}"
            )


def test_every_water_model_has_a_mapping():
    declared = set(opts.BY_ID["water_model"].choice_values())
    mapped = set(opts.WATER_MODEL_XML)
    assert declared == mapped, (
        f"the water model choices and WATER_MODEL_XML have diverged: {declared ^ mapped}"
    )


@pytest.mark.parametrize("name", sorted(opts.WATER_MODEL_XML))
def test_water_model_actually_solvates(name):
    """Really run addSolvent for each water model on a tiny box.

    Asserting the XML file merely *loads* is not enough. addSolvent's `model=`
    argument accepts only five geometry names, and passing an unrecognised one
    is the kind of mistake that produces a system with the wrong water rather
    than an exception. The only check that catches that is doing it.
    """
    from openmm import app, unit

    xml, geometry = opts.WATER_MODEL_XML[name]
    base = "charmm36.xml" if name == "charmm_tip3p" else "amber14-all.xml"
    forcefield = app.ForceField(base, xml)

    # A real 12-residue protein fragment, not a capped dipeptide: CHARMM36
    # matches both ALAD and AANM against a capped dipeptide and refuses to
    # choose, which is a property of the fixture rather than of the mapping.
    pdb = app.PDBFile(str(FIXTURES / "peptide12.pdb"))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model=geometry, padding=0.6 * unit.nanometer)

    waters = [r for r in modeller.topology.residues() if r.name in ("HOH", "WAT")]
    assert waters, f"{name} added no water at all"

    expected_sites = {"tip3p": 3, "spce": 3, "tip4pew": 4}[geometry]
    assert len(list(waters[0].atoms())) == expected_sites, (
        f"{name} placed {len(list(waters[0].atoms()))}-site water via geometry "
        f"{geometry!r}; expected {expected_sites}"
    )

    # And the parameters must actually assign -- this is what would fail if the
    # XML and the geometry disagreed about how many sites a water has.
    forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)


def test_lipid_types_are_supported_by_addmembrane():
    """addMembrane only knows a fixed set of lipids; offering others would fail late."""
    from openmm import app

    ff = app.ForceField("charmm36.xml", "charmm36/water.xml")
    for choice in opts.BY_ID["lipid_type"].choices:
        # addMembrane looks the lipid up in its own residue templates; the
        # cheapest honest check is that the force field knows the residue name.
        assert any(
            choice.value in str(t) for t in ff._templates
        ), f"lipid {choice.value!r} has no CHARMM36 template and addMembrane will reject it"


def test_integrators_exist():
    import openmm

    for choice in opts.BY_ID["integrator"].choices:
        name = f"{choice.value}Integrator"
        assert hasattr(openmm, name), f"openmm has no {name}"


def test_barostats_exist():
    import openmm

    for choice in opts.BY_ID["barostat"].choices:
        if choice.value == "none":
            continue
        name = f"{choice.value}Barostat"
        assert hasattr(openmm, name), f"openmm has no {name}"


def test_nonbonded_methods_exist():
    from openmm import app

    for choice in opts.BY_ID["nonbonded_method"].choices:
        assert hasattr(app, choice.value), f"openmm.app has no {choice.value}"


def test_constraint_types_exist():
    from openmm import app

    for choice in opts.BY_ID["constraints"].choices:
        if choice.value == "None":
            continue
        assert hasattr(app, choice.value), f"openmm.app has no {choice.value}"


def test_state_data_reporter_accepts_every_declared_field():
    """The energy-log field picker must match StateDataReporter's real keyword names."""
    import inspect as py_inspect

    from openmm.app import StateDataReporter

    signature = py_inspect.signature(StateDataReporter.__init__)
    for choice in opts.BY_ID["state_fields"].choices:
        assert choice.value in signature.parameters, (
            f"StateDataReporter has no {choice.value!r} keyword; the energy-log field "
            f"list has drifted from the OpenMM API"
        )
