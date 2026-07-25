"""Structure introspection tests, against real deposited files.

Fixtures are real PDB entries rather than synthetic minimal files, because the
things that break a parser are exactly the things a hand-written fixture leaves
out: alternate conformations, unusual heteroatoms, chain breaks, multiple
identical chains.

  1AKI  hen egg-white lysozyme, 1.5 A. One chain, four disulfides, no ligand.
  4HHB  deoxyhaemoglobin. Four chains, four haems, iron, real quaternary structure.
"""

from __future__ import annotations

import pytest

from flexappeal import options as opts
from flexappeal import structure
from flexappeal.structure import StructureError

from conftest import FIXTURES


@pytest.fixture(scope="module")
def lysozyme():
    return structure.analyse((FIXTURES / "1aki.pdb").read_bytes(), "1aki.pdb")


@pytest.fixture(scope="module")
def haemoglobin():
    return structure.analyse((FIXTURES / "4hhb.pdb").read_bytes(), "4hhb.pdb")


# ---------------------------------------------------------------------------
#  Parsing
# ---------------------------------------------------------------------------


def test_rejects_a_file_that_is_not_a_structure():
    with pytest.raises(StructureError):
        structure.analyse(b"this is not a pdb file at all\n" * 10, "notes.txt")


def test_rejects_empty_input():
    with pytest.raises(StructureError):
        structure.analyse(b"", "empty.pdb")


def test_lysozyme_basics(lysozyme):
    assert lysozyme.models == 1
    protein = lysozyme.protein_chains
    assert len(protein) == 1
    assert protein[0].observed_residues == 129  # hen lysozyme is 129 residues
    assert protein[0].sequence.startswith("KVFGRCELAA")


def test_lysozyme_finds_all_four_disulfides(lysozyme):
    """Lysozyme has exactly four: 6-127, 30-115, 64-80, 76-94."""
    pairs = {(d.resid_a, d.resid_b) for d in lysozyme.disulfides}
    assert pairs == {(6, 127), (30, 115), (64, 80), (76, 94)}
    for d in lysozyme.disulfides:
        assert 1.9 <= d.distance <= 2.5


def test_lysozyme_has_crystallographic_water_but_no_ligand(lysozyme):
    categories = {h.name: h.category for h in lysozyme.heteroatoms}
    assert categories.get("HOH") == "water"
    assert not [h for h in lysozyme.heteroatoms if h.category == "ligand"]


def test_haemoglobin_has_four_chains_and_four_haems(haemoglobin):
    assert len(haemoglobin.protein_chains) == 4
    hem = [h for h in haemoglobin.heteroatoms if h.name == "HEM"]
    assert len(hem) == 1
    assert hem[0].count == 4
    assert hem[0].category == "cofactor"
    assert hem[0].description == "haem b"


def test_haemoglobin_warns_about_ligands_needing_parameters(haemoglobin):
    assert any("parameters" in w for w in haemoglobin.warnings)


def test_each_author_chain_appears_exactly_once(haemoglobin):
    """gemmi splits chains per entity; a user thinks of chain A as one chain.

    4HHB is the case that exposes this: each of its four chains carries a haem,
    so setup_entities() produces eight Chain objects for four author chains.
    """
    ids = [c.id for c in haemoglobin.chains]
    assert len(ids) == len(set(ids)), f"duplicate chain ids: {ids}"
    assert {c.id for c in haemoglobin.protein_chains} == {"A", "B", "C", "D"}


def test_merged_chain_keeps_its_polymer_identity(haemoglobin):
    chain_a = next(c for c in haemoglobin.chains if c.id == "A")
    assert chain_a.kind == "protein"
    assert chain_a.observed_residues == 141
    # The haem's atoms belong to chain A too, so the atom count must exceed
    # what the polymer alone contributes.
    assert chain_a.atoms > 141 * 7


def test_geometry_is_reported_in_nanometres(lysozyme):
    # Lysozyme is a compact ~3 nm globule; anything wildly outside that means
    # units have been mixed up somewhere.
    assert 2.0 < lysozyme.max_extent_nm < 6.0
    assert all(e > 0 for e in lysozyme.extent_nm)


def test_solute_atoms_excludes_water(lysozyme):
    water_atoms = sum(h.atoms for h in lysozyme.heteroatoms if h.category == "water")
    total_atoms = sum(c.atoms for c in lysozyme.chains)
    assert water_atoms > 0
    assert lysozyme.solute_atoms == total_atoms - water_atoms


def test_report_is_json_serialisable(lysozyme):
    import json

    json.dumps(lysozyme.to_dict())


# ---------------------------------------------------------------------------
#  Dynamic choices
# ---------------------------------------------------------------------------


def test_dynamic_choices_cover_the_registry_dynamic_options(haemoglobin):
    choices = structure.dynamic_choices(haemoglobin)
    dynamic_ids = {o.id for o in opts.OPTIONS if o.dynamic}
    assert set(choices) == dynamic_ids, \
        "structure.dynamic_choices and the registry's dynamic options have diverged"


def test_chain_choices_exclude_water_chains(haemoglobin):
    choices = structure.dynamic_choices(haemoglobin)["chains"]
    assert {c["value"] for c in choices} >= {"A", "B", "C", "D"}
    for c in choices:
        assert "water" not in c["help"]


def test_hetero_choices_flag_what_needs_parameters(haemoglobin):
    choices = structure.dynamic_choices(haemoglobin)["keep_heteroatoms"]
    hem = next(c for c in choices if c["value"] == "HEM")
    assert "needs small-molecule parameters" in hem["help"]


# ---------------------------------------------------------------------------
#  Size estimation
# ---------------------------------------------------------------------------


def test_size_estimate_is_in_the_right_ballpark(lysozyme):
    """Lysozyme in a 1.2 nm-padded dodecahedron is roughly 25-40k atoms.

    This is the standard tutorial system, so the expected size is well known --
    which makes it a real check on the estimator rather than a tautology.
    """
    est = structure.estimate_system_size(lysozyme, opts.defaults())
    assert 15_000 < est["total_atoms"] < 50_000, est
    assert est["water_molecules"] > 4_000
    assert est["box_nm"] > lysozyme.max_extent_nm


def test_dodecahedron_needs_fewer_waters_than_a_cube(lysozyme):
    cfg = opts.defaults()
    cube = structure.estimate_system_size(lysozyme, cfg | {"box_shape": "cube"})
    dodec = structure.estimate_system_size(lysozyme, cfg | {"box_shape": "dodecahedron"})
    # The geometric factor is 0.7071, so expect roughly 29% fewer.
    ratio = dodec["water_molecules"] / cube["water_molecules"]
    assert 0.6 < ratio < 0.8, f"expected ~0.71, got {ratio:.2f}"


def test_more_padding_means_more_water(lysozyme):
    cfg = opts.defaults()
    small = structure.estimate_system_size(lysozyme, cfg | {"padding": 1.0})
    large = structure.estimate_system_size(lysozyme, cfg | {"padding": 2.0})
    assert large["water_molecules"] > small["water_molecules"] * 1.5


def test_implicit_solvent_adds_no_water(lysozyme):
    est = structure.estimate_system_size(lysozyme, opts.defaults() | {"solvent_mode": "implicit"})
    assert est["water_molecules"] == 0
    assert est["total_atoms"] == lysozyme.solute_atoms


def test_four_site_water_has_more_atoms_than_three_site(lysozyme):
    cfg = opts.defaults()
    tip3p = structure.estimate_system_size(lysozyme, cfg | {"water_model": "tip3p"})
    opc = structure.estimate_system_size(lysozyme, cfg | {"water_model": "opc"})
    assert opc["total_atoms"] > tip3p["total_atoms"]


def test_membrane_estimate_includes_lipid(lysozyme):
    """Lysozyme is not a membrane protein, but the estimator must still be sane."""
    est = structure.estimate_system_size(
        lysozyme, opts.defaults() | {"use_membrane": True, "lipid_type": "POPC"}
    )
    assert est["lipid_count"] > 0
    assert est["lipid_atoms"] > 0
    assert est["total_atoms"] > lysozyme.solute_atoms + est["lipid_atoms"]


def test_estimate_feeds_the_derived_readouts(lysozyme):
    """The whole point of the estimate: it makes derive() give real numbers."""
    from flexappeal import schema

    cfg = opts.defaults()
    est = structure.estimate_system_size(lysozyme, cfg)
    cfg["_estimated_atoms"] = est["total_atoms"]
    cfg["_solute_atoms"] = lysozyme.solute_atoms

    wall = schema.estimate_wall_time(cfg)
    assert wall["basis"] == "estimated"
    assert wall["hours"] > 0
    assert schema.derive(cfg)["traj_bytes"] > 0
