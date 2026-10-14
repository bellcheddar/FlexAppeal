"""MM/GBSA binding free energies.

Most of this is about the two ways this calculation fails quietly, because both
of them happened while it was being built and neither raised anything:

- A topology written from an OpenMM system has its GB radii all zero, so every
  solvation term comes back `nan` and MMPBSA.py still reports success. The first
  real run here produced a complete, well-formatted results file whose every
  total was `nan`.
- If the complex, receptor and ligand topologies disagree about the molecule, the
  bonded terms stop cancelling and the result is a plausible-looking number that
  is wrong. DeepFlex's own build log records a near-miss topology that loads
  against the wrong trajectory, so this is not hypothetical.

The parsing tests use canned MMPBSA.py output so they run everywhere. The one
end-to-end test needs AmberTools and the committed example, and skips without
them.
"""

from __future__ import annotations

import shutil

import pytest

from flexappeal import energy


# A real difference block, trimmed, from the benzene run.
_GOOD = """
GENERALIZED BORN:

Complex:
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
BOND                       520.4067               20.7223              4.6336
TOTAL                    -2499.5872               93.2336             20.8477


Differences (Complex - Receptor - Ligand):
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
BOND                        -0.0000                0.0000              0.0000
ANGLE                       -0.0000                0.0000              0.0000
DIHED                       -0.0000                0.0001              0.0000
VDWAALS                    -11.2868                7.6385              0.7639
EEL                         -1.1099                0.8477              0.0848
1-4 VDW                      0.0000                0.0001              0.0000
1-4 EEL                     -0.0000                0.0001              0.0000
EGB                          4.3853                2.9658              0.2966
ESURF                       -1.7661                1.1847              0.1185

DELTA G gas                -12.3967                8.3824              0.8382
DELTA G solv                 2.6192                1.7972              0.1797

DELTA TOTAL                 -9.7775                6.6430              0.6643
"""

# The same calculation with no GB radii on the topology: it completes, it is
# well formed, and every solvation number is nan.
_NAN = _GOOD.replace("4.3853                2.9658              0.2966",
                     "nan                   nan                 nan")

# Bonded terms that did not cancel, which means the three topologies are not the
# same molecule.
_UNCANCELLED = _GOOD.replace("BOND                        -0.0000                0.0000              0.0000",
                             "BOND                        12.4000                3.1000              0.3000")


# ===========================================================================
#  The protocol
# ===========================================================================


def test_the_protocol_matches_the_cohort_these_numbers_join():
    """igb=2 with mbondi2, because that is what produced the energies in the
    package this is ranked against.

    A percentile across a cohort measured two ways is not a percentile. igb=2
    paired with the default Bondi radii instead of mbondi2 is a known way to get
    numbers that look fine and are not, so both halves are pinned here.
    """
    assert energy.IGB == 2
    assert energy.RADII_SET == "mbondi2"


def test_the_cancelling_terms_are_the_internal_ones():
    """Only terms that are genuinely internal to the molecules can be expected
    to cancel. A van der Waals or electrostatic term that cancelled would mean
    the ligand was not interacting with anything."""
    assert set(energy.CANCELLING_TERMS) == {
        "BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL"}
    for term in ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA TOTAL"):
        assert term not in energy.CANCELLING_TERMS


def test_the_generator_mapping_agrees_with_the_bundle():
    """This module carries its own copy because it travels inside a bundle and
    `bundle` does not. A copy that can drift is the thing to guard."""
    from flexappeal import bundle

    assert energy.LIGAND_GENERATOR == bundle._LIGAND_GENERATOR


# ===========================================================================
#  Parsing
# ===========================================================================


def test_the_difference_block_is_read_and_the_absolute_blocks_are_not():
    """The three absolute blocks are thousands of kcal/mol and only their
    difference is a binding quantity -- which is also what the package writer
    refuses values the size of."""
    terms = energy._difference_block(_GOOD)
    assert terms["DELTA TOTAL"] == [-9.7775, 6.6430, 0.6643]
    assert terms["VDWAALS"][0] == pytest.approx(-11.2868)
    # The complex block's BOND is 520, the difference block's is 0. Reading the
    # wrong one would be a 500 kcal/mol error that still looks like a number.
    assert terms["BOND"][0] == pytest.approx(0.0)


def test_every_difference_term_is_carried():
    terms = energy._difference_block(_GOOD)
    for name in energy.DIFFERENCE_TERMS:
        assert name in terms, name


def test_a_nan_solvation_term_is_refused_and_named():
    """The failure that produced a complete results file full of nan."""
    with pytest.raises(energy.EnergyError) as caught:
        energy._difference_block(_NAN)
    message = str(caught.value)
    assert "nan" in message
    assert "GB radii" in message, "the message must name the cause"


def test_output_with_no_difference_block_is_refused():
    with pytest.raises(energy.EnergyError, match="no difference block"):
        energy._difference_block("MMPBSA.py started and then stopped")


def test_output_with_no_total_is_refused():
    truncated = _GOOD.split("DELTA G gas")[0]
    with pytest.raises(energy.EnergyError, match="no total"):
        energy._difference_block(truncated)


# ===========================================================================
#  The cancellation check
# ===========================================================================


def test_a_single_trajectory_result_passes_the_cancellation_check():
    energy._check_cancellation(energy._difference_block(_GOOD))


def test_uncancelled_bonded_terms_are_refused():
    """The error that produces a plausible number rather than a failure."""
    terms = energy._difference_block(_UNCANCELLED)
    with pytest.raises(energy.EnergyError) as caught:
        energy._check_cancellation(terms)
    message = str(caught.value)
    assert "did not cancel" in message
    assert "BOND" in message
    assert "same molecule" in message


def test_rounding_in_the_bonded_terms_is_tolerated():
    """MMPBSA.py prints four decimals, so an exact zero can print as -0.0001.
    The tolerance has to admit that and nothing larger."""
    terms = energy._difference_block(_GOOD)
    terms["DIHED"] = [0.009, 0.0, 0.0]
    energy._check_cancellation(terms)                      # under tolerance

    terms["DIHED"] = [0.5, 0.0, 0.0]
    with pytest.raises(energy.EnergyError):
        energy._check_cancellation(terms)


# ===========================================================================
#  Shape and degradation
# ===========================================================================


def test_the_result_is_shaped_as_supplied_energies():
    """So `package.assemble` needs nothing new: computed and supplied energies
    reach the file by one path, which is one path to keep right."""
    result = energy.Energies(
        compound="BNZ", receptor="181L", dg=-9.7775, sd=6.643, sem=0.6643,
        n_frames=200, terms={"VDWAALS": [-11.2868, 7.6385, 0.7639]})
    supplied = result.as_supplied()

    # Both keys, same entry. The receptor token is derived from the structure
    # filename or the accession, so a caller that derives it differently would
    # silently lose the energy; the bare-compound key makes the lookup unable to
    # miss, and `<compound>/<receptor>` is still the specific one.
    assert set(supplied) == {"BNZ/181L", "BNZ"}
    assert supplied["BNZ"] == supplied["BNZ/181L"]
    entry = supplied["BNZ/181L"]
    assert entry["mmpbsa_gb_dg"] == pytest.approx(-9.7775)
    assert entry["mmpbsa_gb_sd"] == pytest.approx(6.643)
    assert entry["terms"]["gb"]["VDWAALS"][0] == pytest.approx(-11.2868)
    # The absolute-block guard in the package writer refuses anything over 1000,
    # so a difference-block number has to be well inside it.
    assert abs(entry["mmpbsa_gb_dg"]) < 1000.0


def test_a_missing_toolchain_is_reported_rather_than_raised(monkeypatch):
    """The energy is one sub-score of six. Losing it must not cost the other
    five, and a Gasteiger-charged run has no AmberTools."""
    monkeypatch.setattr(energy.shutil, "which", lambda name: None)
    ok, reason = energy.available()
    assert ok is False
    assert "MMPBSA.py" in reason
    assert "AmberTools" in reason


def test_a_run_with_no_ligand_chemistry_cannot_be_scored(tmp_path):
    """Without the SDF the run was parameterised from, the charges here would
    not be the charges that were simulated."""
    import json

    (tmp_path / "config.json").write_text(json.dumps({"config": {}}))
    with pytest.raises(energy.EnergyError, match="no .sdf"):
        energy.forcefield_for_run(tmp_path)


def test_a_run_with_no_config_cannot_be_scored(tmp_path):
    with pytest.raises(energy.EnergyError, match="no config.json"):
        energy.forcefield_for_run(tmp_path)


def test_an_unknown_ligand_forcefield_is_refused(tmp_path):
    import json

    (tmp_path / "config.json").write_text(
        json.dumps({"config": {"ligand_ff": "something-made-up"}}))
    (tmp_path / "ligands").mkdir()
    (tmp_path / "ligands" / "LIG.sdf").write_text("")
    with pytest.raises(energy.EnergyError, match="not a ligand force field"):
        energy.forcefield_for_run(tmp_path)


def test_a_source_without_a_run_directory_is_warned_about_not_skipped():
    """A version 1 .fxa has no topology on disk, no ligand SDF and no config, so
    it cannot have an energy computed. That has to be said rather than leaving a
    reader to wonder why one complex has an energy and another does not."""
    from flexappeal import package

    source = package.Source(label="from_fxa1", origin="version 1 .fxa",
                            trajectory=None, config={"keep_heteroatoms": ["BEN"]})
    found, warnings = package._compute_energies([source])
    assert found == {}
    assert any("needs the run directory" in w for w in warnings)
    assert any("stays pending" in w for w in warnings)


# ===========================================================================
#  End to end
# ===========================================================================


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("MMPBSA.py") is None,
                    reason="AmberTools is not on PATH")
def test_the_committed_example_gives_a_physical_binding_energy():
    """Benzene in T4 lysozyme L99A, the canonical buried non-polar cavity.

    Checked for physics rather than for a number. Experiment is about
    -5.2 kcal/mol and single-trajectory MM/GBSA with no entropy term is
    systematically over-negative, so the window is wide -- but the decomposition
    has to say the right thing: van der Waals dominating, electrostatics almost
    absent for a molecule with no polar atom, a positive desolvation penalty and
    a favourable surface term.
    """
    from pathlib import Path

    pytest.importorskip("parmed")
    pytest.importorskip("openmmforcefields")

    run = Path("examples/lysozyme_benzene_100ns/run")
    out = run / "flexappeal_output/lysozyme_benzene_100ns"
    if not (out / "trajectory.xtc").is_file():
        pytest.skip("the example run directory is not present")

    supplied, warnings = energy.for_run_directory(
        run, out, receptor="181L", max_frames=20)
    assert supplied, warnings

    entry = supplied["BNZ/181L"]
    dg = entry["mmpbsa_gb_dg"]
    assert -40.0 < dg < 0.0, f"{dg} kcal/mol is not a physical binding energy"

    gb = entry["terms"]["gb"]
    assert gb["VDWAALS"][0] < gb["EEL"][0], "van der Waals must dominate here"
    assert abs(gb["EEL"][0]) < 5.0, "benzene has no polar atom"
    assert gb["EGB"][0] > 0.0, "desolvation on binding costs energy"
    assert gb["ESURF"][0] < 0.0, "burying non-polar surface is favourable"
    # And the signature of a single-trajectory calculation.
    for term in energy.CANCELLING_TERMS:
        assert abs(gb[term][0]) < energy.CANCELLATION_TOLERANCE, term


# ---------------------------------------------------------------------------
#  The solute, and only the solute
# ---------------------------------------------------------------------------


def _solvated_pdb(path, waters: int = 30):
    """A two-residue peptide, one ligand, and some water. Written here, not
    committed: a fixture is generated in the test."""
    lines = []
    serial = 0

    def atom(record, name, resname, resseq, element, xyz):
        nonlocal serial
        serial += 1
        lines.append(
            f"{record:<6s}{serial:>5d} {name:<4s}{resname:>4s} A{resseq:>4d}"
            f"    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
            f"          {element:>2s}")

    for index, (resname, names) in enumerate(
            [("ALA", ("N", "CA", "C", "O")), ("GLY", ("N", "CA", "C", "O"))], 1):
        for offset, name in enumerate(names):
            atom("ATOM", name, resname, index, name[0],
                 (index * 3.0 + offset * 0.5, 0.0, 0.0))
    for offset, name in enumerate(("C1", "C2", "C3")):
        atom("HETATM", name, "BEN", 3, "C", (10.0 + offset * 0.5, 1.0, 0.0))
    for index in range(waters):
        atom("HETATM", "O", "HOH", 100 + index, "O",
             (20.0 + index * 3.0, 5.0, 0.0))
    # The ligand's bonds, the way a file under 99,999 atoms carries them. A
    # topology that has them is left alone; one that lost them to the serial
    # limit is the case `_bond_the_ligands` exists for.
    lines += ["CONECT    9   10", "CONECT   10    9   11", "CONECT   11   10"]
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_the_solvent_is_stripped_before_the_energy(tmp_path):
    """GB is the solvent model, so explicit water is wrong before it is
    anything else.

    A bundle writes a topology already subset-matched to its trajectory and
    nothing is removed. A run packed from a queue's own solvated system carries
    the whole box, and the first of those failed with "No template found for
    residue 263 (HOH)", which names the water rather than the mistake: 43,862
    of its 48,185 atoms were solvent.
    """
    pytest.importorskip("mdtraj")

    solvated = _solvated_pdb(tmp_path / "topology.pdb")
    out, keep, stripped = energy._solute_only(solvated, ["BEN"], tmp_path)

    assert stripped == 30, "every water should go"
    assert out != solvated, "the original file is never rewritten"
    assert len(keep) == 11, "eight protein atoms and three of the ligand"

    import mdtraj as md

    kept = md.load(str(out))
    assert sorted({r.name for r in kept.topology.residues}) == ["ALA", "BEN", "GLY"]


def test_a_solute_topology_is_passed_through_untouched(tmp_path):
    """The usual case: a bundle's own topology.pdb, with nothing to remove."""
    pytest.importorskip("mdtraj")

    solute = _solvated_pdb(tmp_path / "topology.pdb", waters=0)
    out, keep, stripped = energy._solute_only(solute, ["BEN"], tmp_path)

    assert (out, keep, stripped) == (solute, None, 0)


def test_a_topology_with_no_complex_in_it_says_so(tmp_path):
    """Rather than scoring whatever happens to be left."""
    pytest.importorskip("mdtraj")

    water_only = _solvated_pdb(tmp_path / "water.pdb", waters=5)
    import mdtraj as md

    loaded = md.load(str(water_only))
    water = loaded.atom_slice(loaded.topology.select("water"))
    water.save_pdb(str(tmp_path / "justwater.pdb"))

    with pytest.raises(energy.EnergyError, match="no complex to score"):
        energy._solute_only(tmp_path / "justwater.pdb", ["BEN"], tmp_path)


def test_a_ligand_that_lost_its_bonds_gets_them_back(tmp_path):
    """A PDB atom serial is five columns wide.

    Past 99,999 atoms a CONECT record cannot name the atoms it means, so a
    reader gets nothing from it. The receptor run here is 187,848 atoms and its
    carazolol arrived with no bonds at all, while the kinase run at 48,185 kept
    all 72 of its imatinib's. The template generator matches a residue by its
    graph, so with no edges there was no match and OpenMM reported it as "no
    template found for residue 322 (CAU)", against the nearest amino acid it
    could think of.
    """
    pytest.importorskip("mdtraj")

    import mdtraj as md

    # Three carbons in a row at 1.5 Å, which is a bond, and a fourth atom 5 Å
    # away, which is not.
    pdb = tmp_path / "ligand.pdb"
    rows = []
    for serial, (name, x) in enumerate(
            [("C1", 0.0), ("C2", 1.5), ("C3", 3.0), ("C4", 8.0)], start=1):
        rows.append(
            f"HETATM{serial:>5d} {name:<4s} BEN A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
    pdb.write_text("\n".join(rows) + "\nEND\n")

    loaded = md.load(str(pdb))
    assert len(list(loaded.topology.bonds)) == 0, "the fixture starts bondless"

    added = energy._bond_the_ligands(loaded, ["BEN"])
    assert added == 2, "C1-C2 and C2-C3, and nothing to the atom 5 Å away"

    # And a second pass changes nothing: a ligand whose bonds were read is left
    # exactly as it was read.
    assert energy._bond_the_ligands(loaded, ["BEN"]) == 0


def test_a_complex_whose_halves_never_touched_is_refused():
    """A clean zero is the signature of a ligand in the wrong periodic copy.

    The bonded terms cancel perfectly in that case, so the check that exists
    for a mis-split mask passes it: carazolol sat 101.7 angstrom from its
    receptor in a 122 by 128 by 116 box and MM/GBSA reported every term as
    zero, with no error anywhere. There is no complex with no van der Waals
    and no electrostatics.
    """
    with pytest.raises(energy.EnergyError, match="never in contact"):
        energy._check_contact({"VDWAALS": [-0.0001, 0.0001, 0.0],
                               "EEL": [-0.021, 0.1019, 0.0079]})

    # A real one passes.
    energy._check_contact({"VDWAALS": [-74.0046, 3.36, 0.26],
                           "EEL": [-35.0976, 5.10, 0.39]})


def test_a_ligand_in_the_next_periodic_copy_is_imaged_back(tmp_path):
    """One shift per ligand per frame, from its own centroid."""
    pytest.importorskip("mdtraj")

    import mdtraj as md
    import numpy as np

    solute = _solvated_pdb(tmp_path / "topology.pdb", waters=0)
    traj = md.load(str(solute))
    box = np.array([30.0, 30.0, 30.0])
    traj.unitcell_lengths = np.array([box / 10.0])
    traj.unitcell_angles = np.array([[90.0, 90.0, 90.0]])

    ligand = traj.topology.select("resname BEN")
    protein = traj.topology.select("protein")
    before = float(np.linalg.norm(traj.xyz[0][ligand].mean(axis=0)
                                 - traj.xyz[0][protein].mean(axis=0))) * 10.0
    # One box along x, which is where whole-molecule wrapping leaves it.
    traj.xyz[0, ligand] += np.array([3.0, 0.0, 0.0], dtype="float32")
    assert energy._image_ligand(traj, ["BEN"]) == pytest.approx(30.0)

    after = float(np.linalg.norm(traj.xyz[0][ligand].mean(axis=0)
                                - traj.xyz[0][protein].mean(axis=0))) * 10.0
    assert after == pytest.approx(before, abs=1e-3), "back where it started"

    # And a ligand already in the right copy is left exactly alone.
    assert energy._image_ligand(traj, ["BEN"]) == 0.0
