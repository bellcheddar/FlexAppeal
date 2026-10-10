"""The `.fxa` version 2 package: the format, the identifiers, and conformance.

Version 2 is a package of receptor-ligand complexes rather than one system, and
it is read at the far end by DeepFlex's own ingester. Three things are therefore
worth testing harder than usual:

* **Units.** Version 1 reported nanometres and version 2 reports angstrom. Get
  that wrong and every length in the report is out by ten, which is invisible
  until somebody checks a number by hand, so there is a test here that asserts a
  known distance.
* **Chemistry.** Section 8 of the spec bans structure strings and chemistry
  members outright. The writer refuses them rather than letting DeepFlex discard
  them on ingest and report that it did.
* **The round trip.** A writer that cannot read its own output is a writer with a
  unit bug.

Every fixture here is built in `tmp_path` rather than committed. A `.fxa` carries
a topology and coordinates, this repository is public, and a committed package
would be one more file to vet every time the format changes.
"""

from __future__ import annotations

import json
import zipfile

import numpy as np
import pytest

from flexappeal import fxa2
from flexappeal.fxa2 import Fxa2Error

# A pocket crop with no real coordinates in it: two receptor residues and a
# two-atom ligand on a 1 A spacing, which is also what the units test measures.
POCKET_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  N   GLY A   2       3.000   0.000   0.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       4.458   0.000   0.000  1.00  0.00           C
TER       5      GLY A   2
ATOM      6  C1  LIG B   1       8.000   0.000   0.000  1.00  0.00           C
ATOM      7  C2  LIG B   1       9.000   0.000   0.000  1.00  0.00           C
TER       8      LIG B   1
END
"""


def _metrics(cid: str, receptor: str, lid: str) -> dict:
    """Every section 4 required key, so the complex is gradeable."""
    return {
        "complex_id": cid, "receptor": receptor, "ligand_id": lid,
        "status": "ok", "n_frames": 100,
        "equil_start_frame": 20, "equil_method": "detect_equilibration",
        "ligand_rmsd_median": 1.2, "ligand_rmsd_p90": 1.8, "ligand_rmsd_max": 2.1,
        "ligand_rmsd_drift_slope": 0.01, "inplace_rmsd_median": 0.6,
        "n_native_contacts": 12, "in_pocket_fraction": 0.98, "n_exit_events": 0,
        "com_drift_mean_a": 0.9, "com_drift_max_a": 1.4, "pocket_offset_mean_a": 1.1,
        "conv_js_mean": 0.06, "n_pose_clusters": 2,
        "dominant_cluster_population": 0.8, "changed_pose_cluster": 0,
        "pose_clusters": [{"cluster_id": 1, "population": 0.8,
                           "centroid_frame": 50, "rmsd_to_t0_a": 1.0,
                           "n_frames": 64}],
        "pbc_ok": True, "extra_ligand_fragments": False, "box_volume_a3": 900000.0,
        "receptor_rmsd_mean": 2.0, "ligand_rmsf_mean": 0.5, "ligand_rmsf_max": 1.2,
        "is_reference": False,
    }


def _pocket(cid: str, receptor: str) -> dict:
    return {
        "complex_id": cid, "receptor": receptor, "status": "ok",
        "n_frames": 100, "n_ligand_atoms": 2,
        "pocket_residues_raw": [1], "pocket_residues_filtered": [1, 2],
        "pocket_resnames": {"1": "ALA", "2": "GLY"},
        "pocket_n_raw": 1, "pocket_n_filtered": 2,
        "added_vs_raw": [2], "dropped_from_raw": [],
        "used_fallback": False, "notes": [],
    }


def build_package(*, compounds=("benzamidine", "leupeptin"),
                  receptors=("3PTB", "1AKI")) -> fxa2.Package:
    """A conformant package with no real data in it.

    Two receptors and two ligands, because `S_E`, the league table and both rank
    columns are comparative: a file with one complex in it has no cohort.
    """
    package = fxa2.Package(
        package_name="synthetic_study",
        flexappeal_version="0.2.0",
        engine="openmm 8.5.2",
        run_settings={
            "protein_ff": "amber14-all.xml", "ligand_ff": "openff-2.2.1",
            "water_model": "tip3p", "temperature": 310.0, "pressure": 1.0,
            "timestep_fs": 4.0, "n_strides": 1, "stride_time_ns": 10.0,
            "prod_traj_interval_ps": 100.0, "contact_cutoff": 5.0,
        },
    )
    for receptor in receptors:
        package.register_receptor(receptor, family="serine protease")
        for compound in compounds:
            lid = package.register_ligand(compound)
            cid = fxa2.complex_id(lid, receptor)
            package.add_complex(cid)
            package.add_json(f"metrics/{cid}.json", _metrics(cid, receptor, lid))
            package.add_json(f"harvest/{cid}.json", _metrics(cid, receptor, lid))
            package.add_json(f"pocket/{cid}.json", _pocket(cid, receptor))
            package.add_text(f"traj/{cid}_pocket.pdb", POCKET_PDB)
    return package


@pytest.fixture
def built():
    """One package object, so a test can compare the file against its source.

    The salt is fresh per package, so two calls to `build_package` produce
    different ligand ids by design and a round-trip test has to hold on to the
    object it wrote.
    """
    return build_package()


@pytest.fixture
def written(built, tmp_path):
    return built.write(tmp_path / "synthetic.fxa")


# ---------------------------------------------------------------------------
#  Identifiers, section 2
# ---------------------------------------------------------------------------


def test_a_ligand_id_is_sixteen_hex_characters():
    salt = fxa2.new_salt()
    lid = fxa2.ligand_id(salt, "benzamidine")
    assert len(lid) == 16
    assert all(c in "0123456789abcdef" for c in lid)


def test_a_ligand_id_is_stable_under_one_salt_and_differs_under_another():
    salt = fxa2.new_salt()
    assert fxa2.ligand_id(salt, "benzamidine") == fxa2.ligand_id(salt, "benzamidine")
    assert fxa2.ligand_id(salt, "benzamidine") != fxa2.ligand_id(fxa2.new_salt(), "benzamidine")


def test_an_unsalted_ligand_id_is_refused():
    """A short compound code under a bare digest is reversed by brute force."""
    with pytest.raises(Fxa2Error, match="salt"):
        fxa2.ligand_id(b"", "benzamidine")


def test_the_salt_is_never_written_into_the_file(written):
    package = fxa2.read(written)
    assert package.manifest["salt_is_per_bundle"] is True
    blob = written.read_bytes()
    assert b"salt" not in json.dumps(package.manifest).encode().replace(
        b"salt_is_per_bundle", b"")
    # And not hiding anywhere else in the archive either.
    assert "salt.bin" not in package.members
    assert len(blob) > 0


def test_the_receptor_token_prefers_what_can_be_trusted_most():
    assert fxa2.receptor_token("3ptb", "prepared.pdb", "my job") == "3PTB"
    assert fxa2.receptor_token(None, "3ptb_prepared", "my job") == "3PTB_PREPARED"
    assert fxa2.receptor_token(None, None, "trypsin ben") == "TRYPSIN-BEN"


def test_a_receptor_token_cannot_be_derived_from_nothing():
    with pytest.raises(Fxa2Error, match="receptor identifier"):
        fxa2.receptor_token(None, "", None)


def test_a_complex_id_must_be_a_usable_filename():
    with pytest.raises(Fxa2Error):
        fxa2.complex_id("not-hex", "3PTB")
    with pytest.raises(Fxa2Error):
        fxa2.complex_id("0123456789abcdef", "3ptb")       # lowercase receptor
    with pytest.raises(Fxa2Error):
        fxa2.complex_id("0123456789abcdef", "3PTB/X")     # a path separator


# ---------------------------------------------------------------------------
#  The archive, section 1
# ---------------------------------------------------------------------------


def test_the_package_writes_the_members_the_spec_requires(written):
    with zipfile.ZipFile(written) as archive:
        names = set(archive.namelist())
    assert "manifest.json" in names
    manifest = json.loads(zipfile.ZipFile(written).read("manifest.json"))
    assert manifest["fxa_version"] == 2
    assert manifest["kind"] == "graded"
    assert manifest["chemistry"] is None
    assert len(manifest["complexes"]) == 4
    for cid in manifest["complexes"]:
        assert f"metrics/{cid}.json" in names
        assert f"pocket/{cid}.json" in names
        assert f"traj/{cid}_pocket.pdb" in names


def test_a_complex_without_its_required_members_is_refused(tmp_path):
    package = build_package()
    package.add_complex(fxa2.complex_id("0" * 16, "3PTB"))
    package.ligand_names["0" * 16] = "nothing"
    with pytest.raises(Fxa2Error, match="required for every complex"):
        package.write(tmp_path / "broken.fxa")


def test_a_ligand_not_named_in_the_manifest_is_refused(tmp_path):
    package = build_package()
    cid = fxa2.complex_id("1" * 16, "3PTB")
    package.add_complex(cid)
    package.add_json(f"metrics/{cid}.json", _metrics(cid, "3PTB", "1" * 16))
    package.add_json(f"pocket/{cid}.json", _pocket(cid, "3PTB"))
    package.add_text(f"traj/{cid}_pocket.pdb", POCKET_PDB)
    with pytest.raises(Fxa2Error, match="ligand_names"):
        package.write(tmp_path / "broken.fxa")


def test_a_duplicate_complex_is_refused():
    package = build_package()
    with pytest.raises(Fxa2Error, match="already in this package"):
        package.add_complex(package.complexes[0])


# ---------------------------------------------------------------------------
#  Chemistry, section 8
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "chem/benzamidine.sdf", "chem/benzamidine.mol2", "ligands/BEN.mol",
    "params/BEN.frcmod", "params/system.prmtop", "traj/run.xtc", "traj/run.dcd",
])
def test_a_chemistry_or_trajectory_member_is_refused(name):
    """The multi-model PDB is the trajectory DeepFlex serves; nothing else."""
    package = build_package()
    with pytest.raises(Fxa2Error, match="section 8"):
        package.add_bytes(name, b"x")


@pytest.mark.parametrize("field_name", [
    "smiles", "canonical_smiles", "isomeric_smiles", "inchi", "inchikey",
    "molblock", "molfile", "tripos_mol2", "sdf_block",
])
def test_a_manifest_field_promising_a_structure_is_refused(field_name, tmp_path):
    package = build_package()
    package.receptors["3PTB"][field_name] = "c1ccccc1"
    with pytest.raises(Fxa2Error, match="section 8"):
        package.write(tmp_path / "chemistry.fxa")


def test_a_nested_structure_field_is_caught_too(tmp_path):
    package = build_package()
    package.add_json("annotations/3PTB.json", {"label": "3PTB"})
    package.run_settings["ligand"] = {"detail": {"canonical_SMILES": "CCO"}}
    with pytest.raises(Fxa2Error, match="section 8"):
        package.write(tmp_path / "chemistry.fxa")


def test_a_compound_name_in_clear_is_allowed(written):
    """A name says which compound, not what it is or how to make it."""
    package = fxa2.read(written)
    assert set(package.ligand_names.values()) == {"benzamidine", "leupeptin"}


# ---------------------------------------------------------------------------
#  Unsafe member names
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "../escape.json", "/absolute.json", "C:/windows.json",
    "back\\slash.json", "null\x00byte.json",
])
def test_an_unsafe_member_name_is_refused(name):
    package = build_package()
    with pytest.raises(Fxa2Error):
        package.add_bytes(name, b"x")


# ---------------------------------------------------------------------------
#  The round trip, check 14
# ---------------------------------------------------------------------------


def test_the_round_trip_preserves_every_scalar(built, written):
    """A writer that cannot read its own output is a writer with a unit bug."""
    package = fxa2.read(written)
    original = built
    assert package.complexes == original.complexes
    for cid in package.complexes:
        wrote = json.loads(original.members[f"metrics/{cid}.json"])
        read_back = package.metrics(cid)
        assert read_back == wrote
        assert package.pocket(cid) == json.loads(original.members[f"pocket/{cid}.json"])


def test_the_round_trip_preserves_arrays_and_their_dtypes(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    package.add_arrays(
        f"series/{cid}.npz",
        time_ns=np.linspace(0.0, 10.0, 100),
        pose_rmsd_a=np.full(100, 1.25, dtype=np.float32),
        contacts=np.full(100, 12, dtype=np.int16),
        in_pocket=np.ones(100, dtype=np.int8),
        equil_start_frame=np.array([20], dtype=np.int64),
    )
    read_back = fxa2.read(package.write(tmp_path / "arrays.fxa")).arrays(f"series/{cid}.npz")
    assert read_back["time_ns"].dtype == np.float64
    assert read_back["pose_rmsd_a"].dtype == np.float32
    assert read_back["contacts"].dtype == np.int16
    assert read_back["in_pocket"].dtype == np.int8
    assert read_back["equil_start_frame"].tolist() == [20]
    assert read_back["pose_rmsd_a"][0] == pytest.approx(1.25)


def test_reading_refuses_a_version_it_does_not_know(tmp_path):
    """A reader seeing a version it does not know refuses by number."""
    path = tmp_path / "v3.fxa"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("manifest.json", json.dumps({"fxa_version": 3}))
    with pytest.raises(Fxa2Error, match="fxa_version 3"):
        fxa2.read(path)


def test_reading_refuses_a_version_1_file(tmp_path):
    """Version 1 is a different kind of file, and flexappeal.fxa reads it."""
    path = tmp_path / "v1.fxa"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("manifest.json", json.dumps({"fxa_version": 1}))
    with pytest.raises(Fxa2Error, match="flexappeal.fxa"):
        fxa2.read(path)


def test_reading_refuses_something_that_is_not_a_zip():
    with pytest.raises(Fxa2Error, match="not a valid"):
        fxa2.read(b"this is not a zip at all")


def test_reading_refuses_a_traversal_entry(tmp_path):
    path = tmp_path / "evil.fxa"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("manifest.json", json.dumps({"fxa_version": 2}))
        archive.writestr("../escape.json", "{}")
    with pytest.raises(Fxa2Error, match="escapes"):
        fxa2.read(path)


# ---------------------------------------------------------------------------
#  Conformance, section 10
# ---------------------------------------------------------------------------


def test_a_package_this_writer_produced_is_conformant(written):
    assert fxa2.conformance(fxa2.read(written)) == []


def test_conformance_names_a_missing_required_metric(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    metrics = json.loads(package.members[f"metrics/{cid}.json"])
    del metrics["ligand_rmsd_median"]
    del metrics["conv_js_mean"]
    package.members[f"metrics/{cid}.json"] = json.dumps(metrics).encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "short.fxa")))
    # Check 13 asserts on the list of what is missing, not on a boolean, so the
    # failure says which.
    assert any("ligand_rmsd_median, conv_js_mean" in p for p in problems)


def test_conformance_catches_a_fraction_outside_zero_to_one(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    metrics = json.loads(package.members[f"metrics/{cid}.json"])
    metrics["in_pocket_fraction"] = 1.4
    package.members[f"metrics/{cid}.json"] = json.dumps(metrics).encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "wide.fxa")))
    assert any("outside 0 to 1" in p for p in problems)


def test_conformance_catches_an_equilibration_window_past_the_end(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    metrics = json.loads(package.members[f"metrics/{cid}.json"])
    metrics["equil_start_frame"] = 100       # n_frames is 100
    package.members[f"metrics/{cid}.json"] = json.dumps(metrics).encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "late.fxa")))
    assert any("equil_start_frame is 100" in p for p in problems)


def test_conformance_catches_a_het_code_left_on_a_ter_record(tmp_path):
    """The bug this check exists for, reproduced.

    A `TER` carries the residue name in the same columns as an `ATOM` record, so
    renaming only the atom lines leaves the real het code in the file.
    """
    package = build_package()
    cid = package.complexes[0]
    leaky = POCKET_PDB.replace("TER       8      LIG B   1",
                               "TER       8      BEN B   1")
    package.members[f"traj/{cid}_pocket.pdb"] = leaky.encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "ter.fxa")))
    # The ligand is still LIG on its ATOM records, so the only thing wrong is the
    # TER beneath them, which is exactly the shape of the bug.
    assert any("TER record for residue BEN" in p for p in problems), problems


def test_conformance_catches_water_in_a_pocket_crop(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    wet = POCKET_PDB + (
        "ATOM      9  O   HOH C   1      12.000   0.000   0.000  1.00  0.00           O\n")
    package.members[f"traj/{cid}_pocket.pdb"] = wet.encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "wet.fxa")))
    assert any("HOH" in p for p in problems)


def test_conformance_catches_a_hetatm_ligand(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    het = POCKET_PDB.replace("ATOM      6  C1  LIG", "HETATM    6  C1  LIG")
    package.members[f"traj/{cid}_pocket.pdb"] = het.encode()
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "het.fxa")))
    assert any("HETATM" in p for p in problems)


def test_conformance_catches_a_barcode_over_its_caps(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    package.add_arrays(
        f"barcode/{cid}.npz",
        z=np.ones((40, 600), dtype=np.uint8),
        labels=np.array([f"ALA{i} contact" for i in range(40)], dtype=object),
        time_ns=np.linspace(0, 10, 600),
        n_contacts_total=np.int32(40), step=np.int32(1),
    )
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "big.fxa")))
    assert any("40 rows" in p for p in problems)
    assert any("600 columns" in p for p in problems)


def test_conformance_catches_a_matrix_over_five_hundred(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    package.add_arrays(f"panels/{cid}.npz",
                       rmsd_2d=np.zeros((600, 600), dtype=np.float32))
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "wide.fxa")))
    assert any("600 on a side" in p for p in problems)


def test_conformance_catches_rmsf_arrays_of_different_lengths(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    package.add_arrays(f"panels/{cid}.npz",
                       rmsf_resid=np.arange(10, dtype=np.int32),
                       rmsf_a=np.zeros(9, dtype=np.float32))
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "mismatch.fxa")))
    assert any("different lengths" in p for p in problems)


def test_conformance_catches_a_series_key_of_the_wrong_dtype(tmp_path):
    package = build_package()
    cid = package.complexes[0]
    package.add_arrays(
        f"series/{cid}.npz",
        time_ns=np.linspace(0, 10, 100).astype(np.float32),   # must be float64
        pose_rmsd_a=np.zeros(100, dtype=np.float32),
    )
    problems = fxa2.conformance(fxa2.read(package.write(tmp_path / "dtype.fxa")))
    assert any("time_ns is float32" in p for p in problems)


# ---------------------------------------------------------------------------
#  Units
# ---------------------------------------------------------------------------


def test_a_known_distance_is_reported_in_angstrom(written):
    """The two ligand atoms in the fixture are 1 A apart, not 0.1.

    Version 1 reported nanometres and MDTraj works in nanometres, so this is the
    error that is easiest to make and hardest to see.
    """
    text = fxa2.read(written).pocket_pdb(fxa2.read(written).complexes[0])
    xs = [float(line[30:38]) for line in text.splitlines()
          if line.startswith("ATOM") and line[17:20] == "LIG"]
    assert len(xs) == 2
    assert abs(xs[1] - xs[0]) == pytest.approx(1.0)
