"""`flexappeal package`: collecting finished runs into one `.fxa` version 2.

Every trajectory here is built in `tmp_path` from `tests/fixtures/3ptb.pdb`, which
is bovine trypsin with benzamidine in the S1 pocket: a public structure, and the
same one the existing fixtures use. Nothing is committed, because a `.fxa` carries
a topology and coordinates and this repository is public.

The trajectory is synthetic: one deposited frame, jittered, with a small rigid
drift applied to the ligand so the pose numbers have something to measure. That is
enough to test the shapes, the units and the member layout, which is what this
module is responsible for. It is not a simulation and no number out of it means
anything physically.
"""

from __future__ import annotations

import json
import warnings
import zipfile

import numpy as np
import pytest

from conftest import FIXTURES
from flexappeal import fxa2, measure, package
from flexappeal.package import PackageError

TRYPSIN = FIXTURES / "3ptb.pdb"


@pytest.fixture(scope="module")
def deposited():
    """Trypsin and benzamidine, heavy atoms, no waters or ions."""
    md = pytest.importorskip("mdtraj")
    if not TRYPSIN.is_file():
        pytest.skip("tests/fixtures/3ptb.pdb is missing")
    with warnings.catch_warnings():
        # 3PTB uses chymotrypsin numbering with insertion codes, which MDTraj
        # drops, so it warns about consecutive residues sharing a number.
        warnings.simplefilter("ignore", UserWarning)
        frame = md.load(str(TRYPSIN))
    return frame.atom_slice(frame.top.select("protein or resname BEN"))


def make_trajectory(deposited, *, compound: str = "BEN", n_frames: int = 60,
                    jitter_nm: float = 0.015, drift_a: float = 1.0, seed: int = 0):
    """A synthetic trajectory, with the ligand renamed and given a slow drift."""
    import mdtraj as md

    frame = deposited[0][:]
    topology = frame.topology
    if compound != "BEN":
        for residue in topology.residues:
            if residue.name == "BEN":
                residue.name = compound

    ligand = topology.select(f"resname {compound}")
    rng = np.random.default_rng(seed)
    frames = []
    for index in range(n_frames):
        xyz = frame.xyz[0] + rng.normal(0.0, jitter_nm, frame.xyz[0].shape)
        # A linear drift on the ligand alone, so the pose RMSD rises and the
        # in-place RMSD does not: that is the distinction those two measure.
        xyz[ligand] += np.array([drift_a / 10.0 * index / max(1, n_frames - 1), 0, 0])
        frames.append(xyz)

    return md.Trajectory(
        xyz=np.stack(frames).astype(np.float32), topology=topology,
        time=np.arange(n_frames, dtype=np.float32) * 10.0,
        unitcell_lengths=np.repeat(frame.unitcell_lengths, n_frames, axis=0),
        unitcell_angles=np.repeat(frame.unitcell_angles, n_frames, axis=0),
    )


def make_source(deposited, *, compound: str = "BEN", pdb_id: str = "3PTB",
                job: str = "trypsin", **kwargs) -> package.Source:
    return package.Source(
        label=f"{compound}_{pdb_id}",
        origin="run directory",
        trajectory=make_trajectory(deposited, compound=compound, **kwargs),
        config={
            "job_name": job, "pdb_id": pdb_id, "keep_heteroatoms": [compound],
            "protein_ff": "amber14-all.xml", "water_model": "tip3p",
            "temperature": 310.0, "pressure": 1.0, "timestep": 4.0,
            "production_duration": 10.0, "traj_interval": 10.0,
        },
        structure={"filename": f"{pdb_id.lower()}.pdb"},
        run_manifest={"openmm_version": "8.5.2"},
        frame_spacing_ps=10.0,
    )


@pytest.fixture
def study(deposited):
    """Two receptors by two ligands: the smallest package with a real cohort."""
    return [
        make_source(deposited, compound=compound, pdb_id=receptor, seed=seed)
        for seed, (compound, receptor) in enumerate(
            [(c, r) for r in ("3PTB", "1TPO") for c in ("BEN", "AMD")])
    ]


# ---------------------------------------------------------------------------
#  Identifying a receptor and a ligand across jobs
# ---------------------------------------------------------------------------


def test_the_receptor_comes_from_the_accession_first(deposited):
    """`pdb_id` is what the user asked for, so it is trusted over a filename."""
    source = make_source(deposited, pdb_id="3PTB")
    source.structure["filename"] = "something_else.pdb"
    assert source.receptor() == "3PTB"


def test_without_an_accession_the_prepared_filename_is_used(deposited):
    source = make_source(deposited)
    source.config["pdb_id"] = ""
    source.structure["filename"] = "trypsin_prepared.pdb"
    assert source.receptor() == "TRYPSIN_PREPARED"


def test_without_either_the_job_name_is_used(deposited):
    source = make_source(deposited, job="my trypsin run")
    source.config["pdb_id"] = ""
    source.structure = {}
    assert source.receptor() == "MY-TRYPSIN-RUN"


def test_two_runs_of_the_same_pair_cannot_share_a_file(deposited):
    """A complex is one receptor and one ligand, so a repeat is an error.

    Silently keeping one would drop a run the user asked to include, and
    silently overwriting would put one run's numbers under the other's name.
    """
    sources = [make_source(deposited, pdb_id="3PTB", job="first"),
               make_source(deposited, pdb_id="3PTB", job="second")]
    with pytest.raises(PackageError, match="both describe"):
        package.assemble(sources, package_name="duplicate")


def test_an_apo_run_is_skipped_with_a_reason(deposited):
    apo = make_source(deposited)
    apo.config["keep_heteroatoms"] = []
    holo = make_source(deposited, compound="AMD", pdb_id="1TPO")
    built, warnings_out = package.assemble([apo, holo], package_name="mixed")
    assert len(built.complexes) == 1
    assert any("no ligand" in w for w in warnings_out)


def test_a_package_of_apo_runs_only_is_refused(deposited):
    apo = make_source(deposited)
    apo.config["keep_heteroatoms"] = []
    with pytest.raises(PackageError, match="package of receptor-ligand complexes"):
        package.assemble([apo], package_name="apo")


def test_one_complex_is_packaged_but_says_what_it_cannot_show(deposited):
    """A cohort of one is no cohort, and the warning has to say so."""
    built, warnings_out = package.assemble(
        [make_source(deposited)], package_name="single")
    assert len(built.complexes) == 1
    assert any("no cohort" in w for w in warnings_out)


# ---------------------------------------------------------------------------
#  The members, section 1
# ---------------------------------------------------------------------------


@pytest.fixture
def written(study, tmp_path):
    built, _ = package.assemble(study, package_name="trypsin_study")
    return built.write(tmp_path / "study.fxa")


def test_the_package_carries_every_required_member(written):
    read_back = fxa2.read(written)
    assert len(read_back.complexes) == 4
    assert sorted(read_back.manifest["receptors"]) == ["1TPO", "3PTB"]
    for cid in read_back.complexes:
        assert read_back.metrics(cid) is not None
        assert read_back.pocket(cid) is not None
        assert read_back.pocket_pdb(cid) is not None


def test_the_package_is_conformant(written):
    assert fxa2.conformance(fxa2.read(written)) == []


def test_the_manifest_reports_the_run_settings_by_name(written):
    settings = fxa2.read(written).run_settings
    assert settings["protein_ff"] == "amber14-all.xml"
    assert settings["temperature"] == 310.0
    assert settings["n_strides"] == 1
    # The production length and the frame spacing: the report prints both and uses
    # them to turn a frame index into nanoseconds.
    assert settings["stride_time_ns"] == 10.0
    assert settings["prod_traj_interval_ps"] == 10.0
    assert settings["contact_cutoff"] == measure.CONTACT_CUTOFF_A


def test_compound_names_travel_in_clear_and_ids_are_hashed(written):
    read_back = fxa2.read(written)
    assert sorted(read_back.ligand_names.values()) == ["AMD", "BEN"]
    for lid in read_back.ligand_names:
        assert len(lid) == 16 and all(c in "0123456789abcdef" for c in lid)


def test_no_member_carries_chemistry(written):
    with zipfile.ZipFile(written) as archive:
        names = archive.namelist()
    for name in names:
        suffix = "." + name.rsplit(".", 1)[-1].lower()
        assert suffix not in fxa2.CHEMISTRY_SUFFIXES, name


def test_sources_that_disagree_on_a_setting_are_flagged(deposited):
    """One manifest describes a whole package, so a disagreement has to be said."""
    cold = make_source(deposited, pdb_id="3PTB")
    hot = make_source(deposited, compound="AMD", pdb_id="1TPO")
    hot.config["temperature"] = 300.0
    _, warnings_out = package.assemble([cold, hot], package_name="mixed")
    assert any("disagree on temperature" in w for w in warnings_out)


# ---------------------------------------------------------------------------
#  The pocket, sections 6.1 and 6.5
# ---------------------------------------------------------------------------


def test_the_pocket_is_the_trypsin_s1_site(written):
    """A check that the pocket is the real one, not merely the right shape.

    Benzamidine binds trypsin's S1 pocket, whose floor is ASP189 and whose walls
    include SER190, GLY216 and TRP215. If the pocket selection ever drifts, these
    are the residues that stop appearing.
    """
    read_back = fxa2.read(written)
    pocket = read_back.pocket(read_back.complexes[0])
    filtered = set(pocket["pocket_residues_filtered"])
    assert {189, 190, 215, 216} <= filtered, sorted(filtered)
    assert pocket["pocket_resnames"]["189"] == "ASP"
    assert pocket["pocket_resnames"]["190"] == "SER"


def test_the_filtered_pocket_is_a_subset_of_the_crop_shell(written):
    """The filter drops; it never invents.

    The filtered set is what the ligand RMSD is superposed on, so it is the
    narrower 6 A selection with the chain termini removed. The raw set is the
    10 A crop shell the pocket PDB contains.
    """
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        pocket = read_back.pocket(cid)
        assert pocket["added_vs_raw"] == []
        assert pocket["dropped_from_raw"]
        assert pocket["pocket_n_filtered"] < pocket["pocket_n_raw"]
        assert pocket["pocket_n_filtered"] == len(pocket["pocket_residues_filtered"])


def test_the_filtered_pocket_excludes_the_chain_termini(written):
    """A terminus makes the superposition wander, and then the ligand looks mobile."""
    read_back = fxa2.read(written)
    pocket = read_back.pocket(read_back.complexes[0])
    assert 229 in pocket["pocket_residues_raw"]       # trypsin's last residue
    assert 229 not in pocket["pocket_residues_filtered"]


def test_the_pocket_crop_names_its_ligand_lig_everywhere(written):
    read_back = fxa2.read(written)
    text = read_back.pocket_pdb(read_back.complexes[0])
    assert "LIG" in text
    for line in text.splitlines():
        if line[:6].strip() in ("ATOM", "HETATM", "TER"):
            assert line[17:20].strip() not in ("BEN", "AMD"), line


def test_the_pocket_crop_renames_the_ter_record_too(written):
    """The bug this exists for: a `TER` carries the residue name in the same
    columns as an `ATOM` record, so renaming only the atom lines leaves the real
    het code in the file under a hashed filename."""
    read_back = fxa2.read(written)
    text = read_back.pocket_pdb(read_back.complexes[0])
    ter = [line for line in text.splitlines() if line.startswith("TER")]
    assert ter, "mdtraj should write a TER for each chain"
    assert any(line[17:20].strip() == "LIG" for line in ter), ter


def test_the_pocket_crop_has_no_solvent_or_ions(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        text = read_back.pocket_pdb(cid)
        names = {line[17:20].strip() for line in text.splitlines()
                 if line[:6].strip() in ("ATOM", "HETATM")}
        assert not names & (fxa2.SOLVENT_RESNAMES | fxa2.ION_RESNAMES)


def test_a_waters_only_source_still_crops_dry(deposited):
    """A run directory's topology has explicit solvent in it.

    The version 1 payload tiers have usually dropped it, but a run directory has
    not, and a 10 A shell drawn round a ligand in explicit solvent is mostly
    water.
    """
    import mdtraj as md

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        wet = md.load(str(TRYPSIN))                   # waters included
    source = package.Source(
        label="wet", origin="run directory",
        trajectory=md.Trajectory(
            xyz=np.repeat(wet.xyz, 4, axis=0), topology=wet.topology,
            time=np.arange(4, dtype=np.float32) * 10.0,
            unitcell_lengths=np.repeat(wet.unitcell_lengths, 4, axis=0),
            unitcell_angles=np.repeat(wet.unitcell_angles, 4, axis=0)),
        config={"job_name": "wet", "pdb_id": "3PTB", "keep_heteroatoms": ["BEN"]},
        structure={"filename": "3ptb.pdb"}, frame_spacing_ps=10.0,
    )
    built, _ = package.assemble([source], package_name="wet")
    cid = built.complexes[0]
    text = built.members[f"traj/{cid}_pocket.pdb"].decode()
    assert "HOH" not in text


def test_the_pocket_crop_carries_the_last_frame_as_well(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        assert f"traj/{cid}_pocket_last.pdb" in read_back.members


# ---------------------------------------------------------------------------
#  Honesty about what is and is not there
# ---------------------------------------------------------------------------


def test_a_package_missing_a_required_input_declares_itself_flexibility(written):
    """Section 9: a graded file missing required inputs is refused, and
    `flexibility` always works and is never wrong, only less."""
    read_back = fxa2.read(written)
    metrics = read_back.metrics(read_back.complexes[0])
    kind = read_back.manifest["kind"]
    missing = [k for k in fxa2.REQUIRED_METRIC_KEYS if metrics.get(k) is None]
    assert (kind == "graded") == (not missing), (kind, missing)


def test_a_complex_that_cannot_be_measured_carries_no_grade(deposited):
    """A missing result and a poor result are different findings."""
    broken = make_source(deposited)
    broken.config["keep_heteroatoms"] = ["ZZZ"]       # no such residue
    good = make_source(deposited, compound="AMD", pdb_id="1TPO")
    built, warnings_out = package.assemble([broken, good], package_name="partial")
    failures = [cid for cid in built.complexes
                if json.loads(built.members[f"metrics/{cid}.json"])["status"]
                != "ok"]
    assert len(failures) == 1
    record = json.loads(built.members[f"metrics/{failures[0]}.json"])
    assert record["status"] == "failed:no_ligand"
    assert "ligand_rmsd_median" not in record
    assert any("different findings" in w for w in warnings_out)


def test_no_measurement_is_written_as_zero(written):
    """A zero is read as a measurement; an absent key is read as pending."""
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        metrics = read_back.metrics(cid)
        for key, value in metrics.items():
            assert value != -1, key
        # Nothing the writer cannot produce may be present at all.
        assert "mmpbsa_gb_dg" not in metrics


# ---------------------------------------------------------------------------
#  Supplied energies
# ---------------------------------------------------------------------------


def test_supplied_energies_are_carried_in(deposited, tmp_path):
    """FlexAppeal cannot compute these, so a study that has them can supply them."""
    path = tmp_path / "energies.json"
    path.write_text(json.dumps({
        "BEN": {"mmpbsa_gb_dg": -42.0, "mmpbsa_gb_sd": 2.1, "lie_drift_slope": 0.0},
        "AMD/1TPO": {"mmpbsa_gb_dg": -31.5},
    }))
    sources = [make_source(deposited, compound="BEN", pdb_id="3PTB"),
               make_source(deposited, compound="AMD", pdb_id="1TPO")]
    built, _ = package.assemble(
        sources, package_name="energies", energies=package._read_energies(path))
    values = {}
    for cid in built.complexes:
        record = json.loads(built.members[f"metrics/{cid}.json"])
        values[record["receptor"]] = record.get("mmpbsa_gb_dg")
    assert values == {"3PTB": -42.0, "1TPO": -31.5}


def test_an_absolute_mmpbsa_block_is_refused(deposited, tmp_path):
    """The three absolute blocks are each around -100,000 kcal/mol.

    Only their difference is a binding quantity, and a percentile over absolute
    blocks ranks the size of the system rather than the affinity of the ligand.
    """
    path = tmp_path / "energies.json"
    path.write_text(json.dumps({"BEN": {"mmpbsa_gb_dg": -102345.6}}))
    with pytest.raises(PackageError, match="difference"):
        package.assemble([make_source(deposited)], package_name="absolute",
                         energies=package._read_energies(path))


def test_a_supplied_file_cannot_smuggle_in_another_field(deposited, tmp_path):
    path = tmp_path / "energies.json"
    path.write_text(json.dumps({"BEN": {"mmpbsa_gb_dg": -42.0, "grade": "A",
                                        "score": 99.0}}))
    built, _ = package.assemble(
        [make_source(deposited)], package_name="smuggle",
        energies=package._read_energies(path))
    record = json.loads(built.members[f"metrics/{built.complexes[0]}.json"])
    # The grade is DeepFlex's: the weights and the bands are versioned there, and
    # a second copy of a judgement drifts away from the first.
    assert "grade" not in record
    assert "score" not in record
    assert record["mmpbsa_gb_dg"] == -42.0


# ---------------------------------------------------------------------------
#  Source discovery, and version 1 as a source
# ---------------------------------------------------------------------------


def test_a_version_1_fxa_is_a_usable_source():
    """The existing trypsin fixture is a real run, so it is the honest test."""
    fixture = FIXTURES / "trypsin_ben.fxa"
    if not fixture.is_file():
        pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")
    source = package.load_fxa1(fixture)
    assert source.trajectory.n_frames > 1
    assert source.ligand_resnames == ["BEN"]
    assert source.receptor() == "3PTB"            # from the structure filename
    assert any("payload tier" in note for note in source.notes)


def test_discovery_expands_a_directory_of_fxa_files(tmp_path):
    fixture = FIXTURES / "trypsin_ben.fxa"
    if not fixture.is_file():
        pytest.skip("tests/fixtures/trypsin_ben.fxa is missing")
    (tmp_path / "one.fxa").write_bytes(fixture.read_bytes())
    (tmp_path / "two.fxa").write_bytes(fixture.read_bytes())
    assert len(package.discover([tmp_path])) == 2


def test_discovery_finds_run_directories_under_an_output_tree(tmp_path):
    run = tmp_path / "flexappeal_output" / "job"
    run.mkdir(parents=True)
    (run / "topology.pdb").write_text("END\n")
    (run / "trajectory.xtc").write_bytes(b"")
    assert package.discover([tmp_path]) == [run]


def test_discovery_says_so_when_a_directory_holds_nothing(tmp_path):
    with pytest.raises(PackageError, match="no finished runs"):
        package.discover([tmp_path])


def test_discovery_refuses_a_path_that_is_not_there(tmp_path):
    with pytest.raises(PackageError, match="does not exist"):
        package.discover([tmp_path / "nowhere"])
