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


# ---------------------------------------------------------------------------
#  The scored window, section 4.1
# ---------------------------------------------------------------------------


def test_equilibration_is_detected_where_it_actually_happens():
    """A series that drifts and then settles should lose roughly the drift.

    Chodera's rule maximises the effectively uncorrelated sample count, so the
    truncation it picks is not exact. What it must not do is keep the drift or
    throw away the plateau.
    """
    rng = np.random.default_rng(1)
    drifting = np.linspace(0.0, 4.0, 60)
    settled = 4.0 + rng.normal(0.0, 0.05, 140)
    start, method = measure.detect_equilibration(np.concatenate([drifting, settled]))
    assert 20 <= start <= 100, start
    assert method == "detect_equilibration"


def test_a_series_that_never_settles_reports_the_ceiling():
    """Capped is honest; silently discarding most of a run is not."""
    start, method = measure.detect_equilibration(np.linspace(0.0, 20.0, 200))
    assert "ceiling" in method
    assert start < 200


def test_a_series_too_short_to_judge_is_scored_whole():
    start, method = measure.detect_equilibration(np.arange(5.0))
    assert (start, method) == (0, "all_frames")


def test_the_scored_window_never_runs_past_the_end(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        metrics = read_back.metrics(cid)
        assert 0 <= metrics["equil_start_frame"] < metrics["n_frames"]


def test_the_statistical_inefficiency_of_noise_is_about_one():
    rng = np.random.default_rng(2)
    assert measure.statistical_inefficiency(rng.normal(size=500)) < 1.5


def test_the_statistical_inefficiency_of_a_slow_drift_is_large():
    assert measure.statistical_inefficiency(np.linspace(0, 1, 500)) > 10


# ---------------------------------------------------------------------------
#  S_pose and S_occ, in angstrom
# ---------------------------------------------------------------------------


def _translated(deposited, offsets_a, *, compound: str = "BEN"):
    """The deposited frame repeated, with the ligand rigidly displaced.

    No jitter at all, so every number out of it is exact and a unit error has
    nowhere to hide: the pose RMSD of a rigid translation is the length of the
    translation, in whatever unit the file claims to be using.
    """
    import mdtraj as md

    frame = deposited[0][:]
    ligand = frame.topology.select(f"resname {compound}")
    frames = []
    for offset in offsets_a:
        xyz = frame.xyz[0].copy()
        xyz[ligand] += np.asarray(offset, dtype=np.float32) / 10.0
        frames.append(xyz)
    return md.Trajectory(
        xyz=np.stack(frames).astype(np.float32), topology=frame.topology,
        time=np.arange(len(frames), dtype=np.float32) * 10.0,
        unitcell_lengths=np.repeat(frame.unitcell_lengths, len(frames), axis=0),
        unitcell_angles=np.repeat(frame.unitcell_angles, len(frames), axis=0),
    )


def _measured(traj, compound: str = "BEN"):
    return measure.measure_complex(
        traj, complex_id=fxa2.complex_id("0" * 16, "3PTB"), receptor="3PTB",
        ligand_id="0" * 16, ligand_resnames=[compound])


def test_a_two_angstrom_shift_is_reported_as_two_not_nought_point_two(deposited):
    """The units test. MDTraj works in nanometres and version 1 reported them.

    Get this wrong and every length in the report is out by a factor of ten,
    which is invisible until somebody checks a number by hand.
    """
    traj = _translated(deposited, [[0, 0, 0]] * 4 + [[2.0, 0, 0]] * 16)
    metrics = _measured(traj).metrics
    assert metrics["ligand_rmsd_median"] == pytest.approx(2.0, abs=0.01)
    assert metrics["ligand_rmsd_max"] == pytest.approx(2.0, abs=0.01)
    # The maximum rather than the mean: the mean depends on how many of the four
    # frames at the origin the scored window happens to keep, and this test is
    # about the unit rather than about the window.
    assert metrics["com_drift_max_a"] == pytest.approx(2.0, abs=0.01)


def test_a_rigid_shift_moves_the_pose_and_not_the_shape(deposited):
    """The two RMSDs separate the ligand moving from the ligand changing shape.

    A rigid translation is all pose and no conformation, so the in-place RMSD,
    which superposes on the ligand itself, has to stay at zero.
    """
    traj = _translated(deposited, [[0, 0, 0]] + [[3.0, 0, 0]] * 19)
    metrics = _measured(traj).metrics
    assert metrics["ligand_rmsd_median"] == pytest.approx(3.0, abs=0.01)
    assert metrics["inplace_rmsd_median"] == pytest.approx(0.0, abs=0.01)


def test_the_centroid_drift_is_measured_from_the_ligand_not_the_pocket(deposited):
    """The pocket centroid sits several angstrom from the ligand centre even in a
    perfectly held pose, so measuring drift against it adds a constant offset and
    the 5 A hard gate then fires on ligands that have not moved."""
    traj = _translated(deposited, [[0, 0, 0]] * 20)
    metrics = _measured(traj).metrics
    assert metrics["com_drift_mean_a"] == pytest.approx(0.0, abs=1e-4)
    assert metrics["com_drift_max_a"] == pytest.approx(0.0, abs=1e-4)
    # The offset is still reported, as a diagnostic, and it is not zero.
    assert metrics["pocket_offset_mean_a"] > 0.5


def test_a_ligand_that_leaves_is_counted_out_and_counted_once(deposited):
    """`in_pocket_fraction` is the scored fraction retaining a quarter of the
    frame-0 contacts, and an exit event is a transition from in to out."""
    traj = _translated(deposited, [[0, 0, 0]] * 10 + [[25.0, 0, 0]] * 10)
    metrics = _measured(traj).metrics
    assert metrics["n_native_contacts"] > 4
    assert metrics["n_exit_events"] == 1
    assert 0.0 < metrics["in_pocket_fraction"] < 1.0


def test_a_ligand_that_stays_scores_full_occupancy(deposited):
    traj = _translated(deposited, [[0, 0, 0]] * 20)
    metrics = _measured(traj).metrics
    assert metrics["in_pocket_fraction"] == 1.0
    assert metrics["n_exit_events"] == 0


def test_a_periodic_jump_is_reported_rather_than_defaulted_true(deposited):
    """`pbc_ok: false` is not a failure, but it makes every distance suspect."""
    traj = _translated(deposited, [[0, 0, 0]] * 5 + [[60.0, 0, 0]] * 5
                       + [[0, 0, 0]] * 5)
    measured = _measured(traj)
    assert measured.metrics["pbc_ok"] is False
    assert any("periodic boundary" in note for note in measured.notes)


def test_a_held_ligand_reports_pbc_ok(deposited):
    traj = _translated(deposited, [[0, 0, 0]] * 20)
    assert _measured(traj).metrics["pbc_ok"] is True


def test_the_box_volume_is_cubic_angstrom(deposited):
    """A 64 A dodecahedral cell is about 190,000 cubic angstrom, not 190."""
    traj = _translated(deposited, [[0, 0, 0]] * 4)
    assert _measured(traj).metrics["box_volume_a3"] > 100_000


def test_the_drift_slope_is_per_nanosecond(deposited):
    """A 1 A rise over the scored window, with the window in nanoseconds.

    The frames are 10 ps apart, so twenty of them span 0.19 ns and a 1 A rise is
    a slope of about 5 A/ns. In picoseconds it would read 0.005, which is the
    error this catches.
    """
    traj = _translated(deposited, [[x, 0, 0] for x in np.linspace(0.0, 1.0, 20)])
    slope = _measured(traj).metrics["ligand_rmsd_drift_slope"]
    assert 2.0 < slope < 12.0, slope


def test_every_scored_quantity_uses_the_window_and_not_the_whole_run(deposited):
    """The discarded head must not reach a scored number.

    The ligand sits in the pocket for the first half and 25 A away for the
    second. Equilibration is therefore detected at the move, so the scored window
    is the half where the ligand is gone and `in_pocket_fraction` is 0. Computed
    over the whole run it would be 0.5, and that is the error this catches: a
    number averaged over frames the file says were discarded.
    """
    traj = _translated(deposited, [[0, 0, 0]] * 10 + [[25.0, 0, 0]] * 10)
    metrics = _measured(traj).metrics
    assert metrics["equil_start_frame"] >= 8, metrics["equil_start_frame"]
    assert metrics["in_pocket_fraction"] == pytest.approx(0.0, abs=0.15)
    # And the window is reported, so a reader can see what was thrown away.
    assert "ceiling" in metrics["equil_method"] or metrics["equil_start_frame"] < 10


# ---------------------------------------------------------------------------
#  series/ and the multi-model trajectory, sections 5.1 and 6.2
# ---------------------------------------------------------------------------


def test_the_series_carries_every_key_at_the_stated_dtype(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        series = read_back.arrays(f"series/{cid}.npz")
        assert series is not None
        for key, dtype in fxa2._SERIES_REQUIRED.items():
            assert key in series, key
            assert series[key].dtype.name == dtype, (key, series[key].dtype)


def test_the_series_covers_the_whole_trajectory_not_the_window(written):
    """DeepFlex draws the discarded head in a different colour, which is how a
    reader sees what was thrown away."""
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        series = read_back.arrays(f"series/{cid}.npz")
        metrics = read_back.metrics(cid)
        assert len(series["time_ns"]) == metrics["n_frames"]
        assert series["equil_start_frame"][0] == metrics["equil_start_frame"]


def test_the_series_and_the_scalars_agree(written):
    """The scalar is the series reduced over the window, so they must match.

    Two codepaths that disagree here would mean the chart and the number beside
    it describe different trajectories.
    """
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        series = read_back.arrays(f"series/{cid}.npz")
        metrics = read_back.metrics(cid)
        window = series["pose_rmsd_a"][metrics["equil_start_frame"]:]
        assert float(np.median(window)) == pytest.approx(
            metrics["ligand_rmsd_median"], rel=1e-5)
        assert float(window.max()) == pytest.approx(metrics["ligand_rmsd_max"], rel=1e-5)


def test_the_native_contact_threshold_is_a_quarter_rounded_up(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        series = read_back.arrays(f"series/{cid}.npz")
        metrics = read_back.metrics(cid)
        expected = max(1, int(np.ceil(metrics["n_native_contacts"] * 0.25)))
        assert series["native_contact_threshold"][0] == expected
        assert np.array_equal(series["in_pocket"],
                              (series["contacts"] >= expected).astype(np.int8))


def test_the_ligand_rmsf_is_per_heavy_atom(written):
    read_back = fxa2.read(written)
    cid = read_back.complexes[0]
    series = read_back.arrays(f"series/{cid}.npz")
    metrics = read_back.metrics(cid)
    assert 0 < len(series["ligand_rmsf_a"]) <= read_back.pocket(cid)["n_ligand_atoms"]
    assert float(series["ligand_rmsf_a"].max()) == pytest.approx(
        metrics["ligand_rmsf_max"], rel=1e-5)


def test_the_trajectory_is_one_multi_model_pdb(written):
    """One file, not one per frame.

    The viewer steps the model index, which is a parameter update; a file per
    frame makes it tear the structure down and rebuild it, so the canvas blanks
    and the camera jumps on every step.
    """
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        text = read_back.members[f"traj/{cid}_scrub.pdb"].decode()
        models = text.count("\nMODEL ") + text.startswith("MODEL ")
        assert models > 1
        assert models <= 110, models
        assert text.count("ENDMDL") == models


def test_the_trajectory_crop_matches_the_pocket_crop(written):
    """Same atom selection as the pocket PDB, so the viewer can swap between them."""
    read_back = fxa2.read(written)
    cid = read_back.complexes[0]
    pocket = read_back.pocket_pdb(cid)
    scrub = read_back.members[f"traj/{cid}_scrub.pdb"].decode()

    def atom_names(text):
        return [line[12:26] for line in text.splitlines() if line.startswith("ATOM")]

    first_model = scrub.split("ENDMDL")[0]
    assert atom_names(first_model) == atom_names(pocket)


def test_the_trajectory_names_its_ligand_lig_in_every_model(written):
    read_back = fxa2.read(written)
    text = read_back.members[f"traj/{read_back.complexes[0]}_scrub.pdb"].decode()
    for line in text.splitlines():
        if line[:6].strip() in ("ATOM", "HETATM", "TER"):
            assert line[17:20].strip() not in ("BEN", "AMD"), line
    assert any(line[17:20].strip() == "LIG" for line in text.splitlines()
               if line.startswith("TER"))


def test_a_trajectory_with_no_time_stamps_falls_back_to_the_frame_spacing(deposited):
    """A time axis of all zeros collapses every chart without anything raising."""
    import mdtraj as md

    traj = _translated(deposited, [[0, 0, 0]] * 12)
    blank = md.Trajectory(
        xyz=traj.xyz, topology=traj.topology, time=np.zeros(traj.n_frames),
        unitcell_lengths=traj.unitcell_lengths, unitcell_angles=traj.unitcell_angles)

    axis, note = measure.time_axis(blank, 20.0)
    assert axis[1] == pytest.approx(0.02)          # 20 ps is 0.02 ns
    assert "frame spacing" in note

    axis, note = measure.time_axis(blank, None)
    assert axis[1] == 1.0
    assert "frame index" in note


def test_a_trajectory_with_its_own_time_stamps_keeps_them(deposited):
    traj = _translated(deposited, [[0, 0, 0]] * 12)      # 10 ps apart
    axis, note = measure.time_axis(traj, 999.0)
    assert note is None
    assert axis[1] == pytest.approx(0.01)


# ---------------------------------------------------------------------------
#  S_int, the occupancy table and the barcode, sections 5.3 and 5.5
# ---------------------------------------------------------------------------


def test_the_fingerprint_finds_the_interaction_benzamidine_actually_makes(written):
    """Benzamidine's amidinium makes a salt bridge to trypsin's ASP189.

    That is the interaction the complex is known for, so it is the one worth
    asserting: a fingerprint that reports plausible occupancies and misses this
    is measuring the wrong atoms.
    """
    read_back = fxa2.read(written)
    rows = read_back.summary(read_back.complexes[0])
    salt_bridge = [row for row in rows
                   if row["resid"] == 189 and row["itype"] == "Ionic"]
    assert salt_bridge, sorted({(r["resid"], r["itype"]) for r in rows})
    assert salt_bridge[0]["occupancy"] > 0.5
    assert salt_bridge[0]["present_at_t0"] == 1


def test_the_occupancy_table_lists_only_interactions_that_happened(written):
    """A row at 0% that was also absent at frame 0 claims a contact never seen.

    The candidate set is every pair the topology makes geometrically possible,
    which for a pocket is mostly pairs that never come within their cutoff.
    Keeping those as rows also drags `hbond_occupancy_mean` down to the mean over
    dozens of zeros rather than over the hydrogen bonds that exist. A contact
    present at frame 0 and then lost is kept, because that loss is exactly what
    the pose memory measures.
    """
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        for row in read_back.summary(cid):
            assert row["occupancy"] > 0.0 or row["present_at_t0"] == 1, row


def test_every_occupancy_lies_in_zero_to_one(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        for row in read_back.summary(cid):
            assert 0.0 <= row["occupancy"] <= 1.0
        record = read_back.interactions(cid)
        for key in ("tanimoto_t0_mean", "pose_memory_mean",
                    "key_interaction_occupancy_mean", "hbond_occupancy_mean"):
            if record.get(key) is not None:
                assert 0.0 <= record[key] <= 1.0, key


def test_the_interaction_record_carries_its_status(written):
    """DeepFlex's scorer skips the sub-score without it, and the spec omits it."""
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        record = read_back.interactions(cid)
        assert record["status"] == "ok"
        assert record["complex_id"] == cid


def test_the_record_says_the_typing_is_geometric(written):
    """The type vocabulary is coarser than ProLIF's and the file has to say so."""
    read_back = fxa2.read(written)
    notes = " ".join(read_back.interactions(read_back.complexes[0])["notes"])
    assert "pi-stack" in notes


def test_the_barcode_is_within_its_caps_and_labelled(written):
    read_back = fxa2.read(written)
    for cid in read_back.complexes:
        barcode = read_back.arrays(f"barcode/{cid}.npz")
        rows, cols = barcode["z"].shape
        assert rows <= fxa2.MAX_BARCODE_ROWS
        assert cols <= fxa2.MAX_BARCODE_COLS
        assert len(barcode["labels"]) == rows
        assert barcode["z"].dtype == np.uint8
        assert set(np.unique(barcode["z"])) <= {0, 1}
        assert len(barcode["time_ns"]) == cols
        assert int(barcode["n_contacts_total"]) >= rows


def test_the_barcode_is_busiest_first(written):
    """Twenty-eight rows is already a tall panel and the tail is noise."""
    read_back = fxa2.read(written)
    occupancies = read_back.arrays(
        f"barcode/{read_back.complexes[0]}.npz")["z"].mean(axis=1)
    assert list(occupancies) == sorted(occupancies, reverse=True)


def test_the_barcode_labels_read_resname_resid_type(written):
    read_back = fxa2.read(written)
    for label in read_back.arrays(f"barcode/{read_back.complexes[0]}.npz")["labels"]:
        residue, itype = str(label).split(" ")
        assert itype in ("HBond", "Ionic", "Hydrophobic")
        assert residue[:3].isalpha() and residue[3:].isdigit()


def test_a_ligand_that_leaves_loses_its_pose_memory(deposited):
    """The fingerprint has to register the ligand going, not just the RMSD."""
    held = _measured(_translated(deposited, [[0, 0, 0]] * 20)).interactions
    gone = _measured(
        _translated(deposited, [[0, 0, 0]] + [[25.0, 0, 0]] * 19)).interactions
    assert held["pose_memory_mean"] == pytest.approx(1.0)
    assert gone["pose_memory_mean"] < 0.1
    assert gone["tanimoto_t0_mean"] < 0.1


def test_the_mean_lifetime_separates_one_contact_from_many_flickers(deposited):
    """Two series at the same occupancy, one held and one flickering."""
    held = np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0], dtype=bool)
    flicker = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0], dtype=bool)
    assert held.mean() < flicker.mean() + 0.2        # roughly the same occupancy
    assert measure._mean_lifetime(held, 0.01) > measure._mean_lifetime(flicker, 0.01)


def test_a_departed_ligand_keeps_no_interaction_it_has_lost(deposited):
    """Omit the key; never write a zero. A zero is read as a measurement."""
    traj = _translated(deposited, [[0, 0, 0]] + [[25.0, 0, 0]] * 19)
    measured = _measured(traj)
    assert measured.interactions is not None
    assert measured.interactions["pose_memory_mean"] < 0.1
    assert measured.interactions["n_t0_interactions"] > 0


def test_a_shift_the_size_of_the_box_is_a_periodic_image_and_is_flagged(deposited):
    """Why `pbc_ok` exists at all, made concrete.

    The cell here is about 64 A across, so a 60 A displacement is the same point
    under the minimum image convention and the contacts correctly persist. The
    measurements are not wrong: the coordinates are, and the only honest thing
    the file can do is say so.
    """
    traj = _translated(deposited, [[0, 0, 0]] + [[60.0, 0, 0]] * 11)
    measured = _measured(traj)
    assert measured.metrics["pbc_ok"] is False
    assert any("periodic boundary" in note for note in measured.notes)
    # And the contacts really do survive, which is the point.
    assert measured.metrics["in_pocket_fraction"] > 0.5
