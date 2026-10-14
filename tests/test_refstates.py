"""Tests for the reference-states analysis.

Every fixture here is built in the test. Nothing is downloaded and nothing is
committed, which is the rule the rest of this suite follows: a test that needs
the network is a test that fails when the network does, and a committed
structure is a structure somebody has to check before it is published.

The cases chosen are the ones where this build actually got the answer wrong.
Each bug in FINDINGS section 3 produced plausible numbers rather than an error,
so what is being protected is not "does it run" but "does it still refuse the
specific wrong answer it used to give".
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

from flexappeal.refstates import ensembles, membrane, sse, state, store, traj


# ===========================================================================
#  The state score
# ===========================================================================


def test_the_same_word_means_opposite_states_in_different_metrics():
    """`broken` is inactive for the R-spine and active for the ionic lock.

    A single table keyed on the word alone labelled every broken lock inactive
    and then agreed with the geometry for the wrong reason.
    """
    assert state.rule_state("r_spine", "broken") == "inactive"
    assert state.rule_state("ionic_lock", "broken") == "active"


def test_an_unknown_metric_or_label_does_not_guess():
    assert state.rule_state("not_a_metric", "broken") is None
    assert state.rule_state("r_spine", "sideways") is None
    assert state.rule_state("r_spine", "unresolved") is None


def test_geometry_and_rule_disagreeing_reports_conflict():
    motif = state.score_motif(
        "xDFG", active={"2F4J": 1.0}, inactive={"3PYY": 4.0}, n_residues=4,
        rule_label="broken", rule_metric="r_spine")
    assert motif.geometric_label == "active"
    assert motif.label == "conflict"


def test_geometry_and_rule_agreeing_reports_the_state():
    motif = state.score_motif(
        "TM6", active={"3SN6": 1.0}, inactive={"2RH1": 4.0}, n_residues=38,
        rule_label="broken", rule_metric="ionic_lock")
    assert motif.label == "active"


def test_a_motif_with_no_rule_falls_back_to_geometry_and_says_so():
    motif = state.score_motif("hinge", active={"a": 1.0}, inactive={"b": 4.0},
                              n_residues=3)
    assert motif.label == "active"
    assert "rule_label" not in motif.as_row()


def test_an_unmeasured_motif_reports_nothing_rather_than_zero():
    motif = state.score_motif("absent", active={}, inactive={}, n_residues=0)
    assert motif.s is None
    row = motif.as_row()
    for key in ("s", "d_act_a", "d_inact_a", "label"):
        assert key not in row, f"{key} was written for an unmeasured motif"


def test_the_whole_score_is_weighted_by_residue_count():
    """A twenty-residue loop is more of the protein than one tryptophan."""
    frame = state.FrameState(time_ns=1.0, motifs=[
        state.score_motif("big", active={"a": 1.0}, inactive={"b": 3.0},
                          n_residues=20),
        state.score_motif("small", active={"a": 3.0}, inactive={"b": 1.0},
                          n_residues=1),
    ])
    # Unweighted these cancel exactly; weighted the large motif decides.
    assert frame.whole_s is not None
    assert frame.whole_s > 0.3
    assert frame.whole_label == "active"


def test_smoothing_leaves_gaps_as_gaps():
    """An interpolated point on a state trace reads as a measurement."""
    values = [0.1, 0.2, None, 0.4, 0.5]
    out = state.smooth(values, frame_ps=100.0, window_ns=0.2)
    assert out[2] is None
    assert all(v is not None for i, v in enumerate(out) if i != 2)


def test_transitions_ignore_flickers_but_keep_real_changes():
    labels = (["active"] * 20 + ["inactive"] + ["active"] * 20
              + ["inactive"] * 20)
    times = [i * 0.1 for i in range(len(labels))]
    events = state.transitions(labels, times, min_dwell_ns=1.0)
    assert len(events) == 1, f"expected one real transition, got {events}"
    assert events[0]["to"] == "inactive"


# ===========================================================================
#  Superposition
# ===========================================================================


def test_superposition_recovers_a_known_rotation_and_shift():
    rng = np.random.RandomState(0)
    original = rng.rand(30, 3) * 10.0
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = original @ rotation.T + np.array([5.0, 1.0, 2.0])
    transform = traj.superpose(original, moved)
    assert traj.rmsd(traj.apply_transform(original, *transform), moved) < 1e-8


def test_superposition_refuses_to_fit_a_mirror_image():
    """An improper rotation would report a small RMSD for the wrong handedness."""
    rng = np.random.RandomState(1)
    original = rng.rand(30, 3) * 10.0
    mirrored = original * np.array([1.0, 1.0, -1.0])
    transform = traj.superpose(original, mirrored)
    fitted = traj.apply_transform(original, *transform)
    assert traj.rmsd(fitted, mirrored) > 1.0
    assert np.linalg.det(transform[0]) > 0


def test_references_in_unrelated_frames_are_brought_into_one():
    """The bug that read a confident 0.877 "active" on a nonsense comparison."""
    rng = np.random.RandomState(2)
    core = rng.rand(20, 3) * 10.0
    motif = rng.rand(5, 3) * 10.0

    def rotated(seed):
        angle = 0.7 * seed
        spin = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle), np.cos(angle), 0.0],
                         [0.0, 0.0, 1.0]])
        shift = np.array([seed * 30.0, 0.0, 0.0])
        return core @ spin.T + shift, motif @ spin.T + shift

    built = {}
    for which, seeds in (("active", (0, 1)), ("inactive", (2, 3))):
        found = ensembles.MotifEnsemble(state=which)
        for seed in seeds:
            core_i, motif_i = rotated(seed)
            found.by_code[f"R{seed}"] = {ensembles.CORE_KEY: core_i,
                                         "m": motif_i}
        built[which] = found

    # Before: the same motif sits tens of angstrom from itself across frames.
    spread = [traj.rmsd(built[w].by_code[c]["m"], built["active"].by_code["R0"]["m"])
              for w in built for c in built[w].by_code]
    assert max(spread) > 10.0

    report = ensembles.align_to_common_frame(built)
    assert report["n_aligned"] == 4
    after = [traj.rmsd(built[w].by_code[c]["m"], built["active"].by_code["R0"]["m"])
             for w in built for c in built[w].by_code]
    assert max(after) < 1e-6, "identical motifs should coincide once aligned"


# ===========================================================================
#  Secondary structure
# ===========================================================================


def test_an_unassigned_residue_does_not_shift_the_motif():
    """mdtraj writes "NA", two characters, and a joined string misaligns."""
    row = ["H", "H", sse.UNASSIGNED, "H", "E", "E", " "]
    motif = sse.motif_sse(row, "test", list(range(7)), list(range(100, 107)))
    assert len(motif.string) == 7
    assert motif.labels[2] == sse.UNASSIGNED


def test_unassigned_residues_are_not_counted_as_coil():
    helix, strand = sse.fractions(["H", "H", "H", sse.UNASSIGNED])
    assert helix == pytest.approx(1.0)
    assert strand == pytest.approx(0.0)


def test_sse_agreement_omits_residues_with_no_data():
    frame = ["H", "E", sse.UNASSIGNED]
    active = sse.reference_profile({"a": ["H", "H", "H"], "b": ["H", "E", "H"]})
    inactive = sse.reference_profile({"c": ["E", "E", "E"]})
    rows = sse.agreement(frame, active, inactive)
    assert [r["index"] for r in rows] == [0, 1]
    assert rows[0]["active_fraction"] == pytest.approx(1.0)
    assert rows[0]["inactive_fraction"] == pytest.approx(0.0)


# ===========================================================================
#  Membrane
# ===========================================================================


def test_leaflets_split_at_their_own_midplane_not_at_zero():
    """A bilayer drifts in z, and a split at the origin halves the area."""
    drifted = [40.0, 41.0, 42.0, 60.0, 61.0, 62.0]
    upper, lower = membrane.split_leaflets(drifted)
    assert len(upper) == 3 and len(lower) == 3


def test_area_per_lipid_verdicts_follow_the_spec_windows():
    assert membrane.MembraneChunk(0, 10, apl_a2=63.0).apl_verdict == "expected"
    assert membrane.MembraneChunk(0, 10, apl_a2=57.0).apl_verdict == "outside expected"
    assert membrane.MembraneChunk(0, 10, apl_a2=40.0).apl_verdict == "G6"


def test_gate_six_tolerates_one_bad_chunk_but_not_a_bad_run():
    good = [membrane.MembraneChunk(i, i + 10, apl_a2=64.0) for i in range(0, 40, 10)]
    assert membrane.gate_six(good + [membrane.MembraneChunk(40, 50, apl_a2=40.0)])["open"] is False
    assert membrane.gate_six([membrane.MembraneChunk(0, 10, apl_a2=40.0)])["open"] is True


def test_axis_tilt_is_folded_into_the_first_quadrant():
    along_z = np.array([[0.0, 0.0, float(i)] for i in range(10)])
    assert membrane.axis_tilt_deg(along_z) == pytest.approx(0.0, abs=1e-6)
    flipped = np.array([[0.0, 0.0, -float(i)] for i in range(10)])
    assert membrane.axis_tilt_deg(flipped) == pytest.approx(0.0, abs=1e-6)


# ===========================================================================
#  The store
# ===========================================================================


def _analysis_payload() -> dict:
    return {
        "name": "fixture", "kind": "kinase",
        "core": {"n_residues": 43, "source": "fixture core"},
        "motif_residues": {"xDFG": [380, 381, 382, 383], "hinge": [316, 317, 318]},
        "ensembles": {"active": {"codes": ["2F4J"]},
                      "inactive": {"codes": ["3PYY"]}},
        "frames": [
            {"time_ns": 0.0, "n_conflicts": 0, "motifs": [
                {"motif": "xDFG", "d_act_a": 1.0, "d_inact_a": 4.0, "s": 0.6,
                 "label": "active"},
                {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0, "s": 0.0,
                 "label": "intermediate"}]},
            {"time_ns": 1.0, "n_conflicts": 0, "motifs": [
                {"motif": "xDFG", "d_act_a": 0.5, "d_inact_a": 4.5, "s": 0.8,
                 "label": "active"},
                {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0, "s": 0.0,
                 "label": "intermediate"}]},
        ],
        "notes": [],
    }


def test_the_schema_executes(tmp_path: Path):
    """Section 7 names a table `references`, which is a SQL keyword."""
    store.create(tmp_path / "t.db")
    with sqlite3.connect(tmp_path / "t.db") as connection:
        names = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'")}
    assert "reference_structures" in names
    assert "references" not in names


def test_ingesting_twice_does_not_double_the_rows(tmp_path: Path):
    path = tmp_path / "t.db"
    payload = _analysis_payload()
    first = store.ingest_analysis(path, payload, system_id="s", run_id="r",
                                  family="kinase", environment="water")
    store.ingest_analysis(path, payload, system_id="s", run_id="r",
                          family="kinase", environment="water")
    with sqlite3.connect(path) as connection:
        rows = connection.execute("SELECT COUNT(*) FROM motif_frames").fetchone()[0]
    assert rows == first["frames"]


def test_an_in_flight_run_keeps_its_own_status(tmp_path: Path):
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water",
                          run_fields={"status": "running"})
    with sqlite3.connect(path) as connection:
        status = connection.execute(
            "SELECT status FROM runs WHERE run_id = 'r'").fetchone()[0]
    assert status == "running"


def test_membrane_rows_require_their_run(tmp_path: Path):
    """The foreign keys are a constraint here, not decoration.

    SQLite leaves foreign keys off by default, which would have made the
    REFERENCES clauses in section 7's schema purely documentary. With them on, a
    membrane row naming a run that does not exist is refused.
    """
    path = tmp_path / "t.db"
    store.create(path)
    with pytest.raises(sqlite3.IntegrityError):
        with store.connect(path) as connection:
            store.put_membrane(connection, "no_such_run",
                               [{"start_ns": 0.0, "apl_a2": 63.0}])


def test_membrane_rows_are_replaced_rather_than_appended(tmp_path: Path):
    path = tmp_path / "t.db"
    # Through the real path, which creates the run before its membrane rows.
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    chunks = [{"start_ns": 0.0, "apl_a2": 63.0, "thickness_a": 38.0}]
    for _ in range(3):
        with store.connect(path) as connection:
            store.put_membrane(connection, "r", chunks)
    with sqlite3.connect(path) as connection:
        assert connection.execute(
            "SELECT COUNT(*) FROM membrane").fetchone()[0] == 1


def test_summarise_reports_the_final_tenth(tmp_path: Path):
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    summary = store.summarise(path, "r")
    assert summary["span_ns"] == [0.0, 1.0]
    xdfg = summary["motifs"]["s:xDFG"]
    assert xdfg["start_s"] == pytest.approx(0.6)
    assert xdfg["final_s"] == pytest.approx(0.8)


# ===========================================================================
#  The pages
# ===========================================================================


# ===========================================================================
#  The reuse hook
# ===========================================================================


def test_family_detection_refuses_to_guess(tmp_path: Path):
    """A structure that looks like neither family must not default to one.

    Annotating a kinase as a receptor would resolve every motif to the wrong
    residues and still produce numbers.
    """
    from flexappeal.refstates.cli import detect_family

    # A short peptide with neither the HRD nor the DFG motif, written here
    # rather than committed.
    lines = []
    for index in range(1, 9):
        lines.append(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    "
            f"{index * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
    lines.append("END")
    path = tmp_path / "tiny.pdb"
    path.write_text("\n".join(lines) + "\n")
    assert detect_family(path) == "unknown"


def test_the_subcommand_is_reachable():
    """`FlexAppeal.py refstates annotate --help` has to work, not just import."""
    import subprocess
    import sys as _sys

    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [_sys.executable, str(root / "FlexAppeal.py"), "refstates", "annotate",
         "--help"], capture_output=True, text=True, cwd=root)
    assert result.returncode == 0, result.stderr[-400:]
    assert "--trajectory" in result.stdout
    assert "--family" in result.stdout


def test_a_bilayer_across_the_periodic_boundary_measures_correctly():
    """The case that made thickness read 78 A when it was 37.6 A.

    A trajectory stores wrapped coordinates, and under a membrane barostat the
    bilayer drifts until it straddles the box edge. The phosphorus z values are
    then bimodal at both ends with nothing in the middle, so an arithmetic
    centre lands in the empty gap and the two leaflets look a box-width apart.
    """
    box_z = 115.6
    # Two leaflets 37.6 A apart, placed so they wrap: one at z ~ 15, the other
    # at z ~ 95, which is the real distribution the b2ar run produced.
    lower = [15.0 + 0.5 * i for i in range(-4, 5)]
    upper = [95.0 + 0.5 * i for i in range(-4, 5)]
    z = lower + upper

    centred, shift = membrane.centre_bilayer(z, box_z)
    halves = membrane.split_leaflets(centred)
    assert len(halves[0]) == len(halves[1]) == 9

    import numpy as np

    gap = abs(float(np.mean(centred[halves[0]])) - float(np.mean(centred[halves[1]])))
    separation = min(gap, box_z - gap)
    # 15 to 95 the short way round the box is 115.6 - 80 = 35.6.
    assert separation == pytest.approx(35.6, abs=0.5)


def test_centring_is_a_no_op_for_a_bilayer_already_in_the_middle():
    """The fix must not move a bilayer that does not need moving."""
    import numpy as np

    box_z = 120.0
    z = [40.0, 41.0, 42.0, 78.0, 79.0, 80.0]
    centred, shift = membrane.centre_bilayer(z, box_z)
    assert abs(shift) < 1.0
    assert np.allclose(sorted(centred), sorted(z), atol=1.0)


# ===========================================================================
#  Section 10's named QC: metrics on geometries with known answers
# ===========================================================================


def _chain_from_positions(positions: dict) -> object:
    """A gemmi chain built in the test, with named atoms at chosen positions.

    `positions` maps a residue number to (residue name, {atom name: (x, y, z)}).
    Built rather than loaded so the expected answer is arithmetic.
    """
    import gemmi

    chain = gemmi.Chain("A")
    for number, (name, atoms) in sorted(positions.items()):
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(number, " ")
        for atom_name, (x, y, z) in atoms.items():
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element(atom_name[0])
            atom.pos = gemmi.Position(x, y, z)
            residue.add_atom(atom)
        chain.add_residue(residue)
    return chain


def test_a_salt_bridge_at_a_known_distance_measures_as_that_distance():
    """Section 10's hand-built salt bridge. 3.5 A in, 3.5 A out.

    Measured between the closest heavy atoms of the two side chains, not
    between one named atom of each: both groups are symmetric, and naming one
    atom makes the distance depend on which way the side chain was built.
    """
    from flexappeal.refstates.gpcr import GenericNumbering, min_heavy_distance

    numbering = GenericNumbering(entry_name="toy")
    numbering.by_generic = {"3.50": 131, "6.30": 268}
    numbering.residue_letter = {131: "R", 268: "E"}

    chain = _chain_from_positions({
        131: ("ARG", {"NE": (0.0, 0.0, 0.0), "CZ": (1.0, 0.0, 0.0),
                      "NH1": (2.0, 0.0, 0.0), "NH2": (1.0, 1.0, 0.0)}),
        # Nearest pair is NH1 at x = 2.0 against OE1 at x = 5.5, so 3.5 A.
        268: ("GLU", {"CD": (7.0, 0.0, 0.0), "OE1": (5.5, 0.0, 0.0),
                      "OE2": (8.0, 0.0, 0.0)}),
    })
    distance = min_heavy_distance(
        chain, numbering, "3.50", "6.30",
        first_atoms=("NE", "CZ", "NH1", "NH2"),
        second_atoms=("CD", "OE1", "OE2"))
    assert distance == pytest.approx(3.5, abs=1e-6)


def test_the_ionic_lock_threshold_is_applied_at_the_right_side():
    """4.0 A is the cut: just inside is intact, just outside is broken."""
    from flexappeal.refstates.gpcr import (IONIC_LOCK_A, GenericNumbering,
                                           ionic_lock)

    numbering = GenericNumbering(entry_name="toy")
    numbering.by_generic = {"3.50": 131, "6.30": 268}
    numbering.residue_letter = {131: "R", 268: "E"}

    for gap, expected in ((IONIC_LOCK_A - 0.2, "intact"),
                          (IONIC_LOCK_A + 0.2, "broken")):
        chain = _chain_from_positions({
            131: ("ARG", {"NH1": (0.0, 0.0, 0.0)}),
            268: ("GLU", {"OE1": (gap, 0.0, 0.0)}),
        })
        assert ionic_lock(chain, numbering)["state"] == expected, gap


def test_a_helix_rotated_by_a_known_angle_reports_that_angle():
    """Section 10's rotated helix. The tilt is read back from the rotation.

    `axis_tilt_deg` measures one rod's principal axis against z, so a helix
    built along z and tilted by a known angle in the xz plane must report it.
    """
    import numpy as np

    straight = np.array([[0.0, 0.0, 1.5 * i] for i in range(20)])
    for degrees in (0.0, 15.0, 30.0, 45.0):
        angle = np.radians(degrees)
        rotation = np.array([[np.cos(angle), 0.0, np.sin(angle)],
                             [0.0, 1.0, 0.0],
                             [-np.sin(angle), 0.0, np.cos(angle)]])
        tilted = straight @ rotation.T
        assert membrane.axis_tilt_deg(tilted) == pytest.approx(degrees, abs=1e-4)


def _seven_helix_bundle(tilt_deg: float, radius_a: float = 18.0,
                        per_helix: int = 28, rise_a: float = 1.35):
    """A synthetic 7TM bundle, tilted by a known angle in the xz plane.

    Seven rods on a circle, each running the opposite way to its neighbour the
    way consecutive TM helices do, then the whole assembly rotated about y. The
    radius and the rise are chosen so the Calpha cloud comes out about as wide
    as it is tall, which is the case the bundle axis has to survive.
    """
    import numpy as np

    coords, residues = [], []
    residue = 0
    for helix in range(7):
        angle = 2.0 * np.pi * helix / 7.0
        x, y = radius_a * np.cos(angle), radius_a * np.sin(angle)
        direction = 1.0 if helix % 2 == 0 else -1.0
        for step in range(per_helix):
            z = direction * rise_a * (step - 0.5 * (per_helix - 1))
            coords.append([x, y, z])
            residues.append(residue)
            residue += 1
        # A gap in the numbering between helices, which is what the loops are.
        residue += 10
    rotation_angle = np.radians(tilt_deg)
    rotation = np.array([[np.cos(rotation_angle), 0.0, np.sin(rotation_angle)],
                         [0.0, 1.0, 0.0],
                         [-np.sin(rotation_angle), 0.0, np.cos(rotation_angle)]])
    return np.array(coords) @ rotation.T, residues


def test_the_bundle_tilt_is_the_mean_of_the_helix_axes_not_the_cloud_axis():
    """The regression behind FINDINGS 4n: a round cloud has no principal axis.

    The old measurement took the principal axis of all the membrane-spanning
    Calpha at once. For a 7TM bundle that cloud is nearly isotropic, so the
    direction it returns is noise. This builds a bundle with a tilt that is
    known by construction and asserts the two measurements apart: the per-helix
    mean recovers it, the whole-cloud axis does not.
    """
    import numpy as np

    for tilt in (0.0, 5.0, 12.0):
        coords, residues = _seven_helix_bundle(tilt)
        assert membrane.bundle_tilt(coords, residues) == pytest.approx(
            tilt, abs=0.5), f"bundle tilt wrong at {tilt} degrees"

    # And the cloud really is round, so the old route really was noise.
    coords, residues = _seven_helix_bundle(5.0)
    centred = coords - coords.mean(axis=0)
    values = np.linalg.svd(centred, full_matrices=False)[1]
    assert values[0] / values[1] < 1.5, "the fixture is not a round cloud"
    assert abs(membrane.axis_tilt_deg(coords) - 5.0) > 20.0, (
        "the whole-cloud axis happened to be right, so this fixture no longer "
        "demonstrates the fault it was built for")


def test_helix_axes_running_opposite_ways_do_not_cancel():
    """Consecutive TM helices alternate direction, so the signs must align.

    Without folding each axis towards +z first, a mean over seven alternating
    rods tends towards zero and the tilt it reports is whatever is left, which
    for an untilted bundle is 90 degrees rather than nothing.
    """
    coords, residues = _seven_helix_bundle(0.0)
    assert membrane.bundle_tilt(coords, residues) == pytest.approx(0.0, abs=0.5)


def test_a_bundle_needs_more_than_two_helices_to_have_an_axis():
    """Two rods are not a bundle, and saying nothing beats saying 48 degrees."""
    coords, residues = _seven_helix_bundle(0.0)
    two = [i for i, r in enumerate(residues) if r < 48]
    assert membrane.bundle_tilt(coords[two], [residues[i] for i in two]) is None


def test_an_axis_through_a_round_run_is_refused_rather_than_averaged_in():
    """A run that is not rod shaped contributes no axis.

    `helix_axis` is averaged into the bundle, so a run whose direction is
    undetermined would drag the mean towards a direction nothing measured.
    """
    import numpy as np

    rng = np.random.default_rng(0)
    assert membrane.helix_axis(rng.normal(size=(30, 3))) is None
    rod = np.array([[0.0, 0.0, 1.5 * i] for i in range(20)])
    assert membrane.helix_axis(rod) is not None


def test_the_bundle_tilt_rejects_mismatched_residue_indices():
    coords, residues = _seven_helix_bundle(0.0)
    with pytest.raises(ValueError, match="run in parallel"):
        membrane.bundle_tilt(coords, residues[:-1])


def test_the_slab_keeps_the_calpha_between_the_phosphate_planes():
    """`slab_ca_indices` is what separates the bundle from the loops."""
    protein_z = [0.0, 10.0, -10.0, 25.0, -25.0, 18.4, 18.6]
    keep = membrane.slab_ca_indices(None, protein_z, midplane=0.0,
                                    half_thickness=18.5)
    assert keep == [0, 1, 2, 5]


def test_the_cli_runs_as_a_script_and_not_only_as_an_import():
    """Section 11's deliverable is a CLI, so the script path has to work.

    Everything else here reaches `cli` by importing it, and an import defines
    every function before anything calls one. Running the module as a script
    does not: the `__main__` guard sat directly under `main`, above the
    subcommand functions, so `build_parser` reached `set_defaults(func=controls)`
    before `controls` existed and every invocation died on a `NameError`. The
    one path a user would take was the only one nothing covered.

    A subprocess rather than calling `main`, because calling `main` is the
    import path again and would have passed throughout.
    """
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    for argv in (["--help"], ["build", "--help"], ["controls", "--help"],
                 ["report", "--help"], ["annotate", "--help"], ["run", "--help"]):
        done = subprocess.run(
            [sys.executable, "-m", "flexappeal.refstates.cli", *argv],
            cwd=root, capture_output=True, text=True)
        assert done.returncode == 0, (
            f"`cli {' '.join(argv)}` exited {done.returncode}: {done.stderr}")
        assert "Traceback" not in done.stderr, done.stderr

    # And every subcommand the deliverable names is actually wired up.
    listed = subprocess.run(
        [sys.executable, "-m", "flexappeal.refstates.cli", "--help"],
        cwd=root, capture_output=True, text=True).stdout
    for name in ("build", "run", "annotate", "report"):
        assert name in listed, f"section 11 names `{name}` and it is not there"


def test_the_state_score_is_symmetric_between_two_references():
    """Section 10's two-reference toy: s(active ref) = -s(inactive ref).

    A frame sitting exactly on the active reference must score +1 and one on
    the inactive reference -1. An asymmetry here would mean the score prefers
    one state for a reason that is not in the data.
    """
    separation = 7.3
    on_active = state.score_motif("m", active={"A": 0.0},
                                  inactive={"I": separation}, n_residues=5)
    on_inactive = state.score_motif("m", active={"A": separation},
                                    inactive={"I": 0.0}, n_residues=5)
    assert on_active.s == pytest.approx(1.0)
    assert on_inactive.s == pytest.approx(-1.0)
    assert on_active.s == pytest.approx(-on_inactive.s)

    # And halfway between, for any separation, is exactly zero.
    midway = state.score_motif("m", active={"A": separation / 2},
                               inactive={"I": separation / 2}, n_residues=5)
    assert midway.s == pytest.approx(0.0)
    assert midway.geometric_label == "intermediate"


def test_a_run_resumes_from_its_checkpoint(tmp_path: Path):
    """Section 10's resume test, on a system small enough to be a unit test.

    Resumability is load-bearing rather than a convenience: one replica per row
    at 50 ns means a run that cannot resume costs its whole row, and the
    membrane rows take 20 hours each on a machine that sleeps. This asserts the
    property the long runs depend on, that a reloaded checkpoint continues the
    clock rather than restarting it.
    """
    import openmm as mm
    from openmm import app, unit

    # Three argon-like particles with no forces: the integrator is the thing
    # under test, not the force field.
    system = mm.System()
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("AR", chain)
    for _ in range(3):
        system.addParticle(39.948 * unit.amu)
        topology.addAtom("Ar", app.element.argon, residue)
    positions = [mm.Vec3(0, 0, 0), mm.Vec3(0.5, 0, 0),
                 mm.Vec3(0, 0.5, 0)] * unit.nanometer

    def simulation():
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
        return app.Simulation(topology, system, integrator,
                              mm.Platform.getPlatformByName("Reference"))

    first = simulation()
    first.context.setPositions(positions)
    first.context.setVelocitiesToTemperature(300 * unit.kelvin)
    first.step(50)
    checkpoint = tmp_path / "state.chk"
    first.saveCheckpoint(str(checkpoint))
    before = first.context.getState().getTime().value_in_unit(unit.picoseconds)
    assert before > 0

    second = simulation()
    second.loadCheckpoint(str(checkpoint))
    after = second.context.getState().getTime().value_in_unit(unit.picoseconds)
    assert after == pytest.approx(before), "the clock restarted instead of resuming"

    second.step(50)
    final = second.context.getState().getTime().value_in_unit(unit.picoseconds)
    assert final > after, "stepping after a resume did not advance the clock"


def test_the_runner_treats_a_checkpoint_as_the_resume_point(tmp_path: Path):
    """`run_md.py` must skip equilibration when a checkpoint exists.

    Read from the source rather than by running it, because running it needs a
    solvated system. What matters is that the resume branch comes before the
    equilibration branch, so a resumed run does not heat an already-warm system.
    """
    source = (Path(__file__).resolve().parents[1] / "tools" / "run_md.py").read_text()
    resume = source.index("resuming = checkpoint.is_file()")
    heating = source.index("heating to")
    assert resume < heating, "equilibration is reached before the resume check"
    assert "loadCheckpoint" in source
    assert "already at" in source, "a completed run must report and stop"


def test_a_resumed_run_adds_the_barostat_before_the_context_is_built():
    """The branch that adds the barostat is the branch a resume skips.

    This is read from the source because the fault was an ordering one and the
    cost of getting it wrong is invisible: OpenMM restores `MonteCarloPressure`
    from the checkpoint whether or not a force consumes it, so a resumed run
    reported the barostat's parameters while having none, and produced at
    constant volume. The force has to be on the System before the Context
    exists, so what matters is that the `addForce` precedes the `Simulation`.
    """
    source = (Path(__file__).resolve().parents[1] / "tools" / "run_md.py").read_text()
    added = source.index("system.addForce(barostat(bool(args.membrane)")
    built = source.index("simulation = app.Simulation(")
    loaded = source.index("simulation.loadCheckpoint(")
    assert added < built, "the barostat is added after the context is built"
    assert built < loaded, "the checkpoint is loaded before the context exists"


def _argon_box(with_barostat: bool):
    """Eight particles in a periodic box, with or without a barostat.

    Small enough to run on the Reference platform in milliseconds, which is
    what makes the behaviour below testable at all: the real fault was found on
    a 194,032-atom membrane row.
    """
    import openmm as mm
    from openmm import unit

    system = mm.System()
    for _ in range(8):
        system.addParticle(39.948 * unit.amu)
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(2, 0, 0) * unit.nanometer,
        mm.Vec3(0, 2, 0) * unit.nanometer,
        mm.Vec3(0, 0, 2) * unit.nanometer)
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(0.5 * unit.nanometer)
    for _ in range(8):
        nonbonded.addParticle(0.0, 0.34 * unit.nanometer,
                              0.996 * unit.kilojoule_per_mole)
    system.addForce(nonbonded)
    if with_barostat:
        system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin))
    return system


def _context(system):
    import openmm as mm
    from openmm import unit

    return mm.Context(
        system,
        mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    0.002 * unit.picoseconds),
        mm.Platform.getPlatformByName("Reference"))


def test_a_checkpoint_reports_a_barostat_it_does_not_have():
    """Why the ordering above matters, demonstrated rather than asserted.

    The test above reads the source, because the real path needs a built
    membrane system and a GPU. What it cannot show is why getting the order
    wrong is invisible, and that is the whole reason the fault survived: a
    checkpoint written with a barostat loads without complaint into a system
    that has none, and the context then reports `MonteCarloPressure` as though
    it did. A resumed run printed the barostat's own pressure while
    integrating at constant volume.

    Eight particles are enough to show it.
    """
    import openmm as mm

    source = _argon_box(with_barostat=True)
    written = _context(source)
    written.setPositions([mm.Vec3(0.3 * (i % 2), 0.3 * ((i // 2) % 2),
                                  0.3 * (i // 4)) for i in range(8)])
    checkpoint = written.createCheckpoint()

    # The resume done wrongly: the force is never added back.
    without = _argon_box(with_barostat=False)
    assert not any(isinstance(f, mm.MonteCarloBarostat)
                   for f in without.getForces()), "the fixture is wrong"
    broken = _context(without)
    broken.loadCheckpoint(checkpoint)
    assert "MonteCarloPressure" in dict(broken.getParameters()), (
        "if OpenMM ever stops restoring this parameter the fault becomes "
        "visible on its own and this test should be revisited")

    # The resume done the way `run_md.py` does it: the force first, then the
    # context, then the checkpoint.
    with_force = _argon_box(with_barostat=True)
    fixed = _context(with_force)
    fixed.loadCheckpoint(checkpoint)
    assert any(isinstance(f, mm.MonteCarloBarostat)
               for f in with_force.getForces())
    assert "MonteCarloPressure" in dict(fixed.getParameters())


# ===========================================================================
#  Section 10's page QC, as far as it goes without a browser
# ===========================================================================
#
# Section 10 asks for Playwright at two viewports with no console errors, and
# axe-core with no serious or critical violations. Neither can run here:
# playwright is not in the environment and there is no preinstalled Chromium,
# and the spec says never to run `playwright install`. Recorded as a gap in
# FINDINGS rather than skipped quietly.
#
# What is checkable server-side is the part that causes most console errors in
# practice, a referenced asset that is not there, and the subset of axe-core's
# rules that are properties of the markup rather than of the rendered page.

def test_the_site_has_an_icon_at_every_place_a_browser_looks():
    """The site had none, so every tab showed the default document glyph.

    Both paths are asserted: a browser requests /favicon.ico whatever the link
    tags say, and a 404 per page view is noise in the access log that reads as
    a fault.
    """
    from flexappeal.webapp import create_app

    client = create_app().test_client()
    assert client.get("/favicon.ico").status_code == 200
    for name in ("favicon.ico", "favicon-16x16.png", "favicon-32x32.png",
                 "apple-touch-icon.png", "android-chrome-192x192.png",
                 "android-chrome-512x512.png", "site.webmanifest"):
        response = client.get(f"/static/icons/{name}")
        assert response.status_code == 200, name
        assert len(response.data) > 200, f"{name} is suspiciously small"

    head = client.get("/").get_data(as_text=True)
    for rel in ('rel="icon"', 'rel="apple-touch-icon"', 'rel="manifest"'):
        assert rel in head, rel


def test_the_drift_fit_recovers_a_known_slope_and_knows_when_it_cannot(tmp_path: Path):
    """Section 5.2's drift check, on series whose answers are arithmetic.

    Two properties matter and they pull in opposite directions: a real slope has
    to be measured, and a slope smaller than its own noise has to be reported as
    unresolved rather than quoted. A tool that only did the first would pass a
    run on a figure the data cannot support.
    """
    import importlib.util

    import numpy as np

    script = Path(__file__).resolve().parents[1] / "tools" / "check_drift.py"
    spec = importlib.util.spec_from_file_location("_check_drift", script)
    drift = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(drift)

    # A clean ramp: slope 250 per unit x, recovered exactly.
    x = [float(i) for i in range(40)]
    y = [1000.0 + 250.0 * v for v in x]
    slope, error = drift._fit(x, y)
    assert slope == pytest.approx(250.0, rel=1e-9)
    assert error < 1e-6
    assert abs(slope) > error, "a clean ramp must be resolved"

    # Pure noise about a flat line, across many draws rather than one. A single
    # draw is the wrong test: with one-sigma the first seed tried gave a 1.7
    # sigma fluctuation and looked like a real drift, which is what sent the
    # tool from a one-sigma criterion to three. The property that matters is
    # that false positives are rare, so it is measured as a rate.
    false_positives = 0
    trials = 60
    for seed in range(trials):
        rng = np.random.RandomState(seed)
        flat = [1000.0 + float(rng.normal(0, 500)) for _ in x]
        slope, error = drift._fit(x, flat)
        if abs(slope) > drift.RESOLVED_SIGMA * error:
            false_positives += 1
    assert false_positives <= 1, (
        f"{false_positives} of {trials} noise draws read as a real drift at "
        f"{drift.RESOLVED_SIGMA} sigma")

    # And a drift well above the noise still has to be caught.
    rng = np.random.RandomState(7)
    ramped = [1000.0 + 300.0 * v + float(rng.normal(0, 500)) for v in x]
    slope, error = drift._fit(x, ramped)
    assert abs(slope) > drift.RESOLVED_SIGMA * error
    assert slope == pytest.approx(300.0, rel=0.1)


def test_degrees_of_freedom_subtracts_constraints(tmp_path: Path):
    """3N would overstate the dof threefold for a constrained system."""
    import importlib.util

    import openmm as mm
    from openmm import unit

    script = Path(__file__).resolve().parents[1] / "tools" / "check_drift.py"
    spec = importlib.util.spec_from_file_location("_check_drift", script)
    drift = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(drift)

    system = mm.System()
    for _ in range(10):
        system.addParticle(1.0 * unit.amu)
    system.addConstraint(0, 1, 0.1 * unit.nanometer)
    system.addConstraint(2, 3, 0.1 * unit.nanometer)
    system.addForce(mm.CMMotionRemover())
    path = tmp_path / "system.xml"
    path.write_text(mm.XmlSerializer.serialize(system))

    # 3 * 10 moving - 2 constraints - 3 for the COM remover.
    assert drift.degrees_of_freedom(path) == 25
    assert drift.degrees_of_freedom(tmp_path / "absent.xml") is None


def test_rmsd_to_nearest_never_compares_a_reference_with_itself():
    """Otherwise one state's column is all zeros and the map is two lines."""
    import numpy as np

    from flexappeal.refstates import ensembles as ens

    built = {}
    for state, offsets in (("active", (0.0, 0.5)), ("inactive", (6.0, 6.4))):
        holder = ens.MotifEnsemble(state=state)
        for index, offset in enumerate(offsets):
            holder.by_code[f"{state[:3].upper()}{index}"] = {
                ens.CORE_KEY: np.zeros((4, 3)),
                "NPxxY": np.array([[offset, 0.0, 0.0]] * 4),
            }
        built[state] = holder

    found = ens.rmsd_to_nearest(built, "NPxxY", "inactive")
    # An inactive structure's distance to the inactive set is its distance to
    # the other inactive structure, 0.4, not zero.
    assert found["INA0"] == pytest.approx(0.4, abs=1e-6)
    assert found["INA1"] == pytest.approx(0.4, abs=1e-6)
    # And an active one is genuinely far.
    assert found["ACT0"] == pytest.approx(6.0, abs=1e-6)


# ===========================================================================
#  The invariant behind three separate bugs
# ===========================================================================
#
# Three defects in this build were the same shape: an endpoint returned 200 with
# data nothing could draw. An empty panel whose fetch 404'd, a state map with
# every y null, and a motif series with eight metric rows and no points. None
# shows up in a status code and none raises.
#
# So the invariant is asserted directly: a 200 from any of these endpoints must
# carry arrays that are non-empty and the same length as each other. A panel is
# allowed to be absent; it is not allowed to be present and blank.


def _lopsided_payload() -> dict:
    """A run where one motif resolved only its inactive distance.

    Legal, because the rule is that an unmeasured key is omitted rather than
    written as zero, so every reader has to cope with it.
    """
    return {
        "name": "lopsided", "kind": "kinase", "core": {},
        "motif_residues": {"one_sided": [1, 2, 3], "complete": [4, 5]},
        "ensembles": {"active": {"codes": ["2F4J"]},
                      "inactive": {"codes": ["3PYY"]}},
        "frames": [{"time_ns": t / 10.0, "n_conflicts": 0, "motifs": [
            {"motif": "one_sided", "d_inact_a": 5.0 + t * 0.1},
            {"motif": "complete", "d_act_a": 1.0, "d_inact_a": 4.0,
             "s": 0.6, "label": "active"}]} for t in range(8)],
        "notes": [],
    }


# ===========================================================================
#  The pair, which is the result the build is actually for
# ===========================================================================


def _pair_db(path: Path, active_s: float, inactive_s: float) -> None:
    """Two runs of one system from opposite start states."""
    def payload(value: float) -> dict:
        return {
            "name": "pair", "kind": "gpcr", "core": {},
            "motif_residues": {"TM6": list(range(269, 300)),
                               "NPxxY": [322, 323, 324, 325, 326]},
            "ensembles": {"active": {"codes": ["3P0G"]},
                          "inactive": {"codes": ["2RH1"]}},
            "frames": [{"time_ns": t / 10.0, "n_conflicts": 0, "motifs": [
                {"motif": "TM6", "d_act_a": 1.0, "d_inact_a": 4.0,
                 "s": value, "label": "x"},
                {"motif": "NPxxY", "d_act_a": 1.0, "d_inact_a": 3.0,
                 "s": value, "label": "x"}]} for t in range(10)],
            "notes": [],
        }

    store.ingest_analysis(path, payload(active_s), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC",
                          run_fields={"start_state": "active",
                                      "status": "complete"})
    store.ingest_analysis(path, payload(inactive_s), system_id="b2ar",
                          run_id="b2ar_carazolol", family="gpcr",
                          environment="POPC",
                          run_fields={"start_state": "inactive",
                                      "status": "complete"})


def test_the_whole_score_is_weighted_by_the_motifs_own_residue_counts(tmp_path: Path):
    """Recomputed from the stored per-motif scores, not stored separately.

    One definition of the weighting, so it cannot drift between what the
    analysis wrote and what a page shows. TM6 is 31 residues against NPxxY's
    five, so a disagreement between them must land nearer TM6.
    """
    path = tmp_path / "t.db"

    def payload(tm6: float, npxxy: float) -> dict:
        return {
            "name": "w", "kind": "gpcr", "core": {},
            "motif_residues": {"TM6": list(range(269, 300)),
                               "NPxxY": [322, 323, 324, 325, 326]},
            "ensembles": {"active": {"codes": ["3P0G"]}},
            "frames": [{"time_ns": 0.0, "n_conflicts": 0, "motifs": [
                {"motif": "TM6", "d_act_a": 1.0, "d_inact_a": 4.0, "s": tm6},
                {"motif": "NPxxY", "d_act_a": 1.0, "d_inact_a": 3.0,
                 "s": npxxy}]}],
            "notes": [],
        }

    store.ingest_analysis(path, payload(1.0, -1.0), system_id="b2ar",
                          run_id="r", family="gpcr", environment="POPC")
    score = store.whole_score(path, "r")
    # 31 at +1 against 5 at -1 is +26/36.
    assert score["final_s"] == pytest.approx(26 / 36, abs=1e-6)
    assert score["n_motifs"] == 2


def test_every_name_the_tools_import_from_the_package_exists():
    """A deletion in a module the unit tests touch lightly is invisible.

    This caught nothing when it was written and would have caught a real
    deletion: a rewrite of `traj.load` replaced everything from that function to
    the end of the file and removed `deposited_numbering`, `renumber`,
    `scaffold_core` and `core_agreement` with it. Every unit test still passed,
    because none of them imports those four, and the loss only surfaced on
    running the real tool.

    So the import surface the scripts actually use is asserted directly, by
    reading what they import rather than by listing it here, which is what makes
    it stay true as the scripts change.
    """
    import importlib
    import re

    root = Path(__file__).resolve().parents[1]
    scripts = sorted((root / "tools").glob("*.py")) + [
        root / "flexappeal" / "refstates" / "cli.py",
        root / "flexappeal" / "refstates" / "panels.py",
    ]

    missing: list[str] = []
    checked = 0
    # Both forms the scripts actually use: a submodule (`...refstates.traj
    # import x`) and a package attribute (`...refstates import traj, pocket`),
    # with leading whitespace because most are inside functions.
    # `\.(?!\.)` is a single leading dot, so a parent-relative
    # `from ..webapp import x` does not get read as
    # `flexappeal.refstates.webapp`, which it is not. A lookbehind does not do
    # it: the character before the first dot is a space either way.
    pattern = re.compile(
        r"from\s+(?:flexappeal\.refstates|\.(?!\.))"
        r"(?:\.?([a-z_]+))?\s+import\s+([^\n(]+)")
    for script in scripts:
        for module, names in pattern.findall(script.read_text()):
            target = f"flexappeal.refstates.{module}" if module else "flexappeal.refstates"
            try:
                loaded = importlib.import_module(target)
            except ModuleNotFoundError:
                missing.append(f"{script.name}: module {target}")
                continue
            # A trailing comment is part of the line and not part of a name.
            # `from ... import DECISIVE_S  # noqa: E402` is ordinary Python,
            # and reading the comment as the name reported a real import as
            # absent, which is a fault in this check rather than in the script
            # it was checking.
            names = names.split("#", 1)[0]
            for name in (n.strip() for n in names.split(",")):
                if not name or name == "as":
                    continue
                name = name.split(" as ")[0].strip()
                checked += 1
                if hasattr(loaded, name):
                    continue
                # A submodule is only an attribute of its package once it has
                # been imported, so `from flexappeal.refstates import pocket`
                # is valid while hasattr is False until then.
                try:
                    importlib.import_module(f"{target}.{name}")
                except ModuleNotFoundError:
                    missing.append(f"{script.name}: {target}.{name}")

    assert checked > 20, f"only {checked} imports checked; the pattern missed some"
    assert not missing, f"imported but absent: {missing}"


# ===========================================================================
#  P5 and P6, exercised before the data lands
# ===========================================================================
#
# "Pending" is not evidence that a control works. Both are run against
# fixtures here so that when the rows land the only new thing is the data.


def _cli():
    import importlib

    return importlib.import_module("flexappeal.refstates.cli")


def _inactive_run(path: Path, scores: list[float]) -> None:
    """An inactive-start run whose whole-protein score follows `scores`."""
    payload = {
        "name": "p5", "kind": "gpcr", "core": {},
        "motif_residues": {"TM6": list(range(269, 300))},
        "ensembles": {"inactive": {"codes": ["2RH1"]}},
        "frames": [{"time_ns": float(i), "n_conflicts": 0, "motifs": [
            {"motif": "TM6", "d_act_a": 4.0, "d_inact_a": 1.0, "s": s,
             "label": "inactive"}]} for i, s in enumerate(scores)],
        "notes": [],
    }
    store.ingest_analysis(path, payload, system_id="b2ar",
                          run_id="b2ar_carazolol", family="gpcr",
                          environment="POPC",
                          run_fields={"start_state": "inactive",
                                      "status": "complete"})


def test_p5_passes_a_run_that_stays_inactive(tmp_path: Path):
    path = tmp_path / "t.db"
    # 50 frames at 1 ns, all negative: the final 20 ns are all below zero.
    _inactive_run(path, [-0.8] * 50)
    rows = _cli().control_p5(path)
    assert len(rows) == 1
    assert rows[0]["passed"] is True
    assert rows[0]["fraction_below_zero"] == pytest.approx(1.0)
    assert rows[0]["replicas"] == 1


def test_p5_fails_a_run_that_activates_late(tmp_path: Path):
    """A run that crosses inside the window must fail, not be averaged over.

    The score spends most of the run inactive and ends up positive, so a mean
    over the whole trajectory would pass it. P5 asks about the final 20 ns.
    """
    path = tmp_path / "t.db"
    _inactive_run(path, [-0.8] * 35 + [0.4] * 15)
    rows = _cli().control_p5(path)
    assert rows[0]["passed"] is False
    assert rows[0]["fraction_below_zero"] < 1.0
    # And the window really is the last 20 ns rather than the whole run.
    assert rows[0]["frames_in_window"] <= 21


def test_p5_ignores_active_start_runs(tmp_path: Path):
    """P5 is about the inactive-start runs; an active one is not a failure."""
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s",
                          run_id="active_one", family="kinase",
                          environment="water",
                          run_fields={"start_state": "active"})
    assert _cli().control_p5(path) == []


def test_p6_passes_a_bilayer_inside_the_acceptable_window(tmp_path: Path):
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    with store.connect(path) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": float(i * 10), "apl_a2": 64.0 + i * 0.3,
             "thickness_a": 38.0} for i in range(5)])
    rows = _cli().control_p6(path)
    assert rows[0]["passed"] is True
    assert rows[0]["apl_in_expected"] is True
    assert rows[0]["thickness_in_expected"] is True


def test_p6_fails_on_g6_but_tolerates_outside_expected(tmp_path: Path):
    """Section 5.8 has two windows and P6 is tied to the wider one.

    A chunk outside 60 to 68 is logged; a run outside 55 to 75 invalidates the
    bilayer. Conflating them would fail a usable run.
    """
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    with store.connect(path) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": 0.0, "apl_a2": 71.0, "thickness_a": 38.0},
            {"start_ns": 10.0, "apl_a2": 70.5, "thickness_a": 38.0}])
    rows = _cli().control_p6(path)
    assert rows[0]["apl_in_expected"] is False, "71 is outside 60 to 68"
    assert rows[0]["passed"] is True, "but inside 55 to 75, so usable"

    with store.connect(path) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": 0.0, "apl_a2": 40.0, "thickness_a": 38.0}])
    assert _cli().control_p6(path)[0]["passed"] is False


# --- and the line a reader of P6 actually sees --------------------------------
#
# 4bn was a number rounded onto its own threshold in a summary document. The
# control that produced the number printed it the same way: one decimal, so
# `b2ar_apo`'s widest chunk read "68.0" against an expected window ending at
# 68.0, beside the word PASS. The verdict was right, because P6 is tied to the
# wider window, and the line gave a reader no way to see that one chunk had
# left the narrower one.

def test_p6_prints_a_chunk_on_the_boundary_without_rounding_it_in(
        tmp_path: Path, capsys):
    """The printed evidence, not only the dict behind it."""
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    with store.connect(path) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": 0.0, "apl_a2": 68.05, "thickness_a": 38.42},
            {"start_ns": 10.0, "apl_a2": 66.23, "thickness_a": 39.35}])
    cli = _cli()
    cli.print_p6(cli.control_p6(path))
    said = capsys.readouterr().out

    assert "68.05" in said, "a chunk outside the window must not print as 68.0"
    assert "68.0," not in said and "68.0)" not in said, said
    assert "PASS" in said, "P6 is tied to the acceptable window and this passes"
    assert "1 of 2 outside the expected 60 to 68" in said, said
    assert "inside the acceptable 55 to 75" in said, said
    # Units on a readout, which is section 10's rule and not only about plots.
    assert "\u00c5\u00b2" in said, said


def test_p6_says_nothing_about_windows_a_bilayer_never_left(
        tmp_path: Path, capsys):
    """The note is a finding, so a clean run does not carry one."""
    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_carazolol", family="gpcr",
                          environment="POPC")
    with store.connect(path) as connection:
        store.put_membrane(connection, "b2ar_carazolol", [
            {"start_ns": 0.0, "apl_a2": 65.71, "thickness_a": 38.49},
            {"start_ns": 10.0, "apl_a2": 67.84, "thickness_a": 39.38}])
    cli = _cli()
    cli.print_p6(cli.control_p6(path))
    said = capsys.readouterr().out

    assert "PASS" in said
    assert "outside" not in said, said


def test_the_figures_are_drawn_from_whatever_the_database_holds(tmp_path: Path):
    """Section 11 asks for FINDINGS "with figures", so they have to generate.

    Run against a database built here rather than against the real one, so the
    test says whether the generator works and not whether production has got
    far enough. A membrane figure and a motif figure are both exercised,
    because they take different routes through the store.
    """
    import subprocess
    import sys

    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    with store.connect(database) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": float(i * 10), "apl_a2": 64.0 + i * 0.3,
             "thickness_a": 38.0, "tilt_deg": 7.0 + i} for i in range(4)])

    out = tmp_path / "figures"
    root = Path(__file__).resolve().parents[1]
    done = subprocess.run(
        [sys.executable, str(root / "tools" / "make_figures.py"),
         "--db", str(database), "--out", str(out)],
        cwd=root, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr

    drawn = sorted(f.name for f in out.glob("*.png"))
    assert "b2ar_apo_motifs.png" in drawn, drawn
    assert "b2ar_apo_bilayer.png" in drawn, drawn
    # Every file is a real PNG and not a zero-length one.
    for figure in out.glob("*.png"):
        assert figure.stat().st_size > 1000, f"{figure.name} is empty"
        assert figure.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_a_figure_with_no_data_is_skipped_rather_than_drawn_empty(tmp_path: Path):
    """An axis with nothing on it looks like a result of zero, not an absence."""
    import subprocess
    import sys

    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    out = tmp_path / "figures"
    root = Path(__file__).resolve().parents[1]
    done = subprocess.run(
        [sys.executable, str(root / "tools" / "make_figures.py"),
         "--db", str(database), "--out", str(out)],
        cwd=root, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert not (out / "b2ar_apo_bilayer.png").exists(), (
        "a run with no membrane chunks drew a bilayer figure anyway")
    assert "no membrane chunks" in done.stdout, done.stdout


def test_the_deploy_leaves_the_working_data_behind_and_takes_the_database():
    """What the droplet needs, and nothing that only this machine needs.

    `data/md/` is the raw simulation output: 565 MB with one row a sixth of the
    way through, several gigabytes once four rows have finished. The web app
    reads none of it. The database it does read is gitignored, so this cannot
    be done with an rsync gitignore filter without shipping a droplet whose
    reference-states tab is silently empty, and an explicit list is a thing
    that can rot. Hence a test.
    """
    deploy = (Path(__file__).resolve().parents[1] / "deploy" / "deploy.sh").read_text()
    rsync = deploy[deploy.index("rsync "):]
    rsync = rsync[:rsync.index("DROPLET_SSH")]

    for path in ("data/md/", "data/cache/", "data/refs/"):
        assert f"--exclude '{path}'" in rsync, (
            f"{path} is working data the app never reads and would be synced")
    for needed in ("data/refstates.db", "data/analysis", "data/systems"):
        assert f"--exclude '{needed}" not in rsync, (
            f"{needed} is read at runtime and must not be excluded")
    # The filter that would have been tidier and wrong, kept out deliberately.
    assert ".gitignore" not in rsync, (
        "a gitignore filter would exclude data/refstates.db, which is "
        "gitignored and is exactly what the explorer reads")


def test_the_queue_keeps_an_append_only_md_manifest(tmp_path: Path,
                                                    monkeypatch):
    """Section 0 asks for `data/manifests/<stage>.jsonl`; section 2.3 names
    this one. There was none.

    The only record of what the queue had done was `production.log`, which is
    prose, and the per-row `run.json`, which exists only once a row finishes.
    A row that was running, or that had died, left nothing a later pass could
    read.
    """
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "queue_md", root / "tools" / "queue_md.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)

    manifest = tmp_path / "md.jsonl"
    monkeypatch.setattr(queue, "MANIFEST", manifest)
    queue.note("started", "abl1_apo", production_ns=50.0)
    queue.note("failed", "abl1_apo", returncode=1)

    lines = [json.loads(line) for line in
             manifest.read_text().strip().splitlines()]
    assert [row["event"] for row in lines] == ["started", "failed"], (
        "the manifest is a history, so a row that started and then failed is "
        "two lines and not one overwritten one")
    assert all(row["name"] == "abl1_apo" for row in lines)
    assert all(row.get("at") for row in lines)


def test_a_queue_that_cannot_write_its_notes_still_runs(tmp_path: Path,
                                                        monkeypatch, capsys):
    """Losing the record is bad; losing a day of GPU time would be worse."""
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "queue_md", root / "tools" / "queue_md.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)

    # A path whose parent is a file, so mkdir and open both fail.
    blocked = tmp_path / "wall"
    blocked.write_text("not a directory")
    monkeypatch.setattr(queue, "MANIFEST", blocked / "md.jsonl")
    queue.note("started", "abl1_apo")           # must not raise
    assert "could not write" in capsys.readouterr().out


def test_a_dead_row_is_recorded_as_failed_and_not_left_running(tmp_path: Path):
    """Section 0's ladder ends in `status = "failed:<reason>"`. Nothing wrote it.

    `run_md.py` writes `run.json` on success only, and `tools/ingest.py` reads
    an absent `run.json` as a row still in flight. A row that died therefore
    showed as `running` on the explorer indefinitely, beside a partial
    trajectory that made it look plausible.
    """
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "queue_md", root / "tools" / "queue_md.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)

    run_dir = tmp_path / "abl1_apo"
    queue.record_failure(run_dir, "abl1_apo", 1, 3600.0,
                         ["[03:14:00] step 1200000",
                          "[03:15:00] Particle coordinate is NaN"])
    written = json.loads((run_dir / "run.json").read_text())
    assert written["status"].startswith("failed:"), written
    assert "NaN" in written["status"], "the reason should name what happened"
    assert written["wall_hours"] == 1.0

    # And a row that did finish keeps its own record, which is the better one.
    good = tmp_path / "b2ar_apo"
    good.mkdir()
    (good / "run.json").write_text(json.dumps({"status": "complete",
                                               "ns_per_day": 59.7}))
    queue.record_failure(good, "b2ar_apo", 1, 10.0, ["[03:00:00] Error"])
    assert json.loads((good / "run.json").read_text())["status"] == "complete"


def _tiny_bilayer(tmp_path: Path, n_frames: int = 4, lipids: bool = True):
    """A POPC slab and a couple of protein residues, written as pdb plus dcd.

    Built here rather than committed, and not a real structure: two planes of
    phosphorus in POPC residues is all `phosphorus_indices` looks for. Returns
    the topology and trajectory paths.
    """
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    if lipids:
        chain = topology.add_chain()
        for _ in range(8):
            residue = topology.add_residue("POPC", chain)
            topology.add_atom("P", md.element.phosphorus, residue)
    protein = topology.add_chain()
    for _ in range(3):
        residue = topology.add_residue("ALA", protein)
        topology.add_atom("CA", md.element.carbon, residue)

    rng = np.random.default_rng(3)
    xyz = rng.uniform(0, 4.0, size=(n_frames, topology.n_atoms, 3))
    traj = md.Trajectory(xyz, topology)
    traj.unitcell_lengths = np.full((n_frames, 3), 4.0)
    traj.unitcell_angles = np.full((n_frames, 3), 90.0)
    pdb = tmp_path / "solvated.pdb"
    dcd = tmp_path / "trajectory.dcd"
    traj[0].save_pdb(str(pdb))
    traj.save_dcd(str(dcd))
    return pdb, dcd


def _analyse_run_module():
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "analyse_run", root / "tools" / "analyse_run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_bilayer_is_written_as_its_own_small_viewer_structure(tmp_path: Path):
    """Section 6.2 asks for the membrane as phosphate spheres on a toggle.

    The viewer bundle drops the solvent and the lipids on load, correctly, on
    the grounds that a 190,000-atom system is mostly water. That reasoning does
    not extend to the phosphorus: beta2AR's bilayer is 476 lipids, so one atom
    each is 476 atoms against about 2,600 for the protein. Written as a
    separate structure rather than merged, so the topology every measurement is
    taken against does not grow non-protein atoms.
    """
    pdb, dcd = _tiny_bilayer(tmp_path)
    out = tmp_path / "analysis"
    out.mkdir()
    bundle = _analyse_run_module().write_lipid_bundle(
        pdb, dcd, "b2ar_apo", out, stride=1)

    assert bundle["n_atoms"] == 8, bundle
    assert (out / "b2ar_apo.lipids.pdb").is_file()
    assert (out / "b2ar_apo.lipids.xtc").is_file()
    # Small enough that the size objection to drawing it does not apply.
    assert bundle["bytes"] < 200_000, bundle


def test_a_water_run_gets_no_lipid_bundle(tmp_path: Path):
    """A system with no bilayer returns nothing, and the toggle is not drawn."""
    pdb, dcd = _tiny_bilayer(tmp_path, lipids=False)
    out = tmp_path / "analysis"
    out.mkdir()
    assert _analyse_run_module().write_lipid_bundle(
        pdb, dcd, "abl1_apo", out, stride=1) == {}
    assert not list(out.glob("*.lipids.*"))


def test_the_run_record_only_names_things_that_exist_where_it_is_written():
    """`run.json` is written at the end of `main`, not inside `build`.

    The force field list is chosen in `build` and recorded in `run.json` by
    `main`, which are different scopes. Referencing `build`'s local from `main`
    parses fine and raises `NameError` at the one moment it costs most: after
    the production loop has finished, before anything is written down. A
    seventeen-hour run would complete and then record nothing.

    Checked by reading the source rather than by running it, because running it
    means running the simulation.
    """
    import ast

    root = Path(__file__).resolve().parents[1]
    tree = ast.parse((root / "tools" / "run_md.py").read_text())
    module_level = set()
    for node in tree.body:
        if isinstance(node, ast.Assign):
            # Walked rather than read directly: several of these are tuple
            # unpackings, `HEAT_NS, EQUIL_NS, PRODUCTION_NS = ...`, whose
            # target is a Tuple and not a Name.
            module_level |= {inner.id for t in node.targets
                             for inner in ast.walk(t)
                             if isinstance(inner, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            module_level.add(node.target.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            module_level |= {(a.asname or a.name).split(".")[0]
                             for a in node.names}
    functions = {node.name: node for node in tree.body
                 if isinstance(node, ast.FunctionDef)}
    main = functions["main"]

    # Every bare name `main` reads that it does not bind itself must come from
    # the module, from a builtin, or from an import.
    bound = {node.id for node in ast.walk(main)
             if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store)}
    bound |= {a.arg for a in main.args.args}
    for node in ast.walk(main):
        if isinstance(node, ast.alias):
            bound.add((node.asname or node.name).split(".")[0])
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            bound.add(node.name)
        if isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)
    import builtins

    free = {node.id for node in ast.walk(main)
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)}
    unexplained = (free - bound - module_level - set(functions)
                   - set(dir(builtins)))
    assert not unexplained, (
        f"`main` reads names that are not module level, not its own, and not "
        f"builtins, so they would raise NameError at run time: "
        f"{sorted(unexplained)}")
    assert "FORCEFIELD_FILES" in module_level


def test_the_thumbnail_centres_the_bilayer_before_drawing_it(tmp_path: Path):
    """Wrapped coordinates put the two leaflets at opposite ends of the box.

    Drawn raw, the receptor appears to be sitting half out of the membrane,
    which is exactly the fault `centre_bilayer` exists for and exactly the kind
    of picture that looks plausible enough not to be questioned. The source is
    checked rather than the pixels: asserting on a PNG would be asserting on
    matplotlib.
    """
    root = Path(__file__).resolve().parents[1]
    source = (root / "tools" / "make_figures.py").read_text()
    thumb = source[source.index("def thumbnail("):source.index("def main(")]
    assert "centre_bilayer" in thumb, (
        "the thumbnail draws wrapped coordinates without centring them")
    assert "unitcell_lengths" in thumb


def test_the_state_strip_starts_on_the_frame_the_slider_starts_on(tmp_path: Path):
    """Section 6.2's badge strip, rendered server-side and then kept in step.

    It is drawn at frame 0 because that is where the scrubber starts. It used
    to be drawn with the run's final label, which disagrees with the viewer
    beside it from the moment the page loads and in a way a reader has no way
    to spot. That needed `start_label` in the summary, which was not there.
    """
    database = tmp_path / "t.db"
    payload = _analysis_payload()
    store.ingest_analysis(database, payload, system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    summary = store.summarise(database, "b2ar_apo")
    for key, row in summary["motifs"].items():
        assert "start_label" in row, key
        assert row["start_label"] is not None, key


def test_every_motif_that_becomes_a_row_has_a_colour_token(tmp_path: Path):
    """The token map was written against a motif list the build then changed.

    Nothing read the tokens, so the mismatch never showed: `c.l` for the
    kinase and `TM5` and `TM7` for the receptor had no colour at all, while
    `r_spine`, `ionic_lock`, `PIF` and `toggle` had one and are not motif rows.
    A motif with no token draws in whatever colour the chart library picks
    next, which is the opposite of a shared palette.
    """
    for family, names in (
            ("kinase", {"A-loop", "alphaC", "c.l", "g.l", "hinge", "xDFG"}),
            ("gpcr", {"H8", "NPxxY", "TM5", "TM6", "TM7"})):
        mapped = set(store.MOTIF_TOKENS[family])
        assert names - mapped == set(), (
            f"{family} motifs with no token: {sorted(names - mapped)}")
        assert mapped - names == set(), (
            f"{family} tokens for things that are not motifs: "
            f"{sorted(mapped - names)}")
        # And no two motifs share a colour, which would make the palette
        # useless for exactly the job it has.
        tokens = list(store.MOTIF_TOKENS[family].values())
        assert len(set(tokens)) == len(tokens), family

    # Section 6.1 fixes these in a table, so every assignment it makes for a
    # motif this build renders has to survive. Only the tokens the spec gave to
    # something that is not a motif row here are free to be reused, and the
    # test says which those are rather than allowing any change at all.
    from_the_spec = {
        "kinase": {"alphaC": "--motif-1", "xDFG": "--motif-2",
                   "A-loop": "--motif-3", "g.l": "--motif-4",
                   "hinge": "--motif-5"},
        "gpcr": {"TM6": "--motif-1", "NPxxY": "--motif-3", "H8": "--motif-6"},
    }
    for family, fixed in from_the_spec.items():
        for motif, token in fixed.items():
            assert store.MOTIF_TOKENS[family][motif] == token, (
                f"section 6.1 gives {family} {motif} {token}, and it is a "
                f"motif this build renders, so it is not free to be moved")


def test_the_state_tokens_section_6_1_names_exist_and_are_used():
    """`--state-active` and its three siblings, named in the same table.

    The first version of the state markers put the four hex values straight
    into the rules, which is the fault the motif tokens had, one paragraph
    later in the same stylesheet.
    """
    root = Path(__file__).resolve().parents[1]
    css = (root / "flexappeal" / "static" / "brand.css").read_text()
    for state in ("active", "inactive", "intermediate", "conflict"):
        assert f"--state-{state}:" in css, f"--state-{state} is not defined"
        assert f"var(--state-{state})" in css, (
            f"--state-{state} is defined and nothing uses it")


def test_the_motif_tokens_resolve_to_real_colours():
    """A token is a CSS custom property name, so the property has to exist.

    These were stored on every motif row from the start and the variables were
    never defined: `var(--motif-1)` resolved to nothing. The charts read them
    through `getComputedStyle`, which returns an empty string for a property
    that is not there and silently falls back.
    """
    root = Path(__file__).resolve().parents[1]
    css = (root / "flexappeal" / "static" / "brand.css").read_text()
    used = {token for family in store.MOTIF_TOKENS.values()
            for token in family.values()}
    for token in sorted(used):
        assert f"{token}:" in css, f"{token} is stored on a motif and undefined"


def test_the_figures_use_the_same_palette_as_the_pages():
    """A motif has to be one colour in the document and on the site.

    `tools/make_figures.py` used to assign colours by index from a list of its
    own, so a motif's colour in a FINDINGS figure had nothing to do with its
    colour on the page. That is a fourth place for one fact to live, and it is
    the fact the token system exists to keep in one place.
    """
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "make_figures", root / "tools" / "make_figures.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    palette = module.motif_palette()
    used = {token for family in store.MOTIF_TOKENS.values()
            for token in family.values()}
    assert used <= set(palette), (
        f"tokens with no colour in the stylesheet: {sorted(used - set(palette))}")
    # Read from the stylesheet, not restated here.
    css = (root / "flexappeal" / "static" / "brand.css").read_text()
    for token, colour in palette.items():
        assert f"{token}:           {colour};" in css or f"{token}:" in css


def _hwprobe():
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "hwprobe", root / "tools" / "hwprobe.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_hardware_probe_writes_the_three_files_section_2_names(tmp_path: Path):
    """Section 2.1 and 2.2 name `hardware.toml`, `tuning.toml` and
    `benchmarks.toml`, and none of the three existed.

    The settings were right, because `run_md.py` works the platform out at run
    time and logs what it kept, but nothing recorded what the machine is or
    what it managed. Section 2.2 ends "every compute budget is derived from
    these measurements, never assumed", and a budget cannot be checked against
    a file that is not there.
    """
    import subprocess
    import sys
    import tomllib

    root = Path(__file__).resolve().parents[1]
    out = tmp_path / "config"
    done = subprocess.run(
        [sys.executable, str(root / "tools" / "hwprobe.py"), "--out", str(out)],
        cwd=root, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr

    for name in ("hardware.toml", "tuning.toml", "benchmarks.toml"):
        path = out / name
        assert path.is_file(), f"{name} was not written"
        tomllib.loads(path.read_text())          # parses, or this raises

    tuning = tomllib.loads((out / "tuning.toml").read_text())
    for key in ("io_workers", "concurrent_md", "memory_headroom_gb"):
        assert key in tuning, key
    assert tuning["concurrent_md"] == 1, "one GPU, one context"
    assert tuning["memory_headroom_gb"] == 24


def test_the_probe_omits_what_it_cannot_measure_rather_than_writing_zero():
    """Section 0: omit a key you cannot produce, never write 0 or null.

    It matters most in a file whose only job is to be trusted later. A machine
    with no `hw.perflevel0.physicalcpu` should get no `cpu_workers`, not a
    `cpu_workers` of 2 that reads like a measurement.
    """
    module = _hwprobe()
    assert module.tuning({}, {}) .get("cpu_workers") is None
    assert "md_platform" not in module.tuning({}, {})
    # And with a measurement it derives the rule: performance cores minus two.
    derived = module.tuning({"performance_cores": 16}, {"OpenCL": 3})
    assert derived["cpu_workers"] == 14
    assert derived["md_platform"] == "OpenCL"
    # Two cores is the floor, so a small machine does not get zero workers.
    assert module.tuning({"performance_cores": 2}, {})["cpu_workers"] == 2


def test_the_toml_writer_refuses_what_it_cannot_represent():
    """It is a small writer, so it raises rather than writing something that
    parses differently from what was meant."""
    module = _hwprobe()
    with pytest.raises(ValueError, match="nests deeper"):
        module.as_toml({"outer": {"inner": {"deeper": 1}}}, header="x")


def test_the_membrane_system_integrates_for_five_picoseconds():
    """Section 10's "membrane builder smoke test on a 5 ps run".

    The one test in section 10's pytest list that did not exist. The others do:
    the numbering maps, the rotated helix, the salt bridge at 3.5 angstrom, the
    state-score symmetry and the resume.

    It runs on the system the membrane builder actually produced rather than on
    a fixture, because the thing being smoke-tested is that builder's output.
    That means it needs the GPU, which for most of this build has been carrying
    a 50 ns row, so it skips while a run is in flight rather than competing for
    the device or spending ten minutes on 194,032 atoms on the CPU. Section 5.3
    already smoke-tested the same system at 1 ns, two hundred times longer, and
    that result is in `run.json`; this is the version that lives in the suite.
    """
    import subprocess

    root = Path(__file__).resolve().parents[1]
    built = root / "data" / "md" / "b2ar_carazolol"
    if not (built / "system.xml").is_file():
        pytest.skip("no built membrane system in this checkout")
    if subprocess.run(["pgrep", "-f", "tools/run_md.py"],
                      capture_output=True).returncode == 0:
        pytest.skip("a production run holds the GPU; this needs the device")

    import openmm as mm
    from openmm import app, unit

    pdb = app.PDBFile(str(built / "solvated.pdb"))
    system = mm.XmlSerializer.deserialize((built / "system.xml").read_text())
    integrator = mm.LangevinMiddleIntegrator(
        303.15 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=200)
    simulation.step(2500)                       # 5 ps at 2 fs

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    assert energy == energy, "the potential energy is NaN after 5 ps"
    assert energy < 0, f"a solvated bilayer should be bound, not {energy:.3g}"
    coordinates = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    assert np.isfinite(coordinates).all(), "a coordinate went non-finite"


def _ensemble_of(lengths: dict[str, list[tuple[str, int]]]):
    """Ensembles with cores of given lengths, for the alignment target tests."""
    rng = np.random.default_rng(0)
    out = {}
    for state, items in lengths.items():
        found = ensembles.MotifEnsemble(state=state)
        for code, n in items:
            found.by_code[code] = {ensembles.CORE_KEY: rng.normal(size=(n, 3)),
                                   "motif": rng.normal(size=(4, 3))}
        out[state] = found
    return out


def test_one_odd_reference_does_not_discard_the_whole_ensemble():
    """The alignment target used to be whichever reference sorted first.

    Every other reference is dropped when its core length differs from the
    target's, so a target that is one residue shorter than everyone else takes
    the entire ensemble with it. On CDK2 that is exactly what happened: the
    target was 1F5Q with a 41-residue core, 28 of the other 29 references had
    42, and one reference survived out of twenty-nine. The score that came out
    of it was a distance to a single structure.

    The target is now the first, in sorted order, whose core is the most common
    length in the set.
    """
    found = _ensemble_of({
        # "1F5Q" sorts before every "A..." code, as it did in the real set.
        "active": [("1F5Q", 41)] + [(f"A{i:03d}", 42) for i in range(28)],
        "inactive": [("Z001", 42)],
    })
    report = ensembles.align_to_common_frame(found)
    assert report["n_aligned"] == 29, report
    assert report["n_dropped"] == 1, "only the odd one out should go"
    assert report["target"] != "1F5Q", (
        "the odd reference was chosen as the target again")


def test_a_uniform_ensemble_aligns_all_of_itself():
    """Was: a uniform set keeps the first reference as its target.

    That assertion encoded D025's tiebreak, which was sorted order, and it was
    the right thing to pin while modal length was doing the deciding. Once
    `traj.intersect_core` made every core the same length the modal rule stopped
    discriminating at all and sorted order became the whole rule, which on CDK2
    picked 1F5Q and put the mean residual at 0.94 A against 0.40 onto the
    medoid. So the target is the medoid now and the identity of it is not the
    property worth pinning.

    What is worth pinning is that a uniform set loses nobody, which is the
    thing D025 was protecting. The fixture's cores are generated from one shape,
    so every reference is equidistant and the medoid is decided by the sorted
    tiebreak, which keeps this deterministic without asserting on it.
    """
    found = _ensemble_of({"active": [(f"B{i:03d}", 43) for i in range(5)],
                          "inactive": [(f"C{i:03d}", 43) for i in range(14)]})
    report = ensembles.align_to_common_frame(found)
    assert report["n_aligned"] == 19
    assert report["n_dropped"] == 0
    assert report["target"] in {f"B{i:03d}" for i in range(5)} | {
        f"C{i:03d}" for i in range(14)}


def _payload_with_frames(n: int = 12) -> dict:
    """The shared fixture has two frames, which is below the four `summarise`
    needs before it will report halves at all. That is correct behaviour and
    makes it the wrong fixture for testing them."""
    payload = _analysis_payload()
    payload["frames"] = [
        {"time_ns": i * 0.1, "n_conflicts": 0, "motifs": [
            {"motif": "xDFG", "d_act_a": 1.0, "d_inact_a": 4.0,
             "s": 0.6 + 0.01 * i, "label": "active"},
            {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0,
             "s": 0.0, "label": "intermediate"}]}
        for i in range(n)]
    return payload


def test_the_summary_counts_frames_and_not_rows(tmp_path: Path):
    """`motif_frames` holds one row per metric, and a motif has two.

    Counting rows reports twice as many frames as the run has, which is the
    same rows-against-things error section 4m records for the kinase reference
    ensembles. The means and the scatter survive the duplication; the count
    does not, and it is the number a reader would quote.
    """
    database = tmp_path / "t.db"
    payload = _payload_with_frames(12)
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")
    counted = {motif: values.get("n_scored_frames")
               for motif, values in summary["motifs"].items()}
    assert counted, "no motif was summarised"
    for motif, n in counted.items():
        assert n == 12, f"{motif} counted {n} where the payload has 12 frames"


def test_the_summary_reports_the_halves_and_the_scatter(tmp_path: Path):
    """A start and a final cannot tell a motif that moved from a noisy one.

    Section 5.0b is the case: beta2AR's H8 falls 0.176 between the halves with
    a standard deviation of 0.168 behind it, and reading only the endpoints
    would have called that a drift. Both numbers live in `summarise` so the
    page and the report cannot derive them differently.
    """
    database = tmp_path / "t.db"
    store.ingest_analysis(database, _payload_with_frames(12), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")

    # xDFG rises steadily, so its second half must read higher than its first.
    rising = summary["motifs"]["abl1:xDFG"]
    assert rising["second_half_mean_s"] > rising["first_half_mean_s"]
    assert rising["sd_s"] > 0.0

    # hinge is constant, so it has no scatter and no difference at all.
    flat = summary["motifs"]["abl1:hinge"]
    assert flat["sd_s"] == 0.0
    assert flat["first_half_mean_s"] == flat["second_half_mean_s"]


def test_every_motif_keeps_its_final_tenth_mean(tmp_path: Path):
    """The trend work shadowed the name holding the final-tenth time cut.

    `cut` is computed once before the motif loop as a time in nanoseconds, and
    the halves calculation reused the name for a frame count. The first motif
    was summarised before the overwrite and looked right; every one after it
    asked the database for frames past 73 ns of a 14.6 ns run and got NULL. The
    page still rendered a number for the first motif and then raised a 500 on
    the second, which is the only reason it was noticed at all.

    One motif having a value is not evidence, so this asserts on all of them.
    """
    database = tmp_path / "t.db"
    store.ingest_analysis(database, _payload_with_frames(12), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")
    assert len(summary["motifs"]) >= 2, (
        "a single motif cannot show an error that starts at the second")
    for motif, values in summary["motifs"].items():
        assert values["final_tenth_mean_s"] is not None, (
            f"{motif} lost its final-tenth mean, which is what happens when "
            f"the time cut is overwritten inside the loop")


def test_a_run_too_short_to_halve_reports_no_trend(tmp_path: Path):
    """Two frames cannot give a meaningful half-mean, so none is offered.

    Section 0's rule: omit a key you cannot produce rather than writing a
    number that reads like a measurement.
    """
    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")
    for motif, values in summary["motifs"].items():
        assert "sd_s" not in values, motif
        assert values.get("final_s") is not None, (
            f"{motif} lost its endpoints as well")


def test_the_summary_agrees_with_the_payload_it_was_ingested_from(tmp_path: Path):
    """Two independent paths to the same seven numbers must agree.

    The analysis JSON is what a run writes; the SQLite summary is what the
    pages read. Everything in between is ingest and SQL, and a fault there
    produces numbers that are wrong rather than absent, which no amount of
    rendering will reveal.

    This is the shape of check that would have caught the final-tenth
    shadowing in section 5.0b on the first run rather than when the page
    raised a 500, and unlike that 500 it would also have caught the version
    where the overwritten value still formats.
    """
    import statistics

    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")

    frames = payload["frames"]
    span_hi, span_lo = frames[-1]["time_ns"], frames[0]["time_ns"]
    cut = span_hi - 0.1 * (span_hi - span_lo)
    series: dict[str, list] = {}
    tail: dict[str, list] = {}
    for frame in frames:
        for motif in frame["motifs"]:
            series.setdefault(motif["motif"], []).append(motif["s"])
            if frame["time_ns"] >= cut:
                tail.setdefault(motif["motif"], []).append(motif["s"])

    assert summary["motifs"], "nothing was summarised"
    for key, values in summary["motifs"].items():
        name = key.split(":")[-1]
        mine = series[name]
        half = len(mine) // 2
        expected = {
            "start_s": mine[0],
            "final_s": mine[-1],
            "final_tenth_mean_s": statistics.mean(tail[name]),
            "first_half_mean_s": round(statistics.mean(mine[:half]), 4),
            "second_half_mean_s": round(statistics.mean(mine[half:]), 4),
            "sd_s": round(statistics.pstdev(mine), 4),
            "n_scored_frames": len(mine),
        }
        for field, want in expected.items():
            got = values.get(field)
            assert got is not None, f"{name}.{field} is missing"
            assert abs(float(got) - float(want)) < 5e-4, (
                f"{name}.{field}: the database says {got} and the payload "
                f"says {want}")


def test_verify_summary_catches_a_number_that_drifted(tmp_path: Path):
    """The live version of the round-trip check, for the unattended landing.

    The fixture test above proves the two paths agree on a clean ingest. This
    proves the checker would actually complain, because a verifier that cannot
    fail is worth nothing and the final pass runs this one on four rows nobody
    is watching.
    """
    import json as _json

    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")

    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (analysis / "abl1_apo.json").write_text(_json.dumps(payload))
    assert _cli().verify_summary(database, analysis) == [], (
        "a clean ingest should agree with its own payload")

    # Move one score in the file and the database must now disagree with it.
    payload["frames"][0]["motifs"][0]["s"] = 0.123456
    (analysis / "abl1_apo.json").write_text(_json.dumps(payload))
    problems = _cli().verify_summary(database, analysis)
    assert problems, "a changed number went unnoticed"
    assert any("start_s" in complaint for complaint in problems), problems


def test_verify_summary_is_quiet_about_a_run_with_no_analysis_file(tmp_path: Path):
    """A row ingested without its payload beside it is not a disagreement.

    The analysis directory is a working directory; a missing file there means
    there is nothing to compare against, which is different from the two
    sources saying different things.
    """
    database = tmp_path / "t.db"
    store.ingest_analysis(database, _payload_with_frames(20), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    empty = tmp_path / "nothing"
    empty.mkdir()
    assert _cli().verify_summary(database, empty) == []


def test_a_checkpoint_is_bound_to_the_platform_that_wrote_it():
    """OpenMM refuses a checkpoint written on another platform.

    Found by trying to read `abl1_apo`'s checkpoint on the CPU while the GPU
    was busy: "Checkpoint was created with a different Platform: OpenCL". It
    matters because `run_md.platform()` falls back to CPU when OpenCL is
    unavailable, so a GPU that disappears between a row starting and its resume
    turns into an OpenMM traceback at the top of a 24-hour queue, with the
    row's progress still on disk and nothing saying it is still good.
    """
    import openmm as mm

    source = _argon_box(with_barostat=False)
    reference = mm.Platform.getPlatformByName("Reference")
    cpu = mm.Platform.getPlatformByName("CPU")

    written = mm.Context(
        source, mm.VerletIntegrator(0.001), reference)
    written.setPositions([mm.Vec3(0.3 * (i % 2), 0.0, 0.0) for i in range(8)])
    checkpoint = written.createCheckpoint()

    other = mm.Context(_argon_box(with_barostat=False),
                       mm.VerletIntegrator(0.001), cpu)
    with pytest.raises(Exception, match="Platform"):
        other.loadCheckpoint(checkpoint)


def test_the_resume_explains_a_platform_mismatch_rather_than_raising():
    """The runner has to say what happened, and not redo the work.

    Starting the row again would discard the nanoseconds in the checkpoint
    without being asked. Section 0's ladder is about degrading with the reason
    recorded, so this checks the source for the guard and for the refusal to
    silently restart.
    """
    root = Path(__file__).resolve().parents[1]
    source = (root / "tools" / "run_md.py").read_text()
    resume = source[source.index("    if resuming:"):
                    source.index("    elif args.membrane:")]
    assert "try:" in resume and "loadCheckpoint" in resume, (
        "the checkpoint load is unguarded")
    assert "SystemExit" in resume, (
        "a failed load must stop the row rather than fall through")
    assert "still" in resume and "good" in resume, (
        "the message should say the checkpoint is not the thing that broke")


# ===========================================================================
#  Resuming without duplicating frames
# ===========================================================================


def _run_md():
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "tools" / "run_md.py"
    spec = importlib.util.spec_from_file_location("_run_md", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_checkpoint_interval_matches_the_frame_interval():
    """Otherwise a resume rewrites every frame since the last checkpoint.

    With a 1 ns checkpoint and 100 ps frames the checkpoint sits up to ten
    frames behind, and because frame time here is the index times the interval,
    duplicates shift every later timestamp. abl1_apo had exactly this: 43
    frames against a checkpoint at 4.0 ns.
    """
    module = _run_md()
    assert module.CHECKPOINT_NS == pytest.approx(module.FRAME_PS / 1000.0)


def test_a_resume_truncates_frames_written_past_the_checkpoint(tmp_path: Path):
    """The frames kept must be the originals, and the box must survive."""
    import mdtraj as md
    import numpy as np
    from mdtraj.formats import DCDTrajectoryFile

    module = _run_md()
    path = tmp_path / "trajectory.dcd"
    rng = np.random.RandomState(0)
    xyz = rng.rand(12, 5, 3).astype("float32") * 10.0
    lengths = np.full((12, 3), 30.0, dtype="float32")
    angles = np.full((12, 3), 90.0, dtype="float32")
    with DCDTrajectoryFile(str(path), "w") as handle:
        handle.write(xyz, cell_lengths=lengths, cell_angles=angles)

    dropped = module._truncate_trajectory(path, 8)
    assert dropped == 4
    with DCDTrajectoryFile(str(path)) as handle:
        kept, kept_lengths, _ = handle.read()
    assert kept.shape[0] == 8
    assert np.allclose(kept, xyz[:8], atol=1e-3)
    assert np.allclose(kept_lengths, lengths[:8], atol=1e-3)

    # And it does nothing when the checkpoint is not behind.
    assert module._truncate_trajectory(path, 8) == 0
    assert module._truncate_trajectory(path, 99) == 0


def test_the_state_table_is_truncated_with_the_trajectory(tmp_path: Path):
    """The two have to stay the same length or every later analysis is offset."""
    module = _run_md()
    path = tmp_path / "state_data.csv"
    header = '#"Step","Time (ps)","Total Energy (kJ/mole)"'
    rows = [f"{i * 25000},{i * 100.0},{-1000.0 - i}" for i in range(1, 13)]
    path.write_text("\n".join([header] + rows) + "\n")

    assert module._truncate_state_data(path, 8) == 4
    lines = path.read_text().splitlines()
    assert lines[0] == header
    assert len(lines) - 1 == 8
    assert lines[-1].startswith("200000,800.0")
    assert module._truncate_state_data(path, 8) == 0


def test_truncation_is_safe_when_there_is_nothing_to_truncate(tmp_path: Path):
    module = _run_md()
    assert module._truncate_trajectory(tmp_path / "absent.dcd", 5) == 0
    assert module._truncate_state_data(tmp_path / "absent.csv", 5) == 0
    empty = tmp_path / "empty.csv"
    empty.write_text("")
    assert module._truncate_state_data(empty, 5) == 0


# ===========================================================================
#  Helix axes, and a measurement that was a parameter choice
# ===========================================================================


def test_a_helix_axis_is_recovered_from_a_built_helix():
    """An ideal alpha helix along z must report an axis along z.

    Built from the standard geometry, 2.3 Å radius, 1.5 Å rise and 100 degrees
    per residue, so the answer is known rather than looked up.
    """
    import numpy as np

    from flexappeal.refstates import traj as T

    angles = np.radians(np.arange(18) * 100.0)
    helix = np.column_stack([2.3 * np.cos(angles), 2.3 * np.sin(angles),
                             1.5 * np.arange(18)])
    axis = T.helix_axis(helix)
    assert axis is not None
    assert abs(float(axis[2])) > 0.999, f"axis is not along z: {axis}"


def test_turn_averaging_is_what_makes_a_short_helix_axis_stable():
    """The measurement that showed one window was a parameter choice.

    A principal axis fitted straight to Calpha wanders as the window changes,
    because the Calpha spiral. Turn-averaging removes the spiral. Asserted as a
    comparison between the two, because the point is not that one works but
    that the other does not.
    """
    import numpy as np

    from flexappeal.refstates import traj as T

    angles = np.radians(np.arange(24) * 100.0)
    helix = np.column_stack([2.3 * np.cos(angles), 2.3 * np.sin(angles),
                             1.5 * np.arange(24)])

    def raw_axis(points):
        centred = points - points.mean(axis=0)
        _, _, vectors = np.linalg.svd(centred, full_matrices=False)
        return vectors[0]

    from flexappeal.refstates.gpcr import KINK_MIN_AXIS_POINTS

    raw, averaged = [], []
    for window in range(5, 16):
        segment = helix[:window]
        raw.append(abs(float(raw_axis(segment)[2])))
        # Only windows with enough axis points to fit a line through, which is
        # the condition the production code enforces rather than a convenience
        # for this test: below it, turn-averaging is fitting two points.
        if T.helix_axis_points(segment).shape[0] < KINK_MIN_AXIS_POINTS:
            continue
        found = T.helix_axis(segment)
        if found is not None:
            averaged.append(abs(float(found[2])))

    assert len(averaged) >= 4, "not enough usable windows to judge stability"
    # Turn-averaging holds the axis on z once there are enough axis points; the
    # raw fit does not, at any window, and the spread is the evidence.
    assert max(averaged) - min(averaged) < 0.005, averaged
    assert max(raw) - min(raw) > 10 * (max(averaged) - min(averaged)), (
        f"raw spread {max(raw) - min(raw):.4f} against averaged "
        f"{max(averaged) - min(averaged):.4f}")


def test_an_unresolved_kink_is_reported_as_unresolved():
    """A spread wider than half the value is not a measurement.

    The active start model's P7.50 reads 7.4 degrees with a 17.2 degree spread
    across windows, so the windows disagree about the answer and the flag says
    so rather than letting the median stand alone.
    """
    import gemmi

    from flexappeal.refstates import Catalogue, receptor_residues
    from flexappeal.refstates.gpcr import build_numbering, kink_angles

    start = Path(__file__).resolve().parents[1] / "data" / "systems" / "b2ar_apo_start.pdb"
    if not start.is_file():
        pytest.skip("the receptor start model is not built in this checkout")
    cache = Path(__file__).resolve().parents[1] / "data" / "cache" / "gpcrdb"
    if not cache.is_dir():
        pytest.skip("no GPCRdb cache; this test does not reach the network")

    numbering = build_numbering(receptor_residues(Catalogue(cache)))
    structure = gemmi.read_structure(str(start))
    structure.setup_entities()
    found = kink_angles(structure[0]["A"], numbering)

    for generic, row in found.items():
        if row.get("status") != "ok":
            continue
        assert "resolved" in row, generic
        if row["resolved"]:
            assert row["spread_deg"] < 0.5 * row["kink_deg"]


def test_a_structure_compared_with_itself_reports_no_motion():
    """The invariant that caught a wrong frame of reference.

    Passing the scoring transform where a start-model fit was wanted gave a 78
    angstrom centroid shift and an 82 angstrom connector RMSD. Both are
    obviously wrong, which was luck; this asserts the version that is not a
    matter of luck, that a structure against itself is exactly zero.
    """
    import gemmi

    from flexappeal.refstates import Catalogue, receptor_residues
    from flexappeal.refstates import gpcr as G
    from flexappeal.refstates import kinase as K
    from flexappeal.refstates.gpcr import build_numbering
    from flexappeal.refstates.traj import gpcr_core

    root = Path(__file__).resolve().parents[1]
    start = root / "data" / "systems" / "b2ar_apo_start.pdb"
    cache = root / "data" / "cache" / "gpcrdb"
    if not start.is_file() or not cache.is_dir():
        pytest.skip("the receptor start model or its catalogue cache is absent")

    numbering = build_numbering(receptor_residues(Catalogue(cache)))
    structure = gemmi.read_structure(str(start))
    structure.setup_entities()
    chain = structure[0]["A"]
    core = gpcr_core(numbering, {r.seqid.num for r in chain})

    fit = K.core_fit(chain, chain, core.residues)
    assert fit is not None

    rotation = G.tm6_rotation(chain, chain, numbering, fit)
    assert rotation["centroid_shift_a"] == pytest.approx(0.0, abs=1e-6)
    assert rotation["axis_rotation_deg"] == pytest.approx(0.0, abs=1e-6)

    connector = G.connector_rmsd(chain, chain, numbering, fit)
    assert connector["rmsd_a"] == pytest.approx(0.0, abs=1e-6)


def test_the_core_fit_is_required_for_a_helix_measurement_to_mean_anything():
    """Without it the number includes whole-protein drift.

    Built rather than loaded: one helix, rotated and translated bodily. With a
    core fit on the rest of the structure the helix has not moved relative to
    it, so the measurement must be near zero; without one it reports the bodily
    move.
    """
    import numpy as np

    from flexappeal.refstates import traj as T

    rng = np.random.RandomState(3)
    core = rng.rand(20, 3) * 20.0
    helix = np.column_stack([
        2.3 * np.cos(np.radians(np.arange(14) * 100.0)),
        2.3 * np.sin(np.radians(np.arange(14) * 100.0)),
        1.5 * np.arange(14)]) + np.array([40.0, 0.0, 0.0])

    # About y, not z. The helix axis is along z, so a rotation about z leaves
    # the axis where it was and the test would prove nothing: the first version
    # of this test did exactly that and measured 0.37 degrees.
    angle = np.radians(35.0)
    spin = np.array([[np.cos(angle), 0.0, np.sin(angle)],
                     [0.0, 1.0, 0.0],
                     [-np.sin(angle), 0.0, np.cos(angle)]])
    shift = np.array([11.0, -4.0, 6.0])
    moved_core, moved_helix = core @ spin.T + shift, helix @ spin.T + shift

    transform = T.superpose(moved_core, core)
    fitted = T.apply_transform(moved_helix, *transform)
    assert T.rmsd(fitted, helix) < 1e-6, "the core fit should undo the move"

    # And the axis angle is zero once fitted, where the unfitted axes differ.
    assert T.axis_angle(T.helix_axis(fitted), T.helix_axis(helix)) == pytest.approx(
        0.0, abs=1e-4)
    assert T.axis_angle(T.helix_axis(moved_helix), T.helix_axis(helix)) > 30.0


def test_smoothing_is_skipped_when_the_window_is_wider_than_the_run():
    """A 2 ns centred window on a 1 ns run is one value repeated.

    Which draws as a flat trace, and a flat trace is a claim about the
    trajectory rather than about the window. Asserted on the smoothing function
    itself: over a run shorter than the window every point collapses to the
    same mean, which is the behaviour the caller has to refuse to publish.
    """
    values = [0.1, 0.3, 0.5, 0.7, 0.9]
    # 5 frames at 100 ps is 0.5 ns, against a 2 ns window.
    collapsed = state.smooth(values, frame_ps=100.0, window_ns=2.0)
    assert len(set(round(v, 6) for v in collapsed)) == 1, collapsed
    assert collapsed[0] == pytest.approx(sum(values) / len(values))

    # And over a run comfortably longer than the window it still varies.
    long_run = [i / 100.0 for i in range(100)]
    smoothed = state.smooth(long_run, frame_ps=100.0, window_ns=2.0)
    assert len(set(round(v, 6) for v in smoothed)) > 50


def test_an_sse_lean_on_a_tiny_margin_is_not_called():
    """Secondary structure is far less discriminating than geometry.

    Most residues of a transmembrane helix are helical in both states, so every
    agreement is high and the differences are small. H8 came out 0.864 active
    against 0.888 inactive, and calling that "inactive" on a 0.024 margin is
    reading noise.
    """
    # A frame that matches the active set slightly better, but only slightly.
    frame = ["H"] * 20
    active = sse.reference_profile({
        "a": ["H"] * 20,
        "b": ["H"] * 19 + ["T"],
    })
    inactive = sse.reference_profile({
        "c": ["H"] * 20,
        "d": ["H"] * 18 + ["T", "T"],
    })
    found = sse.summarise(sse.agreement(frame, active, inactive))
    assert found["margin"] < sse.DECISIVE_MARGIN
    assert found["decisive"] is False
    assert found["leans"] == "neither"

    # And a real difference is still called.
    inactive_far = sse.reference_profile({"c": ["E"] * 20, "d": ["E"] * 20})
    decided = sse.summarise(sse.agreement(frame, active, inactive_far))
    assert decided["decisive"] is True
    assert decided["leans"] == "active"


def test_sse_is_scored_only_where_the_reference_sets_disagree():
    """The fault that gave two different structures the same answer.

    A transmembrane helix is helical in both states, so its SSE string is the
    same in each and so is its agreement. Scored across all residues, TM6 read
    0.934 active against 0.868 inactive for a structure built from the active
    reference *and* for one built from the inactive: bit-identical, and the lean
    was reporting that the active reference set is more uniformly helical, which
    is a fact about the catalogue rather than about the frame.
    """
    # Every reference helical everywhere: nothing here can tell the states
    # apart, however lopsided the agreement looks.
    active = sse.reference_profile({"a": ["H"] * 10, "b": ["H"] * 10})
    inactive = sse.reference_profile({"c": ["H"] * 9 + ["T"],
                                      "d": ["H"] * 9 + ["T"]})
    # The consensus differs at exactly one residue, below the minimum.
    assert sse.discriminating_residues(active, inactive) == {9}
    found = sse.motif_agreement(["H"] * 10, "m",
                                {"active": {"m": active},
                                 "inactive": {"m": inactive}})
    assert found["leans"] == "neither"
    assert found["decisive"] is False
    assert "too few" in found["note"]

    # No disagreement at all is reported as carrying no information.
    same = sse.reference_profile({"c": ["H"] * 10, "d": ["H"] * 10})
    none = sse.motif_agreement(["H"] * 10, "m",
                               {"active": {"m": active}, "inactive": {"m": same}})
    assert none["n_discriminating"] == 0
    assert "no state information" in none["note"]


def test_sse_tells_two_states_apart_when_enough_residues_differ():
    """And the same measurement must give opposite answers for opposite frames.

    This is the property the all-residue version did not have: it is not enough
    that a frame matching the active set reads active, the frame matching the
    inactive set has to read inactive.
    """
    # Five residues where the sets genuinely disagree.
    active = sse.reference_profile({f"a{i}": ["H"] * 5 for i in range(3)})
    inactive = sse.reference_profile({f"i{i}": ["T"] * 5 for i in range(3)})
    profiles = {"active": {"m": active}, "inactive": {"m": inactive}}
    assert len(sse.discriminating_residues(active, inactive)) == 5

    looks_active = sse.motif_agreement(["H"] * 5, "m", profiles)
    looks_inactive = sse.motif_agreement(["T"] * 5, "m", profiles)
    assert looks_active["leans"] == "active"
    assert looks_inactive["leans"] == "inactive"
    assert looks_active["decisive"] and looks_inactive["decisive"]


# --- the reference manifests, and the script that rebuilds them --------------
#
# `tools/build_references.py` exists so the committed manifests can be
# re-derived rather than trusted. These tests cover the part of it that does
# not need the network: the invariant it has to preserve, and the diff it
# reports. The invariant is worth a test because getting it wrong is silent --
# a per-chain rebuild still passes P1 and P2, it just counts a different n.

def _build_references():
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_build_references", root / "tools" / "build_references.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_kinase_manifest_holds_one_row_per_deposition():
    """KLIFS lists a row per chain; the manifest must not.

    ABL1 comes back as 61 KLIFS rows over 20 depositions, with 2G1T appearing
    four times. Chains of one asymmetric unit are not independent checks of
    whether our labels agree with KLIFS's, so a duplicate PDB code here would
    mean P1 and P2 were quoting an n larger than what was actually measured.
    """
    import json
    from collections import Counter
    from pathlib import Path

    path = (Path(__file__).resolve().parents[1]
            / "data" / "manifests" / "refs_abl1.jsonl")
    rows = [json.loads(line) for line in path.read_text().splitlines()
            if line.strip()]
    assert rows, "the kinase manifest is empty"
    repeated = {code: n for code, n in
                Counter(r["pdb"] for r in rows).items() if n > 1}
    assert not repeated, f"these depositions appear more than once: {repeated}"


def test_manifest_diff_is_per_record_not_per_file():
    """A reordering is not a change.

    The order a catalogue returns its rows in is not stable, so a diff taken on
    the file would report a shuffle as a difference and the verify mode would
    cry wolf on every run.
    """
    module = _build_references()
    one = [{"pdb": "1ABC", "chain": "A", "our_dfg": "in"},
           {"pdb": "2DEF", "chain": "B", "our_dfg": "out"}]
    shuffled = list(reversed(one))
    assert not module.compare(shuffled, one)["changed"]

    changed = [dict(one[0], our_dfg="out"), one[1]]
    diff = module.compare(changed, one)
    assert [row["record"] for row in diff["changed"]] == [("1ABC", "A")]
    assert diff["changed"][0]["fields"]["our_dfg"] == ("in", "out")

    diff = module.compare(one[:1], one)
    assert diff["removed"] == [("2DEF", "B")] and not diff["added"]
    diff = module.compare(one, one[:1])
    assert diff["added"] == [("2DEF", "B")] and not diff["removed"]


def test_manifest_diff_tolerates_float_noise():
    """A distance that moves in the sixth decimal is not a changed reference."""
    module = _build_references()
    committed = [{"pdb": "1ABC", "chain": "A", "tm3_tm6": 11.3900}]
    jittered = [{"pdb": "1ABC", "chain": "A", "tm3_tm6": 11.3901}]
    assert not module.compare(jittered, committed)["changed"]
    moved = [{"pdb": "1ABC", "chain": "A", "tm3_tm6": 11.5}]
    assert module.compare(moved, committed)["changed"]


def test_row_adapter_exposes_what_agreement_reads():
    """`refs.agreement` takes prepared references; verify re-runs it on dicts.

    Rebuilding PreparedReference objects to re-run a control would mean
    fetching and annotating all 20 structures a second time.
    """
    from flexappeal.refstates import refs

    module = _build_references()
    rows = [{"status": "ok", "dfg_agrees": True, "alphac_agrees": True},
            {"status": "ok", "dfg_agrees": False, "alphac_agrees": True}]
    control = refs.agreement([module._Row(r) for r in rows], "dfg")
    assert control["compared"] == 2 and control["agreed"] == 1


def test_rspine_verdict_reads_the_gap_the_way_the_state_rule_does():
    """D015's threshold check, on numbers rather than on 20 structures.

    This is the half of `rspine_separation` that decides whether the kinase
    verify passes, so it is the half worth testing directly: a verdict that
    inverted would still print a tidy table.
    """
    module = _build_references()

    # The real measurement: active tops out at 5.96, inactive starts at 10.35.
    measured = {"active": [5.33, 5.61, 5.66, 5.92, 5.96],
                "inactive": [10.35, 10.72, 11.11, 12.92]}
    out = module.spine_verdict(measured, 7.7, 6.0)
    assert out["gap_a"] == [5.96, 10.35]
    assert out["gap_width_a"] == 4.39
    assert out["threshold_in_gap"] and out["correct"]
    assert out["n_classified"] == 9
    # The published value also classifies these crystals, by four hundredths,
    # which is the whole argument of D015 rather than an incidental number.
    assert out["published_in_gap"]
    assert out["published_margin_a"] == 0.04

    # A threshold below the active cluster calls an active reference broken.
    low = module.spine_verdict(measured, 5.5, 6.0)
    assert not low["threshold_in_gap"] and not low["correct"]

    # And one above the inactive cluster calls an inactive reference assembled.
    high = module.spine_verdict(measured, 11.0, 6.0)
    assert not high["threshold_in_gap"] and not high["correct"]


def test_rspine_verdict_stays_silent_without_both_clusters():
    """No opinion rather than a default one.

    A run where the structures failed to load must not report a gap it did not
    measure, for the same reason the manifests omit a key they cannot produce.
    """
    module = _build_references()
    out = module.spine_verdict({"active": [5.3, 5.9]}, 7.7, 6.0)
    assert "gap_a" not in out and "correct" not in out
    assert out["threshold_a"] == 7.7


# --- an empty result that is a fault, told apart from one that is not --------
#
# Two faults this build hit were the same shape: a function returned an empty
# result because a precondition was missing, and the caller could not tell that
# from "there was nothing to do". The bilayer QC is the case that matters, since
# it decides whether P6 has anything to say.

def _ingest_module():
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_ingest_for_tests", root / "tools" / "ingest.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bilayer_qc_refuses_audibly_when_a_membrane_row_lacks_its_topology(
        tmp_path: Path, capsys):
    """A POPC row with no `solvated.pdb` beside its trajectory is a fault.

    Before this said anything, the pass returned an empty list and the only
    visible consequence was P6 reading `pending` on a finished membrane row,
    with nothing anywhere explaining it.
    """
    module = _ingest_module()
    run_dir = tmp_path / "b2ar_apo"
    run_dir.mkdir()
    (run_dir / "trajectory.dcd").write_bytes(b"not really a dcd")

    assert module.membrane_chunks(run_dir, "b2ar_apo",
                                  expect_bilayer=True) == []
    said = capsys.readouterr().out
    assert "REFUSING" in said and "solvated.pdb" in said


# --- an atom count belongs to a box, not to a protein -----------------------
#
# Section 7 gives `systems` one `n_atoms` column and section 6.2 asks the
# Examples card to show an atom count. A system here has two runs, and their
# boxes are not the same size: imatinib adds 1,464 atoms to ABL1 and carazolol
# removes 6,184 from the receptor, because a bound ligand changes how many
# waters and lipids fit. The column used to be filled by whichever row was
# ingested last, and the card printed that one row's box as the system's.

def test_the_atom_count_of_every_row_of_a_system_is_gathered(tmp_path: Path):
    """Both rows, from their own manifests, rather than the one being ingested."""
    module = _ingest_module()
    for name, atoms in (("abl1_apo", 46721), ("abl1_imatinib", 48185)):
        run = tmp_path / name
        run.mkdir()
        (run / "run.json").write_text(json.dumps({"atoms": atoms}))

    assert module.system_atom_counts(tmp_path, "abl1") == {
        "abl1_apo": 46721, "abl1_imatinib": 48185}


def test_a_system_whose_rows_disagree_shows_the_range_and_stores_no_count(
        tmp_path: Path):
    """And one whose rows agree keeps the single number in the column."""
    database = tmp_path / "t.db"
    store.create(database)

    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water",
                          system_fields={"name": "ABL1 kinase domain",
                                         "n_atoms": None},
                          atoms_by_run={"abl1_apo": 46721,
                                        "abl1_imatinib": 48185})
    with store.connect(database) as connection:
        row = dict(connection.execute(
            "SELECT * FROM systems WHERE system_id = 'abl1'").fetchone())

    assert row["n_atoms"] is None, (
        "no single number is true of a system whose rows differ, and the "
        "column should say so rather than carry one of them")
    # The card that printed the range is gone with the fifth tab. What the
    # store does with two disagreeing rows is still the claim worth holding:
    # the per-run counts stay in `notes` and the column says nothing.
    assert json.loads(row["notes"])["atoms"] == {"abl1_apo": 46721,
                                                 "abl1_imatinib": 48185}

    store.ingest_analysis(database, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr",
                          environment="membrane",
                          system_fields={"n_atoms": 194032},
                          atoms_by_run={"b2ar_apo": 194032,
                                        "b2ar_replica": 194032})
    with store.connect(database) as connection:
        agreed = dict(connection.execute(
            "SELECT * FROM systems WHERE system_id = 'b2ar'").fetchone())
    assert agreed["n_atoms"] == 194032
    assert json.loads(agreed["notes"])["atoms"] == {"b2ar_apo": 194032,
                                                    "b2ar_replica": 194032}


def test_a_stale_atom_count_is_cleared_rather_than_left_behind(tmp_path: Path):
    """The regression: an upsert leaves out what it is not given.

    Omitting the key where the rows disagree left the column holding the number
    from when they agreed, which reads exactly like a measurement.
    """
    database = tmp_path / "t.db"
    store.create(database)
    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water",
                          system_fields={"n_atoms": 46721},
                          atoms_by_run={"abl1_apo": 46721})
    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_imatinib", family="kinase",
                          environment="water",
                          system_fields={"n_atoms": None},
                          atoms_by_run={"abl1_apo": 46721,
                                        "abl1_imatinib": 48185})
    with store.connect(database) as connection:
        count = connection.execute(
            "SELECT n_atoms FROM systems WHERE system_id = 'abl1'").fetchone()[0]
    assert count is None


def test_a_system_is_named_for_itself_rather_than_for_its_last_run(
        tmp_path: Path):
    """`systems.name` used to be the analysis payload's name, which is a run id.

    Nothing on the Examples card showed it, because `SYSTEM_TITLES` has a title
    for both of this build's systems, but the card falls back to this column
    for a system it does not know, which is every system the reuse path brings
    in, and the verifier prints it beside each system it checks.
    """
    database = tmp_path / "t.db"
    store.create(database)
    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_imatinib", family="kinase",
                          environment="water",
                          system_fields={"name": "ABL1 kinase domain"})
    with store.connect(database) as connection:
        name = connection.execute(
            "SELECT name FROM systems WHERE system_id = 'abl1'").fetchone()[0]
    assert name == "ABL1 kinase domain"


def test_bilayer_qc_stays_quiet_for_a_water_row(tmp_path: Path, capsys):
    """And the same absence in a water run is not worth a word.

    Every kinase row would otherwise print a refusal on every streaming pass,
    which is how a real warning stops being read.
    """
    module = _ingest_module()
    run_dir = tmp_path / "abl1_apo"
    run_dir.mkdir()
    (run_dir / "trajectory.dcd").write_bytes(b"not really a dcd")

    assert module.membrane_chunks(run_dir, "abl1_apo") == []
    assert "REFUSING" not in capsys.readouterr().out


def test_a_missing_manifest_is_reported_rather_than_returned_empty(
        monkeypatch, tmp_path: Path, capsys):
    """The caller only logs its reference count when there were rows.

    So a manifest that had been renamed would drop every reference from the
    explorer and print nothing at all about having done so.
    """
    module = _ingest_module()
    monkeypatch.setitem(module.MANIFESTS, "abl1", tmp_path / "gone.jsonl")
    assert module.reference_rows("abl1") == {}
    assert "no manifest at" in capsys.readouterr().out


# --- resuming a run, which had been losing its barostat ---------------------
#
# A resumed run skips the branch that adds the barostat, and OpenMM restores
# `MonteCarloPressure` from the checkpoint whether or not a force consumes it.
# So a resumed context reported the barostat's own parameters while having no
# barostat in it, and production continued at constant volume with nothing in
# the log to say so.

def _run_md():
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_run_md_for_tests", root / "tools" / "run_md.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_two_run_types_get_the_barostat_each_needs():
    """One helper, because the two call sites had drifted apart.

    The membrane rows need the membrane barostat at zero surface tension so the
    bilayer sets its own area; a plain isotropic one would scale z with xy and
    squeeze the bilayer.
    """
    import openmm as mm

    module = _run_md()
    water = module.barostat(False, 300.0)
    membrane = module.barostat(True, 303.15)
    assert isinstance(water, mm.MonteCarloBarostat)
    assert isinstance(membrane, mm.MonteCarloMembraneBarostat)
    assert not isinstance(water, mm.MonteCarloMembraneBarostat)
    # Section 5.1: the bilayer is not held at an area.
    assert membrane.getDefaultSurfaceTension().value_in_unit(
        mm.unit.bar * mm.unit.nanometer) == 0.0


def test_an_equilibration_checkpoint_is_not_mistaken_for_production(
        tmp_path: Path):
    """The two look identical on disk and are resumed differently.

    Resuming an equilibration checkpoint as production would start with the
    protein still restrained and count the equilibration time against the
    production target, writing a short row without a word.
    """
    module = _run_md()
    # Nothing recorded: every checkpoint written before the marker existed was
    # a production one, so that is the backward-compatible answer.
    assert module.checkpoint_phase(tmp_path) == "production"

    (tmp_path / "checkpoint_phase.json").write_text('{"phase": "equil", "rung": 3}')
    assert module.checkpoint_phase(tmp_path) == "equil"

    (tmp_path / "checkpoint_phase.json").write_text('{"phase": "production"}')
    assert module.checkpoint_phase(tmp_path) == "production"

    # A marker that cannot be read must not be an exception in the middle of a
    # resume; it falls back to the same default as an absent one.
    (tmp_path / "checkpoint_phase.json").write_text("{not json")
    assert module.checkpoint_phase(tmp_path) == "production"


def test_saving_a_checkpoint_records_the_phase_beside_it(tmp_path: Path):
    import json as _json

    module = _run_md()

    class FakeSimulation:
        def __init__(self):
            self.saved_to = None

        def saveCheckpoint(self, path):
            self.saved_to = path
            Path(path).write_bytes(b"checkpoint")

    simulation = FakeSimulation()
    module.save_checkpoint(simulation, tmp_path, "equil", rung=2)
    assert simulation.saved_to == str(tmp_path / "state.chk")
    written = _json.loads((tmp_path / "checkpoint_phase.json").read_text())
    assert written == {"phase": "equil", "rung": 2}

    module.save_checkpoint(simulation, tmp_path, "production")
    written = _json.loads((tmp_path / "checkpoint_phase.json").read_text())
    assert written == {"phase": "production"}


def test_every_queued_row_declares_the_ligand_its_start_model_contains():
    """A row with a bound ligand cannot be built without its chemistry.

    `abl1_imatinib` had no `ligand` field, and amber14 has no template for
    imatinib, so the row would have failed at build time after the three rows
    ahead of it had already spent their GPU hours. Checked against the start
    models rather than against a list, because the list is the thing that was
    wrong.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_queue_md_for_tests", root / "tools" / "queue_md.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)

    # Water, ions and lipids are not ligands: amber14 and lipid17 template them.
    solvent = {"HOH", "WAT", "NA", "CL", "K", "MG", "ZN", "SO4", "PO4",
               "GOL", "EDO", "POPC", "PC", "OL", "PA"}
    # Nor are the modified amino acids PDBFixer rewrites to a standard one.
    # Taken from PDBFixer's own table rather than listed here, so this stays
    # true if that table changes. ABL1's start model carries PTR 393, which is
    # why it matters: it is phosphotyrosine and it becomes plain tyrosine.
    from pdbfixer.pdbfixer import substitutions

    replaced = set(substitutions)
    checked = 0
    for row in queue.QUEUE:
        start = root / "data" / "systems" / row["start"]
        if not start.is_file():
            continue            # the start models are gitignored
        checked += 1
        found = {line[17:20].strip() for line in start.read_text().splitlines()
                 if line.startswith("HETATM")} - solvent - replaced
        declared = {row["ligand"]} if row.get("ligand") else set()
        assert found <= declared, (
            f"{row['name']}: its start model contains {sorted(found - declared)} "
            f"but the row declares {sorted(declared) or 'no ligand'}")
    if checked:
        # One row at least must actually have a ligand, or this test would pass
        # by finding nothing anywhere.
        assert any(r.get("ligand") for r in queue.QUEUE)


def test_the_queue_hands_run_md_a_ccd_code_not_a_path(tmp_path: Path):
    """The two sides of `--ligand` have to mean the same thing.

    `run_md.py` takes a CCD residue code, because the chemistry comes from the
    dictionary and the coordinates from the deposition. I changed it to that and
    left the queue passing a path to the SDF, which nothing would have caught
    until the row ran: the pre-built system is reused, so the ligand argument is
    only read when a system has to be built.
    """
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_queue_for_ligand_test", root / "tools" / "queue_md.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)

    # Read off the row rather than by calling `run_row`, which would start a
    # subprocess. What matters is the value the queue holds.
    row = next(r for r in queue.QUEUE
               if r.get("ligand") and not r["membrane"])
    assert row["ligand"] == row["ligand"].upper()
    assert len(row["ligand"]) <= 3, "a CCD code is at most three characters"
    assert not row["ligand"].endswith(".sdf")

    # And the runner's own option must still be a plain string, not a Path, or
    # a code would be turned into a filename.
    source = (root / "tools" / "run_md.py").read_text()
    option = source[source.index('parser.add_argument("--ligand"'):]
    option = option[:option.index(")\n")]
    assert "type=Path" not in option, "--ligand must take a code, not a path"
    assert "CCD" in option, "the option should say it wants a CCD code"


def test_the_bilayer_qc_streams_chunks_without_losing_frames(tmp_path: Path):
    """The pass is streamed, so the chunk windows must still tile the run.

    It used to load the whole solvated trajectory, which for a finished 50 ns
    membrane row is 500 frames of 194,032 atoms: about 1.2 GB held at once for
    data used one 10 ns window at a time. `iterload` caps that at one chunk, and
    what has to be checked is the bookkeeping that replaced `loaded.n_frames`:
    contiguous windows, every frame counted once.

    The fixture is built here rather than committed, and is not a real
    structure: two planes of phosphorus in POPC residues and some water, which
    is all `phosphorus_indices` and `analyse_chunk` look for.
    """
    import importlib.util

    import mdtraj as md
    import numpy as np

    from flexappeal.refstates import membrane

    n_per_leaflet, n_water, n_frames = 8, 6, 10
    topology = md.Topology()
    chain = topology.add_chain()
    for _ in range(2 * n_per_leaflet):
        residue = topology.add_residue("POPC", chain)
        topology.add_atom("P", md.element.phosphorus, residue)
    water_chain = topology.add_chain()
    for _ in range(n_water):
        residue = topology.add_residue("HOH", water_chain)
        topology.add_atom("O", md.element.oxygen, residue)

    rng = np.random.default_rng(11)
    box_nm = 4.0
    xyz = np.zeros((n_frames, topology.n_atoms, 3))
    for frame in range(n_frames):
        flat = rng.uniform(0, box_nm, size=(topology.n_atoms, 2))
        xyz[frame, :, :2] = flat
        # Two headgroup planes about 3.8 nm apart, and water above the top one.
        xyz[frame, :n_per_leaflet, 2] = 2.9 + rng.normal(0, 0.02, n_per_leaflet)
        xyz[frame, n_per_leaflet:2 * n_per_leaflet, 2] = (
            1.0 + rng.normal(0, 0.02, n_per_leaflet))
        xyz[frame, 2 * n_per_leaflet:, 2] = rng.uniform(3.4, 3.9, n_water)

    trajectory = md.Trajectory(xyz, topology)
    trajectory.unitcell_lengths = np.tile([box_nm, box_nm, box_nm],
                                          (n_frames, 1))
    trajectory.unitcell_angles = np.tile([90.0, 90.0, 90.0], (n_frames, 1))

    run_dir = tmp_path / "b2ar_apo"
    run_dir.mkdir()
    trajectory[0].save_pdb(str(run_dir / "solvated.pdb"))
    trajectory.save_dcd(str(run_dir / "trajectory.dcd"))

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_ingest_stream_test", root / "tools" / "ingest.py")
    ingest = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ingest)

    # 0.2 ns windows over 1 ns of 100 ps frames: five chunks, so the
    # cross-chunk bookkeeping is actually exercised.
    original = membrane.CHUNK_NS
    membrane.CHUNK_NS = 0.2
    try:
        rows = ingest.membrane_chunks(run_dir, "b2ar_apo",
                                      expect_bilayer=True)
    finally:
        membrane.CHUNK_NS = original

    assert len(rows) > 1, "the chunking was not exercised"
    assert sum(r["n_frames"] for r in rows) == n_frames, "frames lost or doubled"
    end = 0.0
    for row in rows:
        assert row["start_ns"] == pytest.approx(end), "windows are not contiguous"
        end = row["end_ns"]
        assert row["n_lipids"] == 2 * n_per_leaflet
        assert row["leaflet_counts"] == [n_per_leaflet, n_per_leaflet]
    assert end == pytest.approx(n_frames * 0.1)


# --- G-REUSE: the module has to work for a protein it has never seen ---------

def _analyse_run():
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_analyse_run_for_tests", root / "tools" / "analyse_run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_annotate_passes_the_protein_to_the_catalogue_lookup():
    """`--protein` was accepted by the parser and then dropped.

    So CDK2 was annotated against ABL1's reference ensembles and reported a
    clean `intermediate` with no conflicts, which is what it would have looked
    like had it worked. Checked on the source because the alternative needs
    KLIFS, and what went wrong was an argument not being forwarded.
    """
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "refstates" / "cli.py").read_text()
    block = source[source.index('if family == "kinase":'):]
    block = block[:block.index("motif_residues, core, built")] \
        if "motif_residues, core, built" in block else block[:600]
    assert "protein=args.protein" in block, "the kinase setup is not given the protein"
    assert block.count("protein=args.protein") == 2, (
        "both families must be given the protein, not just one")
    # One shared default cannot be right for both: ABL1 is not a receptor.
    assert '"ABL1"' in block and '"adrb2_human"' in block

    # And the setups must accept it.
    module = _analyse_run()
    import inspect

    for function in (module.kinase_setup, module.gpcr_setup):
        parameters = inspect.signature(function).parameters
        assert "protein" in parameters, f"{function.__name__} takes no protein"
        assert "max_references" in parameters


def test_reference_cap_takes_the_best_and_leaves_the_examples_alone():
    """A state score is the distance to the nearest reference.

    CDK2 yields 337 active structures, so the three hundredth costs a download
    and changes almost nothing. The cap applies to the reuse path only, which is
    why no limit has to mean every entry: the committed results for both worked
    examples are computed over all of theirs.
    """
    module = _analyse_run()
    rows = [{"pdb": "AAAA", "quality_score": 4.0, "resolution": 2.5},
            {"pdb": "BBBB", "quality_score": 8.0, "resolution": 1.5},
            {"pdb": "CCCC", "quality_score": 6.0, "resolution": 3.0}]

    # No limit is every row, in catalogue order.
    assert module._best_entries(rows, None, key="quality_score") == rows
    assert module._best_entries(rows, 0, key="quality_score") == rows
    # A limit at or above the count changes nothing either.
    assert module._best_entries(rows, 5, key="quality_score") == rows

    # Quality: highest first, because a better structure scores higher.
    best = module._best_entries(rows, 2, key="quality_score")
    assert [r["pdb"] for r in best] == ["BBBB", "CCCC"]

    # Resolution: lowest first, because a better structure resolves finer.
    best = module._best_entries(rows, 2, key="resolution", ascending=True)
    assert [r["pdb"] for r in best] == ["BBBB", "AAAA"]


def test_reference_cap_puts_an_unscored_entry_last_either_way():
    """A row with no quality or resolution must not win the sort.

    It sorted last for a descending sort already; ascending needed its own
    fallback, because a missing value read as -inf would be the best resolution
    there is.
    """
    module = _analyse_run()
    rows = [{"pdb": "GOOD", "quality_score": 7.0, "resolution": 1.8},
            {"pdb": "BLANK"},
            {"pdb": "ALSO", "quality_score": 5.0, "resolution": 2.2}]

    kept = module._best_entries(rows, 2, key="quality_score")
    assert "BLANK" not in [r["pdb"] for r in kept]
    kept = module._best_entries(rows, 2, key="resolution", ascending=True)
    assert "BLANK" not in [r["pdb"] for r in kept]
    assert [r["pdb"] for r in kept] == ["GOOD", "ALSO"]


def test_a_thin_reference_ensemble_is_reported_in_the_output():
    """How many structures a score rests on is part of the score.

    CDK2's reuse run came back with one usable active reference after the core
    alignment, and without this the output looked the same as one computed
    against forty-five.
    """
    from flexappeal.refstates import cli

    assert cli.THIN_ENSEMBLE >= 2, "a threshold of one would never fire"
    source = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "refstates" / "cli.py").read_text()
    assert "THIN_ENSEMBLE" in source
    block = source[source.index("for state, ensemble in built.items():"):]
    block = block[:800]
    assert 'payload.setdefault("notes"' in block, (
        "the warning must reach the payload, not only the log")


def _queue_md():
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "tools" / "queue_md.py"
    spec = importlib.util.spec_from_file_location("_queue_md", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _written(tmp_path: Path, **over):
    """A `run.json` as `run_md.py` writes one, with fields overridden."""
    payload = {"name": "b2ar_apo", "production_ns": 50.0,
               "produced_ns": 50.0, "status": "complete"}
    payload.update(over)
    for key in [k for k, v in payload.items() if v is None]:
        del payload[key]
    run_dir = tmp_path / "b2ar_apo"
    run_dir.mkdir(exist_ok=True)
    (run_dir / "run.json").write_text(json.dumps(payload))
    return tmp_path


def test_a_short_row_is_not_certified_by_the_target_it_aimed_at(tmp_path: Path):
    """`complete` compared the target against itself, so it could not fail.

    `run.json` recorded `production_ns` as the target the row was asked for, and
    `complete` then checked that number against the same target. A row that
    stopped early still wrote 50 and was certified by it. Section 4ag: a resume
    that misread its clock would have ended 10 ns short with every log line
    looking ordinary.
    """
    queue = _queue_md()
    assert queue.complete(_written(tmp_path, produced_ns=50.0),
                          "b2ar_apo", 50.0), "a row at its target is complete"
    assert not queue.complete(_written(tmp_path, produced_ns=40.0),
                              "b2ar_apo", 50.0), (
        "40 ns of a 50 ns target is short, whatever production_ns says")


def test_a_short_row_is_judged_on_what_it_produced_not_what_it_claimed(
        tmp_path: Path):
    """The case the old check got wrong: target intact, row short."""
    queue = _queue_md()
    assert not queue.complete(
        _written(tmp_path, production_ns=50.0, produced_ns=39.9,
                 status="short:39.90_of_50ns"), "b2ar_apo", 50.0)


def test_rounding_at_the_last_step_does_not_make_a_row_short(tmp_path: Path):
    """An integer step count lands sub-picosecond off, which is not a shortfall."""
    queue = _queue_md()
    assert queue.complete(_written(tmp_path, produced_ns=49.9999),
                          "b2ar_apo", 50.0)


def test_a_row_recorded_before_produced_ns_existed_is_taken_at_its_word(
        tmp_path: Path):
    """Rows finished before the field existed have only the target.

    Same reasoning as `checkpoint_phase`: the fallback has to describe what is
    already on disk, not what the current code would write.
    """
    queue = _queue_md()
    assert queue.complete(_written(tmp_path, produced_ns=None),
                          "b2ar_apo", 50.0)
    assert not queue.complete(
        _written(tmp_path, produced_ns=None, production_ns=None),
        "b2ar_apo", 50.0), "neither field present is not a finished row"


def test_the_short_row_tolerance_agrees_between_the_runner_and_the_queue():
    """Two literals that have to match, held together by this test.

    `queue_md` keeps its own copy rather than importing `run_md`, which builds
    an OpenMM system at import on some paths. The duplication is deliberate and
    this is what stops it drifting.
    """
    import ast

    root = Path(__file__).resolve().parents[1]
    values = {}
    for name in ("run_md", "queue_md"):
        tree = ast.parse((root / "tools" / f"{name}.py").read_text())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Assign)
                    and any(getattr(t, "id", None) == "SHORT_ROW_TOLERANCE_NS"
                            for t in node.targets)):
                values[name] = ast.literal_eval(node.value)
    assert set(values) == {"run_md", "queue_md"}, values
    assert values["run_md"] == values["queue_md"], values


def test_a_row_already_at_its_target_still_records_itself(tmp_path: Path):
    """The branch that returned 0 and wrote nothing.

    A row resumed to exactly its target took the "nothing to do" exit, which
    announced itself and wrote no `run.json`. `queue_md.complete` reads an
    absent `run.json` as unfinished and so does `tools/ingest.py`, so a finished
    row would have shown as in flight for good.
    """
    import types

    module = _run_md()
    args = types.SimpleNamespace(name="b2ar_apo", membrane=True, equil_ns=None)
    module.write_run_json(tmp_path, args, target_ns=50.0, produced_ns=50.0,
                          platform_name="OpenCL", atoms=61234,
                          temperature=303.15)
    written = json.loads((tmp_path / "run.json").read_text())
    assert written["status"] == "complete"
    assert written["produced_ns"] == 50.0
    assert "ns_per_day" not in written and "wall_hours" not in written, (
        "section 0.3: omit a key you cannot produce rather than writing a "
        "stand-in, and a row that did no stepping has no speed")
    queue = _queue_md()
    assert queue.complete(tmp_path.parent, tmp_path.name, 50.0)


def test_a_short_row_records_itself_as_short_rather_than_complete(
        tmp_path: Path):
    """The status has to say so, because that is what the explorer reads."""
    import types

    module = _run_md()
    args = types.SimpleNamespace(name="b2ar_apo", membrane=True, equil_ns=None)
    module.write_run_json(tmp_path, args, target_ns=50.0, produced_ns=40.0,
                          platform_name="OpenCL", atoms=61234,
                          temperature=303.15, ns_per_day=59.7,
                          wall_hours=16.1)
    written = json.loads((tmp_path / "run.json").read_text())
    assert written["status"].startswith("short:"), written
    assert "40" in written["status"] and "50" in written["status"], (
        "the status should name both numbers, as failed:<reason> does")
    assert written["ns_per_day"] == 59.7, "a short row still measured a speed"


def test_a_resumed_row_does_not_erase_its_own_log(tmp_path: Path, monkeypatch):
    """`row_log.open("w")` truncated the row's history on every restart.

    b2ar_apo's five hours of staged equilibration went when its resume opened
    the same path, and survived only because the queue echoes row lines into its
    own log. Nothing reads these files, so appending after a separator costs
    nothing and keeps a row's account of itself in one place. Section 4ag.
    """
    import subprocess as sp

    queue = _queue_md()
    out_dir = tmp_path
    (out_dir / "abl1_apo.log").write_text(
        "[14:47:00] minimising\n[14:52:00] equilibrating 1.0 ns (NPT)\n")

    class FakeProcess:
        def __init__(self, *a, **k):
            self.stdout = iter(["[17:35:00] resumed from checkpoint at 4.000 ns\n"])

        def wait(self):
            return 0

    monkeypatch.setattr(sp, "Popen", FakeProcess)
    monkeypatch.setattr(queue, "MANIFEST", tmp_path / "md.jsonl")
    queue.run_row({"name": "abl1_apo", "start": "abl1_apo_start.pdb",
                   "membrane": False, "ligand": None, "note": "DFG-in, apo"},
                  out_dir=out_dir, production_ns=50.0, equil_ns=None,
                  smoke=False)

    written = (out_dir / "abl1_apo.log").read_text()
    assert "minimising" in written, (
        "the resume must not take the row's earlier log with it")
    assert "resumed from checkpoint at 4.000 ns" in written
    assert "started again" in written, (
        "a reader needs to see where one attempt ends and the next begins")


def test_the_production_clock_is_read_from_the_row_not_assumed(tmp_path: Path):
    """A row runs the code it started with, so only the row can be believed.

    The clock reset at production start and the phase marker were both added
    while rows were already in flight. A resume that mistook production time for
    production plus equilibration would have ended b2ar_apo 10 ns short with
    every log line looking ordinary, and what settled it was the row's own
    `state_data.csv`: step 4,300,000 at 4 fs is 17,200 ps exactly. Section 4ag.
    """
    queue = _queue_md()
    run_dir = tmp_path / "b2ar_apo"
    run_dir.mkdir()
    assert queue.production_ns_on_disk(run_dir) is None, (
        "no state file is not a reading of zero")
    (run_dir / "state_data.csv").write_text(
        '#"Step","Time (ps)","Temperature (K)"\n'
        '25000,100.0,303.1\n'
        '4300000,17200.0,304.4\n')
    assert queue.production_ns_on_disk(run_dir) == 17.2
    # 4,300,000 steps at 4 fs is 17,200 ps, so the clock had been reset: had it
    # carried the 10 ns equilibration the two would disagree by exactly that.
    assert 4_300_000 * 4e-6 == 17.2


def test_a_state_file_with_a_header_and_no_rows_reads_as_unknown(tmp_path: Path):
    """Section 0.3: no reading is None, never a stand-in zero."""
    queue = _queue_md()
    run_dir = tmp_path / "abl1_apo"
    run_dir.mkdir()
    (run_dir / "state_data.csv").write_text('#"Step","Time (ps)"\n')
    assert queue.production_ns_on_disk(run_dir) is None


def test_preflight_blocks_a_row_it_can_neither_reuse_nor_build(
        tmp_path: Path, monkeypatch, capsys):
    """The point of the check is to fail before the GPU hours, not after."""
    queue = _queue_md()
    monkeypatch.setattr(queue, "QUEUE", [
        {"name": "nothing_here", "start": "no_such_system.pdb",
         "membrane": False, "ligand": None, "note": "a row with no system"}])
    assert queue.preflight(tmp_path, 50.0) == 1, (
        "a row that cannot run should give a non-zero exit, so this can gate "
        "a queue"
    )
    assert "BLOCKED" in capsys.readouterr().out


def test_preflight_reads_and_writes_nothing(tmp_path: Path, monkeypatch):
    """Safe to run beside a production row, which is the only time it is wanted.

    Nothing here may build a system, step an integrator or open an OpenMM
    context: there is one GPU and section 2.2's `concurrent_md: 1` exists
    because two contexts on it helped neither row.
    """
    queue = _queue_md()
    run_dir = tmp_path / "abl1_apo"
    run_dir.mkdir()
    (run_dir / "system.xml").write_text("<System/>")
    (run_dir / "solvated.pdb").write_text("END\n")
    monkeypatch.setattr(queue, "QUEUE", [
        {"name": "abl1_apo", "start": "abl1_apo_start.pdb",
         "membrane": False, "ligand": None, "note": "DFG-in, apo"}])
    before = sorted(p.name for p in run_dir.iterdir())
    assert queue.preflight(tmp_path, 50.0) == 0
    assert sorted(p.name for p in run_dir.iterdir()) == before, (
        "preflight must not leave anything behind in a row's directory")


def _stream_producing(tmp_path: Path, *, queue_log: str, pgrep_script: str,
                      grace_checks: int = 2) -> int:
    """Run `stream_analysis.sh`'s `producing` with a stubbed `pgrep`.

    The script is sourced rather than executed, which its own guard allows, so
    the half-hour loop does not start inside the suite.
    """
    import subprocess

    root = Path(__file__).resolve().parents[1]
    log = tmp_path / "production.log"
    log.write_text(queue_log)
    stub_dir = tmp_path / "bin"
    stub_dir.mkdir()
    stub = stub_dir / "pgrep"
    stub.write_text(pgrep_script)
    stub.chmod(0o755)
    env = {
        "PATH": f"{stub_dir}:/usr/bin:/bin",
        "QUEUE_LOG": str(log),
        "GRACE_CHECKS": str(grace_checks),
        "GRACE_SLEEP": "0",
        "COUNTER": str(tmp_path / "calls"),
    }
    return subprocess.run(
        ["bash", "-c",
         f'source "{root}/tools/stream_analysis.sh"; producing'],
        cwd=root, env=env, capture_output=True).returncode


def test_a_gap_between_rows_does_not_end_the_streaming(tmp_path: Path):
    """The queue starts each row as a fresh subprocess, so there is a gap.

    At every row boundary there is a moment with no runner alive and nothing yet
    saying the queue is done. Looking once in that window would have ended the
    streaming for every row after the first and fired the final pass, the one
    that computes the bilayer QC and G6, on a partial trajectory. The same
    window opened for about a minute when the queue was restarted by hand at
    04:19 on 11 October, and the loop survived only because it was asleep.
    Section 4ag.
    """
    # Nothing on the first look, a queue process on the second.
    stub = ('#!/bin/bash\n'
            'n=$(cat "$COUNTER" 2>/dev/null || echo 0)\n'
            'echo $((n + 1)) > "$COUNTER"\n'
            '[ "$n" -ge 1 ] && echo 4242\n'
            'exit 0\n')
    assert _stream_producing(tmp_path, queue_log="[21:22] producing\n",
                             pgrep_script=stub) == 0, (
        "a transient gap is a row boundary, not the end of production")


def test_production_really_stopping_still_ends_the_streaming(tmp_path: Path):
    """The grace must not become a loop that cannot stop.

    If an absent process never counted, the loop would stream for ever and the
    final pass would never run, which is the fault the liveness check was
    rewritten for in the first place.
    """
    stub = '#!/bin/bash\nexit 1\n'
    assert _stream_producing(tmp_path, queue_log="[21:22] producing\n",
                             pgrep_script=stub) == 1


def test_the_queue_marker_ends_the_streaming_without_waiting(tmp_path: Path):
    """The queue saying so is definitive and needs no grace at all."""
    # A stub that would report production, so only the marker can be what
    # ends it.
    stub = '#!/bin/bash\necho 4242\nexit 0\n'
    assert _stream_producing(tmp_path, grace_checks=12,
                             queue_log="[21:22] producing\n[23:00] queue done\n",
                             pgrep_script=stub) == 1


def test_a_live_production_process_keeps_the_streaming_going(tmp_path: Path):
    stub = '#!/bin/bash\necho 4242\nexit 0\n'
    assert _stream_producing(tmp_path, queue_log="[21:22] producing\n",
                             pgrep_script=stub) == 0


# Classes used in a template and deliberately not defined in brand.css.
#
# `md-reanalyse` is a JavaScript hook: `analysis.js` selects on it and nothing
# styles it, which is a legitimate use of a class. The other four are dangling
# the way `md-stats` was, but they belong to pages outside the reference-state
# work (the analysis, prepare, run-detail and example-intro templates), so they
# are recorded here rather than fixed in passing. Anything new that dangles
# fails this test.
TEMPLATE_CLASSES_WITHOUT_CSS = {
    "md-reanalyse",
    "md-example-terminal", "md-source-form", "md-tile-note",
}


def test_every_class_a_template_uses_is_defined_in_the_stylesheet():
    """`md-tag` and five others were used in templates and defined nowhere.

    brand.css is the only stylesheet the pages load, and none of these
    templates carries a `<style>` block, so an undefined class is simply
    unstyled markup. It was not cosmetic: `md-stats` is the figures block at the
    top of every index card and every run page, and as an undefined `<dl>` it
    rendered as a default indented list rather than a row of labelled numbers.
    The Mol* container escaped only because its height is an inline style.

    Nothing was going to notice. The pages render, the page tests pass, and the
    difference is only visible against the design it was meant to have.
    """
    import re

    root = Path(__file__).resolve().parents[1]
    css = (root / "flexappeal" / "static" / "brand.css").read_text()
    defined = set(re.findall(r"\.(md-[a-z0-9-]+)", css))
    assert "md-stats" in defined and "md-tag" in defined, (
        "the two this test was written for")

    undefined: dict[str, set[str]] = {}
    for template in sorted((root / "flexappeal" / "templates").rglob("*.html")):
        body = template.read_text()
        for attribute in re.finditer(r'class="([^"]*)"', body, re.S):
            for token in re.findall(r"\bmd-[a-z0-9-]+\b", attribute.group(1)):
                if (token not in defined
                        and token not in TEMPLATE_CLASSES_WITHOUT_CSS):
                    undefined.setdefault(token, set()).add(template.name)

    assert not undefined, (
        "used in a template and defined in no stylesheet: "
        + "; ".join(f"{k} ({', '.join(sorted(v))})"
                    for k, v in sorted(undefined.items())))


def test_the_allowed_undefined_classes_are_all_still_used():
    """An allow-list that outlives its entries starts hiding things.

    If one of these is deleted from the templates, the entry should go too, so
    that it cannot quietly cover a future class of the same name.
    """
    root = Path(__file__).resolve().parents[1]
    bodies = "".join(p.read_text() for p in
                     (root / "flexappeal" / "templates").rglob("*.html"))
    stale = sorted(c for c in TEMPLATE_CLASSES_WITHOUT_CSS
                   if c not in bodies)
    assert not stale, f"no longer used in any template: {stale}"


def test_requiring_a_platform_refuses_rather_than_falling_back():
    """The fallback ladder is right on a machine you own and wrong on a rented one.

    `run_md.platform` tries CUDA, then OpenCL, then CPU, and the measured cost
    of the last step is 29.5 ns/day against 206. On a pod billed by the hour
    that is a 1 µs row dropping to a fiftieth of the speed and going on
    billing, with one line in a log nobody is reading at the time. So
    `deploy/runpod/job.sh` passes `--require-platform CUDA`.

    CUDA is not registered on this machine, so asking for it here exercises the
    refusal for real rather than through a mock.
    """
    module = _run_md()
    with pytest.raises(SystemExit) as refused:
        module.platform(require="CUDA")
    assert "Refusing to fall back" in str(refused.value)

    with pytest.raises(SystemExit) as unknown:
        module.platform(require="Metal")
    assert "not a platform this runner knows" in str(unknown.value)


def test_the_unconstrained_ladder_is_unchanged():
    """The guard must not have cost the Mac its fallback."""
    module = _run_md()
    found, props = module.platform()
    assert found.getName() in {"CUDA", "OpenCL", "CPU"}


def test_the_queue_passes_the_required_platform_to_each_row(tmp_path: Path,
                                                            monkeypatch):
    """Decided by the caller, so one queue can be pinned and another not."""
    import subprocess as sp

    queue = _queue_md()
    captured: list[list[str]] = []

    class FakeProcess:
        def __init__(self, command, *a, **k):
            captured.append(command)
            self.stdout = iter([])

        def wait(self):
            return 0

    monkeypatch.setattr(sp, "Popen", FakeProcess)
    monkeypatch.setattr(queue, "MANIFEST", tmp_path / "md.jsonl")
    row = {"name": "b2ar_apo", "start": "b2ar_apo_start.pdb",
           "membrane": True, "ligand": None, "note": "active, apo"}
    queue.run_row(row, out_dir=tmp_path, production_ns=1000.0, equil_ns=None,
                  smoke=False, require_platform="CUDA")
    queue.run_row(row, out_dir=tmp_path, production_ns=50.0, equil_ns=None,
                  smoke=False)
    assert "--require-platform" in captured[0]
    assert captured[0][captured[0].index("--require-platform") + 1] == "CUDA"
    assert "--require-platform" not in captured[1], (
        "the Mac keeps its fallback ladder")


def test_the_pod_pushes_its_checkpoint_even_when_the_row_fails(tmp_path: Path):
    """`set -e` made the push unreachable on exactly the path that needed it.

    `deploy/runpod/job.sh` ran the queue, then `status=$?`, then pushed. Under
    the `set -euo pipefail` at the top of that file a failing command exits
    immediately, so the assignment was never reached and the push never ran: a
    row that crashed, or a pod reclaimed mid-run, lost the checkpoint it had
    just spent hours earning. A failed row's checkpoint is the most valuable
    thing on the pod's disk, because the next pod resumes from it.

    Both shapes are run here rather than asserted about, since the fault is in
    how bash behaves and not in how the line reads.
    """
    import subprocess

    def outcome(body: str) -> str:
        script = tmp_path / "shape.sh"
        script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + body)
        done = subprocess.run(["bash", str(script)], capture_output=True,
                              text=True)
        return done.stdout.strip()

    assert outcome("false\nstatus=$?\necho pushed\n") == "", (
        "the old shape: set -e exits at the failing command and the push, "
        "which is the next line, never happens")
    assert outcome("status=0\nfalse || status=$?\necho pushed\n") == "pushed", (
        "making the call part of a condition is what stops set -e firing")

    job = (Path(__file__).resolve().parents[1]
           / "deploy" / "runpod" / "job.sh").read_text()
    assert "|| status=$?" in job
    assert "--require-platform CUDA" in job
    # The push must not be inside a success branch.
    push = job[job.index("sync.sh push"):]
    assert "status" not in push.split("\n")[0]


def test_the_runpod_scripts_say_they_have_never_been_run():
    """Section 5.3: written during the build, not executed.

    The claim belongs in the files themselves, because that is where someone
    about to run one of them will be looking.
    """
    root = Path(__file__).resolve().parents[1] / "deploy" / "runpod"
    for name in ("job.sh", "sync.sh", "README.md"):
        body = (root / name).read_text().lower()
        assert "never" in body and ("executed" in body or "been run" in body), (
            f"{name} should say it has not been run")
    readme = (root / "README.md").read_text()
    for gate in ("G1", "G2"):
        assert gate in readme, f"{gate} stands in front of this and should be named"
    # The env example names variables and carries no values.
    for line in (root / ".env.example").read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, _, value = line.partition("=")
            assert not value or value.startswith(("https://", "/")), (
                f"{key} should not carry a value: {line}")


def test_a_metric_that_is_nan_in_every_frame_is_caught(tmp_path: Path):
    """Section 10: "any metric with all-NaN for a motif is a test failure".

    This is the one failure the round-trip verifier cannot see. An all-NaN
    metric agrees with its own analysis file perfectly, because the file says
    NaN too. What makes it worth catching is how it renders: an empty column
    and a flat line, which on the page cannot be told apart from a quantity
    that genuinely did not move, and "did not move" is the headline this build
    is looking for.
    """
    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    assert _cli().verify_live_metrics(database) == [], (
        "a clean ingest has a reading for every metric")

    # Blank one metric across every frame, which is what a missing atom or an
    # unresolved residue produces.
    with store.connect(database) as connection:
        connection.execute(
            "UPDATE motif_frames SET value = NULL "
            "WHERE motif_id LIKE '%xDFG' AND metric = 'd_act_a'")
    problems = _cli().verify_live_metrics(database)
    assert problems, "a metric that was never measured went unnoticed"
    assert any("xDFG" in p and "d_act_a" in p for p in problems), problems
    assert any("all 20 frames" in p for p in problems), problems


def test_some_nan_frames_are_not_a_failure(tmp_path: Path):
    """A gap is honest; the whole thing missing is not.

    A salt bridge whose atoms are absent from a handful of frames has a real
    reading in the rest, and reporting that would make the check noise rather
    than a gate.
    """
    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    with store.connect(database) as connection:
        connection.execute(
            "UPDATE motif_frames SET value = NULL "
            "WHERE motif_id LIKE '%xDFG' AND metric = 'd_act_a' AND frame < 5")
    assert _cli().verify_live_metrics(database) == []


def test_a_measurement_with_no_state_score_is_caught(tmp_path: Path):
    """Half a reading is its own fault, and a different one.

    The distance is there and the score derived from it is not, so the motif
    has a measurement and no reading of it. The page's score column would be
    empty while its distance column was full, which looks like two different
    bugs and is one.
    """
    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    with store.connect(database) as connection:
        connection.execute(
            "UPDATE motif_frames SET s_score = NULL WHERE motif_id LIKE '%hinge'")
    problems = _cli().verify_live_metrics(database)
    assert problems and any("no state score" in p for p in problems), problems


def test_a_motif_with_no_frames_at_all_is_caught(tmp_path: Path):
    """Grouping cannot see a motif that contributed no rows to group.

    A motif whose every frame failed to ingest leaves nothing behind for a
    `GROUP BY` to find, so it would pass a check built only on groups while
    being as absent as one that is all NaN.
    """
    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    with store.connect(database) as connection:
        connection.execute("DELETE FROM motif_frames WHERE motif_id LIKE '%hinge'")
    problems = _cli().verify_live_metrics(database)
    assert problems and any("none of them" in p for p in problems), problems


def test_the_state_threshold_has_one_definition():
    """It was written out three times and named once.

    `state.DECISIVE_S` is the boundary between "active" and "intermediate".
    `store.summarise` wrote `0.2` as a literal, `web` wrote it twice more in the
    pair verdict, and the run template shaded its bands at a fourth copy. Three
    chances for the summary label to disagree with the per-frame label it
    summarises, silently, with nothing on the page to say which one had moved.
    """
    from flexappeal.refstates import panels as panels_module
    from flexappeal.refstates import state as state_module

    assert store.DECISIVE_S == state_module.DECISIVE_S
    # `web.py` was the third reader of this constant and is gone with the tab.
    # The panels at the end of a run page are the fourth, and they draw the
    # bands and write the sentence beside them from the same name.
    assert panels_module.DECISIVE_S == state_module.DECISIVE_S

    root = Path(__file__).resolve().parents[1]
    for name in ("flexappeal/refstates/store.py",
                 "flexappeal/refstates/panels.py",
                 "tools/make_figures.py"):
        body = (root / name).read_text()
        assert "DECISIVE_S" in body, f"{name} should use the constant"
        for line in body.splitlines():
            if ("0.2" in line and "final_s" in line) or (
                    "0.2" in line and "final >" in line):
                raise AssertionError(f"{name} still writes the threshold out: {line}")

    # The FINDINGS figure shaded from zero while the label rule cut at +0.2, so
    # a motif at +0.1 was drawn inside the active band and called intermediate
    # on the page beside it. Section 10 wants the line labelled as well.
    figures = (root / "tools" / "make_figures.py").read_text()
    assert "axis.axhspan(DECISIVE_S, 1" in figures
    assert "axis.axhspan(-1, -DECISIVE_S" in figures
    assert 'f"s = {edge:+.2f}"' in figures, (
        "the threshold line should carry its value")


def _payload_held_at(value: float) -> dict:
    """A payload where every motif sits at `value` for the whole run.

    `whole_score` labels the residue-weighted mean of the motifs' final-tenth
    means, so holding all of them at one value makes that mean the value and
    the label a direct reading of the threshold.
    """
    payload = _payload_with_frames(20)
    for frame in payload["frames"]:
        for motif in frame["motifs"]:
            motif["s"] = value
    return payload


def test_the_whole_protein_label_flips_exactly_at_the_threshold(tmp_path: Path):
    """Probed either side rather than read off the source.

    Asserting that a file contains a constant proves the spelling. This proves
    the behaviour, which is the thing that would actually drift if a fourth
    copy of the number appeared somewhere.
    """
    from flexappeal.refstates import state as state_module

    edge = state_module.DECISIVE_S
    for value, expected in ((edge + 0.01, "active"),
                            (edge - 0.01, "intermediate"),
                            (-edge - 0.01, "inactive"),
                            (-edge + 0.01, "intermediate")):
        database = tmp_path / f"t{value:+.3f}.db"
        store.ingest_analysis(database, _payload_held_at(value),
                              system_id="abl1", run_id="abl1_apo",
                              family="kinase", environment="water")
        whole = store.whole_score(database, "abl1_apo")
        assert whole.get("label") == expected, (value, whole)
        assert abs(whole["final_s"] - value) < 1e-9, whole


def _jumping_trajectory(n_frames: int = 6, box=(12.5, 13.1, 11.4)):
    """A small molecule that crosses a periodic boundary half the time.

    Modelled on what beta2AR actually does: whole in every frame, with its
    centre of mass alternating between the bottom of the box and the top, a
    full box length apart.
    """
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    for name in ("N", "CA", "C", "O"):
        topology.add_atom(name, md.element.carbon, residue)

    base = np.array([[0.0, 0.1, 0.0], [0.3, 0.1, 0.0],
                     [0.3, 0.4, 0.0], [0.0, 0.4, 0.3]])
    frames = []
    for index in range(n_frames):
        xyz = base.copy()
        # Every other frame is written in the other image along y.
        if index % 2:
            xyz[:, 1] += box[1]
        frames.append(xyz)
    traj = md.Trajectory(np.array(frames), topology)
    traj.unitcell_lengths = np.tile(np.array(box, dtype=float), (n_frames, 1))
    traj.unitcell_angles = np.tile(np.array([90.0, 90.0, 90.0]), (n_frames, 1))
    return traj


def test_a_molecule_on_a_periodic_boundary_is_brought_back_to_one_image():
    """beta2AR sits on the boundary and is written into either image.

    Its centre of mass alternates between y = 0.5 A and y = 129 A in a
    131.2 A box while the protein stays whole in every frame. Nothing measured
    within a frame cares, which is why the analysis is unaffected. The viewer
    does: it fits the camera to the whole trajectory, so a receptor occupying
    54 A of a 187 A bounding box is drawn small and off centre and teleports
    when the slider is dragged.
    """
    import numpy as np

    traj = _jumping_trajectory()
    shifts = _analyse_run_module().periodic_shifts(traj)
    fixed = traj.xyz + shifts[:, None, :]

    before = traj.xyz.reshape(-1, 3).ptp(axis=0) if hasattr(
        np.ndarray, "ptp") else (traj.xyz.reshape(-1, 3).max(axis=0)
                                 - traj.xyz.reshape(-1, 3).min(axis=0))
    after = (fixed.reshape(-1, 3).max(axis=0)
             - fixed.reshape(-1, 3).min(axis=0))
    assert before[1] > 13.0, "the fixture should span more than one box in y"
    assert after[1] < 1.0, f"the images were not collapsed: {after}"

    # And nothing about any single frame moved, which is what makes this safe
    # to do to a trajectory the page draws.
    assert np.allclose(traj.xyz.max(axis=1) - traj.xyz.min(axis=1),
                       fixed.max(axis=1) - fixed.min(axis=1))


def test_a_trajectory_that_does_not_jump_is_left_alone():
    """The correction has to be a no-op on the common case."""
    import numpy as np

    traj = _jumping_trajectory(n_frames=4)
    traj.xyz[1::2, :, 1] -= traj.unitcell_lengths[0, 1]   # undo the jumps
    shifts = _analyse_run_module().periodic_shifts(traj)
    assert np.allclose(shifts, 0.0), shifts


def test_a_non_orthorhombic_box_is_not_guessed_at():
    """A whole-box translation is the wrong correction for a triclinic cell.

    Section 0.3's shape: where the right answer is not available, do nothing
    and leave the data as it is rather than apply something that looks like a
    fix.
    """
    import numpy as np

    traj = _jumping_trajectory()
    traj.unitcell_angles[:, 0] = 109.47
    assert np.allclose(_analyse_run_module().periodic_shifts(traj), 0.0)


def test_a_trajectory_with_no_box_is_not_guessed_at():
    import numpy as np

    traj = _jumping_trajectory()
    traj.unitcell_lengths = None
    traj.unitcell_angles = None
    assert np.allclose(_analyse_run_module().periodic_shifts(traj), 0.0)


def test_the_bilayer_is_imaged_around_the_receptor(tmp_path: Path):
    """Not shifted with it, which was the first thing tried and was wrong.

    A bilayer spans the box, so its centroid does not jump when the protein
    crosses a boundary: it was already in one piece. Applying the protein's
    whole-box shift to it moved half its frames a box length and took its y
    span from 145 A to 269, which is a worse picture than the one being fixed.

    The minimum-image convention about the receptor's centre is the right
    operation: every phosphorus goes in whichever periodic copy sits nearest
    the receptor, giving a contiguous patch around the protein however either
    of them happened to be wrapped.
    """
    import mdtraj as md
    import numpy as np

    module = _analyse_run_module()
    pdb, dcd = _tiny_bilayer(tmp_path)
    out = tmp_path / "analysis"
    out.mkdir()

    plain_out = tmp_path / "plain"
    plain_out.mkdir()
    plain = module.write_lipid_bundle(pdb, dcd, "b2ar_apo", plain_out, stride=1)
    loose = md.load(str(plain_out / plain["trajectory"]),
                    top=str(plain_out / plain["topology"]))
    box = np.asarray(loose.unitcell_lengths)

    # A receptor one whole box away in y from where the lipids were written.
    centres = loose.xyz.mean(axis=1) + np.array([0.0, box[0, 1], 0.0])
    imaged = module.write_lipid_bundle(pdb, dcd, "b2ar_apo", out, stride=1,
                                       centres=centres)
    moved = md.load(str(out / imaged["trajectory"]),
                    top=str(out / imaged["topology"]))

    # Every atom should now sit within half a box of the receptor's centre,
    # which is what "imaged around it" means and what makes the patch
    # contiguous on the page.
    offset = moved.xyz - centres[:, None, :]
    assert np.all(np.abs(offset) <= box[:, None, :] / 2 + 1e-6), (
        np.abs(offset).max(axis=(0, 1)), box[0])


def test_a_frame_count_mismatch_leaves_the_bilayer_where_it_was(tmp_path: Path,
                                                                capsys):
    """Rather than imaging against centres that do not line up frame for frame."""
    import numpy as np

    module = _analyse_run_module()
    pdb, dcd = _tiny_bilayer(tmp_path)
    out = tmp_path / "analysis"
    out.mkdir()
    module.write_lipid_bundle(pdb, dcd, "b2ar_apo", out, stride=1,
                              centres=np.zeros((99, 3)))
    assert "may be drawn away from it" in capsys.readouterr().out


def test_the_analysis_payload_does_not_carry_the_centres_array():
    """It is a numpy array for the lipid bundle's use, not a result.

    Left in the viewer dict it would reach `json.dumps`, which cannot
    serialise it, and the run that found this out would be one that had just
    finished.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "analyse_run.py").read_text()
    assert 'viewer.pop("centres", None)' in source, (
        "the caller should take it out of the payload")


def test_the_readme_lists_every_positive_control_including_the_undecided_one():
    """Section 9 defines six. The public table listed five.

    The missing two were P6, the bilayer windows, which passes, and P5, an
    inactive start staying inactive, which is pending because the run that
    would decide it is third in a queue that is still draining. A table of five
    passes with the undecided one left out reads as six, which is the kind of
    omission that is not a lie in any single row.

    Matched on what each control is about rather than on its figures, which
    move as rows land. A control that is dropped from the table fails this; a
    control whose number changes does not.
    """
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text()
    table = readme[readme.index("### Positive controls"):]
    table = table[:table.index("### ", 1)] if "### " in table[1:] else table

    for control, needle in (
            ("P1 DFG label", "DFG label against KLIFS"),
            ("P2 alphaC label", "label against KLIFS"),
            ("P3 TM3-TM6 separation", "TM3"),
            ("P4 start scores as its own state", "scores as its own state"),
            ("P5 inactive start stays inactive", "stays inactive"),
            ("P6 bilayer windows", "Bilayer inside its published windows")):
        assert needle in table, f"{control} is missing from the README table"

    assert "pending" in table, (
        "a control that has not been decided should say so rather than be "
        "left out")


def _build_references_module():
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "build_references", root / "tools" / "build_references.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_resolution_each_catalogue_already_returned_is_kept():
    """It was read to filter on and then dropped.

    KLIFS applies a 2.8 A cut and GPCRdb a 3.6 A one, so both filters parse
    `resolution` out of the response. Neither carried it any further, and
    section 7's schema has a `resolution REAL` column that the reference
    browser renders a heading for: all 136 rows showed a dash under it, filled
    from nothing while the number sat in a response already in hand.
    """
    module = _build_references_module()
    assert module._resolution({"resolution": "2.79"}) == {"resolution": 2.79}
    assert module._resolution({"resolution": 3.2}) == {"resolution": 3.2}


def test_a_resolution_that_cannot_be_read_is_omitted_not_nulled():
    """Section 0.3: omit a key you cannot produce, never a stand-in.

    An EM entry with no resolution, or a null, or a string that is not a
    number, must leave the column alone rather than write a zero into it. Zero
    angstrom would sort to the top of a column people read to judge quality.
    """
    module = _build_references_module()
    for entry in ({}, {"resolution": None}, {"resolution": ""},
                  {"resolution": "n/a"}, {"resolution": 0},
                  {"resolution": -1}):
        assert module._resolution(entry) == {}, entry


def test_a_reference_keeps_its_resolution_through_the_store(tmp_path: Path):
    """And a reference without one leaves the column null rather than zero."""
    database = tmp_path / "t.db"
    store.create(database)
    with store.connect(database) as connection:
        store.put_system(connection, "b2ar", name="beta2AR", family="gpcr",
                         environment="membrane")
        store.put_reference(connection, "b2ar:4LDE", "b2ar", pdb="4LDE",
                            chain="A", source="GPCRdb", source_state="active",
                            our_state="active", resolution=2.79, metrics={})
        store.put_reference(connection, "b2ar:3SN6", "b2ar", pdb="3SN6",
                            chain="R", source="GPCRdb", source_state="active",
                            our_state="active", metrics={})
        rows = {r["pdb"]: r["resolution"] for r in connection.execute(
            "SELECT pdb, resolution FROM reference_structures")}
    assert rows["4LDE"] == 2.79
    assert rows["3SN6"] is None, (
        "no reading is null, which the page draws as a dash, not a zero")


def test_the_ingest_omits_a_resolution_it_does_not_have():
    """An older manifest has no such key and must ingest unchanged.

    The manifests are rebuilt rather than migrated, but they are gitignored, so
    the first ingest after a fresh clone reads whatever is on that machine.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "ingest.py").read_text()
    assert 'if row.get("resolution") is not None else {}' in source


def _section_key(number: str, suffix: str) -> tuple:
    """Sort key for a FINDINGS heading like `4`, `4b`, `4z`, `4aa`.

    Shorter suffixes come first, so `4z` precedes `4aa` the way the sections
    were actually written, and an empty suffix precedes them all.
    """
    return (int(number), len(suffix), suffix)


def test_the_findings_sections_are_in_order():
    """Three sections landed in the wrong place tonight before this existed.

    New sections are inserted by anchoring on the heading that should follow
    them, and the anchor was wrong three times running: 4o through 4af were
    written after section 5 and had to be moved back each time. The document is
    long enough that nobody reads it end to end to notice, and a section in the
    wrong place is invisible from inside the section.
    """
    import re

    root = Path(__file__).resolve().parents[1]
    headings = re.findall(r"^## (\d+)([a-z]*)\.",
                          (root / "FINDINGS.md").read_text(), re.M)
    assert len(headings) > 40, f"only found {len(headings)} sections"

    keys = [_section_key(n, s) for n, s in headings]
    for before, after, (n, s) in zip(keys, keys[1:], headings[1:]):
        assert before < after, (
            f"section {n}{s} follows "
            f"{headings[keys.index(before)][0]}{headings[keys.index(before)][1]}")
    assert len(set(keys)) == len(keys), "a section number is used twice"


def test_the_decision_numbers_run_in_order_without_gaps():
    """A document that is only ever appended to still drifts.

    Two decisions sharing a number, or a number skipped, both make a reference
    from `FINDINGS.md` or a code comment point at the wrong thing, and both are
    invisible unless something counts.
    """
    import re

    root = Path(__file__).resolve().parents[1]
    numbers = [int(n) for n in re.findall(
        r"^## D0*(\d+):", (root / "DECISIONS.md").read_text(), re.M)]
    assert numbers, "no decisions found"
    assert numbers == sorted(numbers), f"out of order: {numbers}"
    assert len(set(numbers)) == len(numbers), (
        f"a number is used twice: "
        f"{[n for n in numbers if numbers.count(n) > 1]}")
    assert numbers == list(range(numbers[0], numbers[0] + len(numbers))), (
        f"a number is missing between {numbers[0]} and {numbers[-1]}")


def test_every_decision_referred_to_in_the_findings_exists():
    """A reference to D0NN that names nothing is worse than no reference."""
    import re

    root = Path(__file__).resolve().parents[1]
    body = (root / "DECISIONS.md").read_text()
    decisions = set(re.findall(r"^## (D0*\d+):", body, re.M))
    # One heading form. D001 to D011 were written `## D001 - Title` and
    # everything after `## D012: title`, which is invisible until something
    # tries to find them all: a reference check written against the colon form
    # reported D004, D009 and D010 as decisions that did not exist.
    assert not re.search(r"^## D\d+ - ", body, re.M), (
        "a decision heading still uses the old separator")
    cited = set(re.findall(r"\b(D0\d\d)\b",
                           (root / "FINDINGS.md").read_text()))
    missing = sorted(cited - decisions)
    assert not missing, f"FINDINGS cites decisions that do not exist: {missing}"


def test_the_motif_figure_draws_each_half_mean():
    """The figure without them misled its own author.

    beta2AR's TM5 oscillates about zero with a scatter of 0.12, and reading the
    trace I called it a drift towards inactive, twice, when its halves are
    -0.034 and +0.016. A trace that crosses a line the eye is drawn to looks
    like it is moving across it. The halves are the check that says otherwise,
    they were already in the summary, and the figure is where anyone who is
    going to be misled is looking.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "make_figures.py").read_text()
    body = source[source.index("for motif, points in sorted(series.items())"):]
    body = body[:body.index("axis.set_xlabel")]
    assert "midpoint = len(points) // 2" in body, (
        "the halves should split where `summarise` splits them")
    assert "mean = sum(p[1] for p in half) / len(half)" in body
    assert "thick segments" in source, "the caption should say what they are"


def test_the_motif_figure_counts_frames_and_not_rows():
    """`motif_frames` holds one row per metric, and a motif has two.

    Without DISTINCT the figure plots every frame twice and averages over the
    duplicates. The means survive that, since both rows of a frame carry the
    same score, which is exactly what makes it the kind of fault that waits:
    it would start biasing the halves the day the two metrics carried different
    scores. Section 4m is the same error in the reference ensembles.
    """
    source = (Path(__file__).resolve().parents[1]
              / "tools" / "make_figures.py").read_text()
    assert "SELECT DISTINCT motif_id, time_ns, s_score" in source


def test_the_half_means_agree_with_the_stored_summary(tmp_path: Path):
    """Two paths to the same two numbers, so the figure cannot drift from the table.

    The figure computes its halves from the frames it plots; `summarise`
    computes them in SQL. A reader who checks one against the other should find
    them equal, and on the live run they agree to three decimals.
    """
    payload = _payload_with_frames(20)
    database = tmp_path / "t.db"
    store.ingest_analysis(database, payload, system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    summary = store.summarise(database, "abl1_apo")

    with store.connect(database, read_only=True) as connection:
        rows = list(connection.execute(
            "SELECT DISTINCT motif_id, time_ns, s_score FROM motif_frames "
            "WHERE run_id = ? AND s_score IS NOT NULL ORDER BY time_ns",
            ("abl1_apo",)))
    series: dict[str, list[tuple[float, float]]] = {}
    for row in rows:
        series.setdefault(row["motif_id"], []).append(
            (row["time_ns"], row["s_score"]))

    checked = 0
    for motif, points in series.items():
        points.sort()
        midpoint = len(points) // 2
        if midpoint < 2:
            continue
        first = sum(p[1] for p in points[:midpoint]) / midpoint
        second = sum(p[1] for p in points[midpoint:]) / len(points[midpoint:])
        stored = summary["motifs"][motif]
        assert abs(first - stored["first_half_mean_s"]) < 1e-9, motif
        assert abs(second - stored["second_half_mean_s"]) < 1e-9, motif
        checked += 1
    assert checked >= 2, "both motifs should have been compared"


# The constants a page might be tempted to write out, and where they live.
# Extend this when a new threshold appears rather than when it is found in a
# sentence: the whole point is that the finding should be automatic.
PAGE_THRESHOLDS = {
    "state.DECISIVE_S": ("state", "DECISIVE_S"),
    "membrane.APL_EXPECTED_A2": ("membrane", "APL_EXPECTED_A2"),
    "membrane.APL_ACCEPTABLE_A2": ("membrane", "APL_ACCEPTABLE_A2"),
    "membrane.THICKNESS_EXPECTED_A": ("membrane", "THICKNESS_EXPECTED_A"),
    "gpcr.IONIC_LOCK_A": ("gpcr", "IONIC_LOCK_A"),
    "kinase.ALPHAC_SALT_BRIDGE_A": ("kinase", "ALPHAC_SALT_BRIDGE_A"),
    "pocket.R_SPINE_CONTACT_A": ("pocket", "R_SPINE_CONTACT_A"),
    "pocket.R_SPINE_CONTACT_CRYSTAL_A": ("pocket", "R_SPINE_CONTACT_CRYSTAL_A"),
    "sources.KINASE_MAX_RESOLUTION_A": ("sources", "KINASE_MAX_RESOLUTION_A"),
    "sources.GPCR_MAX_RESOLUTION_A": ("sources", "GPCR_MAX_RESOLUTION_A"),
    "kinase.DFG_PUBLISHED_A": ("kinase", "DFG_PUBLISHED_A"),
}


def _template_prose(path: Path) -> str:
    """A template with everything that is not prose taken out.

    Scripts, Jinja comments, Jinja expressions and statements, and tags with
    their attributes. What is left is what a reader sees, which is the only
    place a hard-coded threshold does damage that nothing can catch.
    """
    import re

    body = path.read_text()
    for pattern in (r"<script.*?</script>", r"\{#.*?#\}", r"\{\{.*?\}\}",
                    r"\{%.*?%\}"):
        body = re.sub(pattern, " ", body, flags=re.S)
    return re.sub(r"<[^>]+>", " ", body)


def test_no_page_writes_a_threshold_out_in_prose():
    """Four instances in three templates, and the pattern had already cost three.

    A threshold in prose cannot be executed, so it cannot be wrong in a way
    anything notices: the index told a reader the pass band was +/-0.2, the
    About page gave the R-spine contact as 7.7 A, the reference browser said
    the kinase set was filtered at 2.8 A. All four were right. All four would
    have stayed on the page unchanged if the constant behind them moved, and a
    reader checking a verdict against a sentence would have been checking it
    against a number nothing uses.

    Only forms carrying a decimal point are searched. A bare integer is too
    common in prose to tell apart from a residue offset: `alphaC-Glu+4` is not
    `IONIC_LOCK_A`, and a check that cannot see the difference would be turned
    off within a week.
    """
    import importlib
    import re

    root = Path(__file__).resolve().parents[1]
    values: dict[str, list[float]] = {}
    for label, (module, attribute) in PAGE_THRESHOLDS.items():
        loaded = importlib.import_module(f"flexappeal.refstates.{module}")
        value = getattr(loaded, attribute)
        values[label] = list(value) if isinstance(value, tuple) else [value]

    found = []
    for path in sorted((root / "flexappeal" / "templates").glob("refstates*.html")):
        prose = _template_prose(path)
        for label, numbers in values.items():
            for number in numbers:
                form = f"{number:.1f}"
                if "." not in form:
                    continue
                for match in re.finditer(
                        r"(?<![\w.])" + re.escape(form) + r"(?![\w.])", prose):
                    context = " ".join(
                        prose[max(0, match.start() - 40):match.end() + 20].split())
                    found.append(f"{path.name}: {label} as {form} ... {context}")

    assert not found, (
        "a page writes a threshold out instead of rendering the constant:\n  "
        + "\n  ".join(sorted(set(found))))


def test_the_threshold_scan_would_find_one_if_it_were_there():
    """A check that cannot fail is worth nothing, and this one strips a lot out.

    Written against a template of its own rather than by breaking a real one,
    so it proves the scan sees prose and ignores the four places a number may
    legitimately sit: a script, a Jinja expression, a Jinja comment and a tag
    attribute.
    """
    import re

    sample = (
        '<div class="md-card" data-threshold="0.2">\n'
        "{# a comment mentioning 0.2 #}\n"
        '<p>rendered {{ "%.1f"|format(decisive_s) }} here</p>\n'
        "<script>const DECISIVE = 0.2;</script>\n"
        "<p>and written out as 0.2 in a sentence</p>\n")
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "sample.html"
        path.write_text(sample)
        prose = _template_prose(path)
    hits = re.findall(r"(?<![\w.])0\.2(?![\w.])", prose)
    assert len(hits) == 1, (
        f"the scan should see the sentence and nothing else, saw {len(hits)}: "
        f"{prose!r}")


def _app_with(database: Path):
    """A Flask app pointed at one database, for the view-level checks."""
    from flexappeal.webapp import create_app

    app = create_app()
    app.config["REFSTATES_DB"] = str(database)
    return app


def test_the_committed_figures_are_not_older_than_the_code_that_draws_them():
    """They were, by nine hours and three changes to the drawing.

    `data/figures/` is committed, so it is what anyone browsing the repository
    sees, and the final pass that regenerates it does not run until the queue
    drains. Between the two, the repository can show a figure drawn by code
    that no longer exists: the committed set shaded its bands from zero rather
    than from `DECISIVE_S` and carried no half-means, both of which the
    document beside them describes.

    Compared on commit dates and not on mtimes. A checkout gives every file the
    time it was written to disk, so an mtime comparison here would mean
    whatever order git happened to unpack them in, which is a test that passes
    or fails for no reason. What git does keep is when each path was last
    committed, and "the drawing changed after the drawings did" is exactly the
    question.
    """
    import subprocess

    root = Path(__file__).resolve().parents[1]
    if not (root / ".git").exists():
        pytest.skip("not a git checkout")
    if not list((root / "data" / "figures").glob("*.png")):
        pytest.skip("no figures committed in this checkout")

    def last_commit(path: str) -> int:
        done = subprocess.run(
            ["git", "log", "-1", "--format=%ct", "--", path],
            cwd=root, capture_output=True, text=True)
        return int(done.stdout.strip() or 0)

    script = last_commit("tools/make_figures.py")
    figures = last_commit("data/figures")
    if not script or not figures:
        pytest.skip("no history for one of the two paths")
    assert figures >= script, (
        "tools/make_figures.py was committed after the figures it draws, so "
        "data/figures/ is showing output from code that no longer exists; "
        "run it and commit the result")


def _tiny_kinase_topology():
    """A topology with the residue numbers ABL1's motifs actually use."""
    import mdtraj as md

    topology = md.Topology()
    chain = topology.add_chain()
    for number in list(range(249, 256)) + list(range(316, 319)) + \
            list(range(380, 390)):
        residue = topology.add_residue("ALA", chain, resSeq=number)
        for name in ("N", "CA", "C"):
            topology.add_atom(name, md.element.carbon, residue)
    return topology


def test_the_cartoon_number_is_the_motif_token_number():
    """Section 6.2's cartoon colouring, carried in the B-factor column.

    The vendored Mol* bundle is minified and exports no MolScript, so a
    component per motif is out of reach from the page; `uncertainty` reads the
    B-factor, so the number rides in the structure. The column's meaning is
    exactly the motif's `--motif-N` number, which is what lets the page use the
    same palette it uses on every chart.
    """
    module = _analyse_run_module()
    topology = _tiny_kinase_topology()
    values = module.motif_bfactors(
        topology, {"g.l": list(range(249, 255)), "hinge": [316, 317, 318]},
        "kinase")
    assert values is not None
    by_residue = {atom.residue.resSeq: int(values[atom.index])
                  for atom in topology.atoms}
    assert by_residue[250] == 4, "g.l is --motif-4"
    assert by_residue[317] == 5, "hinge is --motif-5"
    assert by_residue[385] == 0, "nothing covers it"


def test_an_overlap_goes_to_the_smaller_motif():
    """xDFG is four residues inside a twenty-six residue activation loop.

    Colouring the overlap by the larger one would hide the motif the whole
    DFG-in/out question is about.
    """
    module = _analyse_run_module()
    topology = _tiny_kinase_topology()
    values = module.motif_bfactors(
        topology,
        {"xDFG": [380, 381, 382, 383], "A-loop": list(range(381, 390))},
        "kinase")
    by_residue = {atom.residue.resSeq: int(values[atom.index])
                  for atom in topology.atoms}
    for number in (380, 381, 382, 383):
        assert by_residue[number] == 2, f"{number} should be xDFG"
    for number in (384, 385, 389):
        assert by_residue[number] == 3, f"{number} should be the A-loop"


def test_nothing_to_colour_by_leaves_the_cartoon_alone():
    """A column of zeroes would paint the whole protein one colour."""
    module = _analyse_run_module()
    topology = _tiny_kinase_topology()
    assert module.motif_bfactors(topology, {}, "kinase") is None
    assert module.motif_bfactors(topology, {"g.l": [250]}, "nonsense") is None
    assert module.motif_bfactors(topology, None, "kinase") is None
    # A motif the token table does not name contributes nothing, and if that
    # leaves nothing at all the answer is None rather than a zeroed column.
    assert module.motif_bfactors(topology, {"not-a-motif": [250]},
                                 "kinase") is None


def test_the_motif_palette_is_reversed_for_the_uncertainty_theme():
    """Unreversed it looked entirely deliberate, which is why it needs a test.

    `uncertainty` runs from certain to uncertain, so it consumes the colour
    list backwards: the background took motif 6's blue and every motif took the
    colour of its opposite number. Nothing errored and the structure was still
    a structure.
    """
    script = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "static" / "analysis.js").read_text()
    assert "colours.slice().reverse()" in script
    assert "name: 'uncertainty'" in script
    assert "domain: [-0.5, MOTIF_SLOTS + 0.5]" in script, (
        "each integer should sit in the middle of its own bin")


def test_the_viewer_colours_come_from_the_stylesheet_as_numbers():
    """Mol* wants 0xRRGGBB. Handed "#0072b2" it drew a near-black cartoon.

    The worst available outcome: the theme applied, nothing errored, and the
    structure was there and unreadable.
    """
    script = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "static" / "analysis.js").read_text()
    assert "function cssColourToInt" in script
    assert "getPropertyValue('--motif-' + i)" in script, (
        "the palette should stay single-sourced in brand.css")
    assert "if (value === null) return null;" in script, (
        "an unparseable colour should leave the default theme alone")


def _overlay_fixture(tmp_path: Path, shift=(0.0, 0.0, 0.0)):
    """A one-frame run and a reference PDB over the same residue numbers.

    Public coordinates: a short synthetic chain, generated here and not
    committed, which is what the fixture rule asks for.
    """
    import mdtraj as md
    import numpy as np

    numbers = list(range(10, 22))
    topology = md.Topology()
    chain = topology.add_chain()
    for number in numbers:
        residue = topology.add_residue("ALA", chain, resSeq=number)
        topology.add_atom("N", md.element.nitrogen, residue)
        topology.add_atom("CA", md.element.carbon, residue)
        topology.add_atom("C", md.element.carbon, residue)

    # A gentle helix-ish path, so the superposition is not degenerate.
    coords = []
    for i, _ in enumerate(numbers):
        base = np.array([np.cos(i * 0.6) * 5.0, np.sin(i * 0.6) * 5.0, i * 1.5])
        for offset in ((-1.2, 0.0, 0.0), (0.0, 0.0, 0.0), (1.2, 0.0, 0.0)):
            coords.append(base + np.array(offset))
    xyz = np.array([coords], dtype=np.float32) / 10.0        # mdtraj is in nm
    frame = md.Trajectory(xyz, topology)

    reference = tmp_path / "ref.pdb"
    moved = md.Trajectory(xyz + np.array(shift, dtype=np.float32) / 10.0,
                          topology)
    moved.save_pdb(str(reference))
    return frame, reference, numbers


def _built_many(references: dict, core: list[int], chain: str = "A"):
    """An aligned ensemble holding a stored core for each named reference.

    The overlay recovers each reference's own-to-common transform by fitting
    its own core onto the core stored here, so every code it may choose needs
    an entry, and they all have to be in **one** frame. Giving each its own
    frame looks simpler and is incoherent: the overlay takes the first code
    alphabetically as the target, so the others would be reconstructed into a
    frame nothing else is in.

    So this does what `align_to_common_frame` does: fits every core onto the
    first code's, which is then the common frame by construction.
    """
    import types

    import numpy as np

    from flexappeal.refstates import ensembles as ens
    from flexappeal.refstates.traj import apply_transform, superpose

    own = {code: np.asarray(ens.motif_coordinates(
        ens.ca_by_residue(path, chain), core), dtype=float)
        for code, path in references.items()}
    target = own[min(own)]
    return {"active": types.SimpleNamespace(by_code={
        code: {ens.CORE_KEY: [tuple(p) for p in
                              apply_transform(coords,
                                              *superpose(coords, target))]}
        for code, coords in own.items()})}


def _built_from(reference: Path, core: list[int], chain: str = "A"):
    """A minimal aligned ensemble whose target core is this reference's own.

    The overlay reconstructs the analysis's frame in two hops, own to common
    and common to run. Making the common frame this reference's own frame makes
    both hops the identity, so the fixture tests the bookkeeping rather than
    the arithmetic, which `test_the_two_hop_composition_matches_the_long_way`
    covers separately.
    """
    import types

    from flexappeal.refstates import ensembles as ens

    coords = ens.motif_coordinates(ens.ca_by_residue(reference, chain), core)
    return {"active": types.SimpleNamespace(
        by_code={"TEST": {ens.CORE_KEY: coords}})}


def test_an_overlay_that_disagrees_with_the_analysis_is_not_written(tmp_path: Path):
    """D023 refused to build this because a ghost in the wrong place is worse
    than no ghost. The user overruled that, so the wrong place has to be loud.

    RMSD is invariant under a rigid transform applied to both sets, so
    superposing the reference onto the run must reproduce exactly the distance
    the analysis got superposing the run onto the reference. Here the stored
    distance is a lie, and the overlay must refuse rather than draw.
    """
    module = _analyse_run_module()
    frame, reference, numbers = _overlay_fixture(tmp_path)
    rows = {"active": [{"pdb": "TEST", "path": str(reference), "chain": "A"}]}
    core = numbers[:8]
    motifs = {"m1": numbers[8:]}
    out = tmp_path / "out"
    out.mkdir()

    written = module.write_overlay_bundle(
        frame, "run", out, rows_by_state=rows,
        motifs_for=lambda code: motifs, core_for=lambda code: core,
        core_residues=core, motif_residues=motifs,
        motif_records={"m1": {"nearest_active": "TEST", "d_act_a": 9.99}},
        built=_built_from(reference, core))
    assert written == {}, "a disagreeing overlay should not be written"
    assert not list(out.glob("*.pdb"))


def test_an_overlay_that_agrees_is_written(tmp_path: Path):
    """The same fixture with the distance the geometry actually gives."""
    module = _analyse_run_module()
    frame, reference, numbers = _overlay_fixture(tmp_path)
    rows = {"active": [{"pdb": "TEST", "path": str(reference), "chain": "A"}]}
    core = numbers[:8]
    motifs = {"m1": numbers[8:]}
    out = tmp_path / "out"
    out.mkdir()

    # Reference and run are the same coordinates, so after the fit the motif
    # sits exactly on itself.
    written = module.write_overlay_bundle(
        frame, "run", out, rows_by_state=rows,
        motifs_for=lambda code: motifs, core_for=lambda code: core,
        core_residues=core, motif_residues=motifs,
        motif_records={"m1": {"nearest_active": "TEST", "d_act_a": 0.0}},
        built=_built_from(reference, core))
    assert "active" in written, written
    assert written["active"]["code"] == "TEST"
    assert written["active"]["motifs_checked"] == 1
    assert (out / "run.ref_active.pdb").is_file()


def test_an_overlay_with_nothing_to_check_against_is_not_written(tmp_path: Path):
    """A reference nearest for no motif cannot be verified, so it is refused.

    This is what made the selection rule change: picking by lowest mean motif
    RMSD chose a reference that was nearest for nothing, and the guard then
    deleted the inactive overlay rather than protecting it.
    """
    module = _analyse_run_module()
    frame, reference, numbers = _overlay_fixture(tmp_path)
    rows = {"active": [{"pdb": "TEST", "path": str(reference), "chain": "A"}]}
    core = numbers[:8]
    motifs = {"m1": numbers[8:]}
    out = tmp_path / "out"
    out.mkdir()

    written = module.write_overlay_bundle(
        frame, "run", out, rows_by_state=rows,
        motifs_for=lambda code: motifs, core_for=lambda code: core,
        core_residues=core, motif_residues=motifs,
        motif_records={"m1": {"nearest_active": "SOMEONE-ELSE",
                              "d_act_a": 0.0}},
        built=_built_from(reference, core))
    assert written == {}


def test_the_overlay_prefers_a_reference_it_can_actually_check(
        tmp_path: Path):
    """Verifiable first, then best mean. Both halves were learned from failing.

    Minimising the mean alone chose a reference nearest for nothing, so there
    was nothing to verify and the guard deleted the overlay. Maximising the
    weight a reference is nearest for chose beta2AR's 9CHU, nearest for TM5 at
    1.6 A and 27 A away averaged over its motifs: a motif read in the wrong
    place never wins a minimum, so it never reaches the analysis, and an
    overlay draws the whole structure.
    """
    module = _analyse_run_module()
    frame, near, numbers = _overlay_fixture(tmp_path)
    # A second reference, displaced, so it is further away on every motif.
    far_dir = tmp_path / "far"
    far_dir.mkdir()
    _, far, _ = _overlay_fixture(far_dir, shift=(3.0, 0.0, 0.0))

    core = numbers[:8]
    motifs = {"m1": numbers[8:]}
    rows = {"active": [{"pdb": "FAR", "path": str(far), "chain": "A"},
                       {"pdb": "NEAR", "path": str(near), "chain": "A"}]}
    out = tmp_path / "out"
    out.mkdir()
    written = module.write_overlay_bundle(
        frame, "run", out, rows_by_state=rows,
        motifs_for=lambda code: motifs, core_for=lambda code: core,
        core_residues=core, motif_residues=motifs,
        motif_records={"m1": {"nearest_active": "NEAR", "d_act_a": 0.0}},
        built=_built_many({"NEAR": near, "FAR": far}, core))
    assert written["active"]["code"] == "NEAR", written


def test_the_ghost_is_coloured_after_loading_not_by_the_preset():
    """`representationParams` is ignored by this bundle.

    The ghosts came back in Mol*'s default chain colours at full opacity, which
    is two bright proteins on top of the run rather than one behind it. The
    component refs are differenced across the load to find the new structure's
    representations.
    """
    script = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "static" / "analysis.js").read_text()
    assert "function ghostComponents" in script
    assert "componentRefs(viewer.plugin)" in script
    assert "old.type.params.alpha = 0.3" in script
    assert "getPropertyValue('--state-' + state)" in script, (
        "the ghost should take the state colour the badges use")


def test_the_viewport_is_hidden_until_its_camera_has_settled():
    """One picture, not two.

    Mol* frames the camera on what it has loaded and animates into that
    framing on its own schedule, so the camera a run was packed with can only
    be applied after that has finished. With the viewport visible throughout,
    a reader watched the structure appear in one view, hold still, and then
    swing into another as the playback started.

    The stylesheet hides it and the viewer reveals it once the camera has
    stopped moving. A timer reveals it too, because a load that fails must not
    leave a blank rectangle on the page.
    """
    root = Path(__file__).resolve().parents[1] / "flexappeal" / "static"
    script = (root / "analysis.js").read_text()
    css = (root / "brand.css").read_text()

    assert "#molstar-viewer.is-settling { opacity: 0; }" in css
    assert "host.classList.add('is-settling')" in script
    assert "host.classList.remove('is-settling')" in script
    assert "setTimeout(reveal, 8000)" in script, "a failed load must still reveal"
    assert "orientToMembraneSoon(viewer, host, revealNow)" in script


def test_the_camera_is_pointed_before_the_bilayer_is_fetched():
    """The bilayer is 17 MB, and the camera cannot wait for it.

    The orientation used to be applied in the branch that runs when there is no
    bilayer to load, so on a membrane run it happened only once the bilayer had
    arrived: seconds of the receptor in Mol*'s default view, and then a swing.
    Both now happen from the same place, the camera first.
    """
    script = (Path(__file__).resolve().parents[1]
              / "flexappeal" / "static" / "analysis.js").read_text()
    body = script[script.index("orientToMembraneSoon(viewer, host, revealNow)"):]
    body = body[: body.index("const scrubbing")]
    assert "lipid-toggle" in body, "the bilayer is still dispatched after it"
    assert "} else {" not in body, (
        "the camera should not be in the branch for runs without a bilayer")


def _fake_pocket(resolved: dict, regions: dict):
    """A stand-in PocketMap: position to residue, and position to region."""
    import types

    return types.SimpleNamespace(by_index=dict(resolved),
                                 region_by_index=dict(regions))


def test_the_core_is_intersected_over_positions_not_residue_numbers():
    """D025 deferred this because residue numbers cannot be intersected.

    ABL1's references are split between the ABL1a and ABL1b numbering nineteen
    residues apart, so intersecting numbers across them collapses to three
    residues. KLIFS positions are numbering independent, which is the whole
    reason this package addresses residues through the catalogue, so the
    intersection is taken there and resolved back per structure.
    """
    from flexappeal.refstates.traj import intersect_core

    regions = {1: "alphaE", 2: "alphaE", 3: "alphaE", 4: "hinge"}
    # Same three scaffold positions, numbered nineteen apart, as ABL1a/ABL1b.
    a = _fake_pocket({1: 100, 2: 101, 3: 102, 4: 103}, regions)
    b = _fake_pocket({1: 119, 2: 120, 3: 121, 4: 122}, regions)
    shared = intersect_core({"A": a, "B": b})

    assert shared["positions"] == [1, 2, 3], "the hinge is mobile and excluded"
    assert shared["residues_by"]["A"] == [100, 101, 102]
    assert shared["residues_by"]["B"] == [119, 120, 121]
    # The point of it: no residue number is shared, and the cores correspond.
    assert not set(shared["residues_by"]["A"]) & set(shared["residues_by"]["B"])


def test_a_structure_missing_a_position_narrows_everyone():
    """Which is what makes every core the same length by construction.

    CDK2's references resolve 42 positions where one resolves 41, and the
    alignment used to drop whichever did not match the modal length.
    """
    from flexappeal.refstates.traj import intersect_core

    regions = {i: "alphaE" for i in range(1, 5)}
    full = _fake_pocket({1: 10, 2: 11, 3: 12, 4: 13}, regions)
    short = _fake_pocket({1: 20, 2: 21, 4: 23}, regions)
    shared = intersect_core({"FULL": full, "SHORT": short})

    assert shared["positions"] == [1, 2, 4]
    assert shared["residues_by"]["FULL"] == [10, 11, 13]
    assert shared["residues_by"]["SHORT"] == [20, 21, 23]
    assert shared["dropped_by"] == {"FULL": 1, "SHORT": 0}


def test_an_empty_set_of_pockets_intersects_to_nothing():
    """Rather than raising inside a production run."""
    from flexappeal.refstates.traj import intersect_core

    assert intersect_core({})["positions"] == []


def test_the_alignment_target_is_the_medoid_not_whoever_sorts_first():
    """Modal length stopped deciding anything once every core was the same.

    With `intersect_core` making the lengths uniform, the modal rule fell
    through to sorted order, which on CDK2 is 1F5Q: far enough from the rest
    that fitting all 29 onto it put the mean residual at 0.94 A. Onto the
    medoid it is 0.40.
    """
    import numpy as np

    from flexappeal.refstates import ensembles as ens

    # Three cores on a line: B sits between A and C, so B is the medoid and A
    # sorts first.
    def holder(points):
        return {ens.CORE_KEY: [tuple(p) for p in points]}

    base = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0]])
    built = {"active": ens.MotifEnsemble(state="active", by_code={
        "AAAA": holder(base + np.array([0.0, 0.0, 0.0])),
        "BBBB": holder(base * 1.10),
        "CCCC": holder(base * 1.20),
    })}
    result = ens.align_to_common_frame(built)
    assert result["target"] == "BBBB", result["target"]
    assert result["n_aligned"] == 3


def test_the_two_hop_composition_matches_the_long_way():
    """Own frame to common, then common to run, as one transform.

    `superpose` returns a rotation and two centroids rather than a matrix, and
    the inverse of `(x - a) R^T + b` is `(y - b) R + a`. Composing those by
    hand is four centroids and two rotations to keep straight, so the result is
    checked against doing it in two separate steps.
    """
    import numpy as np

    from flexappeal.refstates.traj import apply_transform, superpose

    module = _analyse_run_module()
    rng = np.random.default_rng(11)
    own = rng.normal(size=(24, 3))
    target = rng.normal(size=(24, 3)) * 1.3 + 4.0
    run = rng.normal(size=(24, 3)) * 0.8 - 2.5

    own_to_common = superpose(own, target)
    run_to_common = superpose(run, target)

    long_way = apply_transform(own, *own_to_common)
    rotation, centre, offset = run_to_common
    long_way = (long_way - offset) @ rotation + centre

    short_way = apply_transform(own, *module._common_to_run(own_to_common,
                                                            run_to_common))
    assert np.allclose(long_way, short_way, atol=1e-10)


# ===========================================================================
#  Section 5.7's bound rows: ligand pose drift and contact persistence
# ===========================================================================


def _bound_fixture(n_frames: int = 4, drift_a: float = 0.0, split: bool = True):
    """A short protein with a small ligand beside it, generated here.

    Public, synthetic coordinates: a helix of alanines and a four-atom ligand.
    `split` writes the ligand as two residues of the same name, which is what
    `abl1_imatinib`'s solvated system does with STI, 37 heavy atoms in one and
    31 hydrogens in the next.
    """
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    protein = []
    for number in range(10, 30):
        residue = topology.add_residue("ALA", chain, resSeq=number)
        for name in ("N", "CA", "C"):
            protein.append(topology.add_atom(name, md.element.carbon, residue))

    ligand_chain = topology.add_chain()
    heavy_res = topology.add_residue("LIG", ligand_chain, resSeq=1)
    for name in ("C1", "C2", "C3", "C4"):
        topology.add_atom(name, md.element.carbon, heavy_res)
    tail = topology.add_residue("LIG", ligand_chain,
                                resSeq=2 if split else 1)
    topology.add_atom("H1", md.element.hydrogen, tail)

    base = []
    for i in range(20):
        anchor = np.array([np.cos(i * 0.5) * 6.0, np.sin(i * 0.5) * 6.0, i * 1.5])
        for offset in ((-1.2, 0, 0), (0, 0, 0), (1.2, 0, 0)):
            base.append(anchor + np.array(offset))
    for offset in ((0, 0, 0), (1.4, 0, 0), (0, 1.4, 0), (0, 0, 1.4), (2.2, 0, 0)):
        base.append(np.array([8.0, 0.0, 14.0]) + np.array(offset))
    base = np.array(base)

    frames = []
    for index in range(n_frames):
        xyz = base.copy()
        # The ligand walks away from the pocket; the protein does not move.
        xyz[len(protein):] += np.array([drift_a * index, 0.0, 0.0])
        frames.append(xyz)
    traj = md.Trajectory(np.array(frames, dtype=np.float32) / 10.0, topology)
    core = list(range(10, 30))
    return traj, core


def test_an_apo_run_reports_no_ligand():
    """Both apo rows go through the same call, so it has to return nothing."""
    from flexappeal.refstates import ligand

    traj, core = _bound_fixture()
    protein_only = traj.atom_slice(
        [a.index for a in traj.topology.atoms if a.residue.name != "LIG"])
    assert ligand.pose_and_contacts(protein_only, core) == {}


def test_a_ligand_split_across_residues_is_one_molecule():
    """`abl1_imatinib` writes STI as two residues, heavy atoms then hydrogens.

    Taking the largest residue returned imatinib's heavy atoms and happened to
    be right, which is the kind of correct that stops being correct when the
    writer changes. Grouped by residue name instead.
    """
    from flexappeal.refstates import ligand

    traj, core = _bound_fixture(split=True)
    out = ligand.pose_and_contacts(traj, core)
    assert out["resname"] == "LIG"
    assert out["n_residues"] == 2, "both residues belong to the ligand"
    assert out["n_heavy_atoms"] == 4, "the hydrogen is not a heavy atom"


def test_a_ligand_that_holds_its_pose_reads_as_zero_drift():
    """And one that walks away does not."""
    from flexappeal.refstates import ligand

    still, core = _bound_fixture(drift_a=0.0)
    assert ligand.pose_and_contacts(still, core)["max_pose_rmsd_a"] == 0.0

    moving, core = _bound_fixture(n_frames=4, drift_a=1.0)
    out = ligand.pose_and_contacts(moving, core)
    assert out["max_pose_rmsd_a"] > 2.5, out["pose_rmsd_a"]
    assert out["final_pose_rmsd_a"] == out["max_pose_rmsd_a"]


def test_the_pose_is_measured_after_superposing_on_the_core():
    """Otherwise the protein's tumbling is reported as the ligand moving.

    The whole complex is rotated rigidly, which changes every coordinate and
    nothing about where the ligand sits in its pocket.
    """
    import numpy as np

    from flexappeal.refstates import ligand

    traj, core = _bound_fixture(n_frames=3, drift_a=0.0)
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                         [np.sin(angle), np.cos(angle), 0.0],
                         [0.0, 0.0, 1.0]])
    for index in range(1, traj.n_frames):
        traj.xyz[index] = traj.xyz[index] @ rotation.T + np.array([3.0, 1.0, 0.0])
    out = ligand.pose_and_contacts(traj, core)
    assert out["max_pose_rmsd_a"] < 1e-3, out["pose_rmsd_a"]


def test_a_contact_counts_a_residue_once_however_many_atoms_are_near():
    """Residues, not atom pairs, which is what "contact persistence" means."""
    from flexappeal.refstates import ligand

    traj, core = _bound_fixture(n_frames=3, drift_a=0.0)
    out = ligand.pose_and_contacts(traj, core)
    contacts = out["persistent_contacts"]
    assert contacts, "the ligand is placed beside the helix, so something touches it"
    numbers = [row["residue"] for row in contacts]
    assert len(numbers) == len(set(numbers)), "a residue should appear once"
    assert all(row["fraction_of_frames"] >= 0.5 for row in contacts)
    assert max(out["n_contacts_per_frame"]) == len(contacts)


def test_a_ligand_written_as_two_residues_still_gets_a_graph():
    """RDKit's proximity bonding does not cross a residue boundary.

    STI written as sequence 1 and 2 gave 41 of imatinib's 72 bonds, every
    hydrogen unbonded, and bond-order assignment failed at every candidate
    charge. The caller saw `None`: no symmetry mapping for the pose RMSD and no
    strain sub-score, with nothing saying the molecule had never been typed.
    """
    from flexappeal.measure import _one_residue

    block = "\n".join([
        "ATOM      1  C1  LIG B   1      10.000  10.000  10.000  1.00  0.00           C",
        "ATOM      2  H1  LIG B   2      10.500  10.500  10.500  1.00  0.00           H",
        "TER",
    ])
    out = _one_residue(block).splitlines()
    assert out[0][22:26] == "   1"
    assert out[1][22:26] == "   1", "the hydrogen joins the heavy atoms"
    # Only the sequence number moves: name, chain and coordinates stay put.
    assert out[1][17:21] == "LIG " and out[1][21] == "B"
    assert "10.500  10.500  10.500" in out[1]
    assert out[2] == "TER", "a record that is not an atom is left alone"


def test_the_analysis_trajectory_still_has_no_ions_in_it():
    """The ions are loaded beside it, not put into it.

    `traj.load` drops the solvent, the ions and the bilayer, and every
    selection downstream runs against the topology it returns. Section 5.5's
    sodium question is answered by loading those positions separately, so this
    asserts the drop is still a drop: the alternative fix, keeping 114
    counter-ions in the analysis topology, would put non-protein atoms in front
    of every other measurement.
    """
    from flexappeal.refstates.traj import NOT_THE_SYSTEM, SODIUM_NAMES

    for name in SODIUM_NAMES:
        if name in ("NA",):
            assert name in NOT_THE_SYSTEM, (
                "sodium must still be dropped from the analysis trajectory")


def test_the_ion_loader_finds_sodium_and_nothing_else(tmp_path: Path):
    """Loaded beside the analysis trajectory, not kept in it.

    Keeping 114 counter-ions in the topology so one measurement can see them
    would put non-protein atoms in front of every other selection. The bilayer
    is handled the same way, as its own structure beside the protein.
    """
    import mdtraj as md
    import numpy as np

    from flexappeal.refstates import traj as traj_mod

    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain, resSeq=1)
    topology.add_atom("CA", md.element.carbon, residue)
    for index, name in enumerate(("NA", "CL", "NA"), start=2):
        ion = topology.add_residue(name, chain, resSeq=index)
        topology.add_atom(name, md.element.sodium if name == "NA"
                          else md.element.chlorine, ion)
    xyz = np.zeros((3, 4, 3), dtype=np.float32)
    xyz[:, 1, 0] = 1.0
    xyz[:, 3, 0] = 2.0
    traj = md.Trajectory(xyz, topology)
    pdb = tmp_path / "s.pdb"
    dcd = tmp_path / "t.dcd"
    traj[0].save_pdb(str(pdb))
    traj.save_dcd(str(dcd))

    found = traj_mod.ion_positions(pdb, dcd)
    assert found is not None
    assert found.shape == (3, 2, 3), "two sodiums, three frames, no chloride"
    assert np.allclose(found[0, 0], [10.0, 0.0, 0.0]), "angstrom, not nanometre"


def test_a_system_with_no_sodium_loads_nothing(tmp_path: Path):
    """And the measurement then omits the ion, exactly as it did before."""
    import mdtraj as md
    import numpy as np

    from flexappeal.refstates import traj as traj_mod

    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain, resSeq=1)
    topology.add_atom("CA", md.element.carbon, residue)
    traj = md.Trajectory(np.zeros((2, 1, 3), dtype=np.float32), topology)
    pdb = tmp_path / "s.pdb"
    dcd = tmp_path / "t.dcd"
    traj[0].save_pdb(str(pdb))
    traj.save_dcd(str(dcd))
    assert traj_mod.ion_positions(pdb, dcd) is None


def _spread_ensembles(active: int, inactive: int):
    """Ensembles whose one motif's references sit at known spacings."""
    import types

    from flexappeal.refstates import ensembles as ens

    def holder(n, offset):
        by_code = {}
        for index in range(n):
            # Three atoms on a line, each reference shifted along it, so the
            # distance to the nearest other reference is the step size.
            step = offset * index
            by_code[f"R{offset}{index:02d}"] = {
                ens.CORE_KEY: [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                               (0.0, 1.0, 0.0)],
                "m1": [(step, 0.0, 0.0), (step + 1.0, 0.0, 0.0),
                       (step, 1.0, 0.0)],
            }
        return types.SimpleNamespace(by_code=by_code)

    return {"active": holder(active, 1), "inactive": holder(inactive, 2)}


def test_the_reference_bands_are_the_spread_each_ensemble_already_has():
    """Section 6.2 wants the 5th to 95th percentile of each ensemble.

    The useful reading is how far apart a state's own references sit for that
    motif, because it says what "close" means before the trajectory is
    compared with it: a line inside the band is as close to the set as the
    set's members are to each other.
    """
    module = _analyse_run_module()
    built = _spread_ensembles(active=5, inactive=5)
    bands = module.reference_spread(built, {"m1": [1, 2, 3]})

    assert set(bands) == {"m1"}
    assert set(bands["m1"]) == {"active", "inactive"}
    # Active references step by 1, inactive by 2, so the inactive set is the
    # wider of the two whatever the percentiles come out at exactly.
    assert bands["m1"]["inactive"]["median_a"] > bands["m1"]["active"]["median_a"]
    assert bands["m1"]["active"]["n"] == 5
    assert bands["m1"]["active"]["low_a"] <= bands["m1"]["active"]["high_a"]


def test_a_small_ensemble_gets_no_band_rather_than_a_narrow_one():
    """Two points do not make a percentile, and an absent band says so.

    Drawn from two references the band would look like a tight ensemble
    instead of a thin one, which is the wrong thing to tell a reader.
    """
    module = _analyse_run_module()
    built = _spread_ensembles(active=2, inactive=5)
    bands = module.reference_spread(built, {"m1": [1, 2, 3]})
    assert "active" not in bands["m1"], "two references cannot give a range"
    assert "inactive" in bands["m1"]


def test_a_motif_no_reference_resolves_gets_no_entry():
    module = _analyse_run_module()
    built = _spread_ensembles(active=5, inactive=5)
    assert module.reference_spread(built, {"absent": [9]}) == {}


def test_the_overlay_uses_the_row_the_ensemble_kept(tmp_path: Path):
    """KLIFS lists a row per chain, and `build` keeps the last one.

    ABL1's fourteen inactive structures arrive as forty-five rows: 6NPV alone
    appears four times, as chains A, B, B, A. `ensembles.build` assigns
    `by_code[code]` once per row so the last wins, while this loop took
    whichever row ranked best, which could be a different chain of the same
    deposition. Its core was then fitted onto a stored core built from the
    other chain, and the ghost landed 0.05 A from where the analysis had
    measured it, which the check refused.

    The two chains here are the same coordinates shifted, so picking the wrong
    one is a transform that is wrong by the shift.
    """
    module = _analyse_run_module()
    first_dir = tmp_path / "chain_a"
    first_dir.mkdir()
    frame, first, numbers = _overlay_fixture(first_dir)
    second_dir = tmp_path / "chain_b"
    second_dir.mkdir()
    _, second, _ = _overlay_fixture(second_dir, shift=(4.0, 0.0, 0.0))

    core = numbers[:8]
    motifs = {"m1": numbers[8:]}
    # Same code twice, as KLIFS would list two chains. The last is what
    # `build` would have kept, so the stored core is built from it.
    rows = {"active": [{"pdb": "TEST", "path": str(first), "chain": "A"},
                       {"pdb": "TEST", "path": str(second), "chain": "A"}]}
    out = tmp_path / "out"
    out.mkdir()
    written = module.write_overlay_bundle(
        frame, "run", out, rows_by_state=rows,
        motifs_for=lambda code: motifs, core_for=lambda code: core,
        core_residues=core, motif_residues=motifs,
        motif_records={"m1": {"nearest_active": "TEST", "d_act_a": 0.0}},
        built=_built_many({"TEST": second}, core))
    assert "active" in written, (
        "the overlay should use the row whose core the ensemble stored")


# ===========================================================================
#  A finished row's bilayer QC does not wait for the queue
# ===========================================================================


def _row_finished(tmp_path: Path, run: str, run_json: str | None) -> bool:
    """`analyse_all.sh`'s `row_finished`, against a run directory we control.

    Sourced rather than executed, which the script's own guard allows, so a
    full analysis pass does not start inside the suite.

    `FLEXAPPEAL_ROOT` is the point of the exercise. The script `cd`s to the
    repository root on the way in, so without the override these tests ran
    against the real `data/md` and their own fixtures were never read: every
    one of them passed against a predicate mutated to return "finished" for
    any run.json at all, failed or not. A test that cannot fail is worse than
    no test, so this passes a root the test owns and the mutation is caught.
    """
    import os
    import subprocess

    root = Path(__file__).resolve().parents[1]
    directory = tmp_path / "data" / "md" / run
    directory.mkdir(parents=True)
    if run_json is not None:
        (directory / "run.json").write_text(run_json)
    return subprocess.run(
        ["bash", "-c",
         f'source "{root}/tools/analyse_all.sh"; row_finished "{run}"'],
        cwd=tmp_path,
        env={**os.environ, "FLEXAPPEAL_ROOT": str(tmp_path)},
        capture_output=True).returncode == 0


def test_a_row_still_producing_is_not_called_finished(tmp_path: Path):
    """No run.json at all, which is every row while it is still going."""
    assert not _row_finished(tmp_path, "b2ar_carazolol", None)


def test_a_failed_row_is_not_called_finished(tmp_path: Path):
    """A row that wrote run.json and failed must not get a bilayer QC.

    The file's presence is not the test: a row that fell over writes one too,
    and computing a bilayer QC on a partial trajectory would put numbers on the
    page that look like a result.
    """
    assert not _row_finished(
        tmp_path, "b2ar_carazolol",
        '{"name": "b2ar_carazolol", "status": "failed:NaN in the forces"}')


def test_a_complete_row_is_finished_whatever_the_spacing(tmp_path: Path):
    """`json.dump` writes `"status": "complete"`, but the match is on the
    pair rather than on one spelling of the whitespace between them."""
    for index, text in enumerate(('{"status": "complete"}',
                                  '{"status":"complete"}',
                                  '{"name": "b2ar_apo", "status"  :   "complete"}')):
        case = tmp_path / f"case{index}"
        case.mkdir()
        assert _row_finished(case, "b2ar_apo", text), text


def test_the_real_finished_row_reads_as_finished():
    """Against this checkout's own b2ar_apo, if it has run.

    The unit tests above build their own run.json, so they would pass against
    a predicate that matched a string this build never writes. This one reads
    what `tools/run_md.py` actually wrote.
    """
    import subprocess

    root = Path(__file__).resolve().parents[1]
    written = root / "data" / "md" / "b2ar_apo" / "run.json"
    if not written.is_file():
        pytest.skip("no finished b2ar_apo in this checkout")
    if '"complete"' not in written.read_text():
        pytest.skip("b2ar_apo did not complete in this checkout")
    assert subprocess.run(
        ["bash", "-c",
         f'source "{root}/tools/analyse_all.sh"; row_finished b2ar_apo'],
        cwd=root, capture_output=True).returncode == 0


# ===========================================================================
#  Section 10's "units everywhere"
# ===========================================================================


def test_the_g6_reason_names_its_unit_in_angstrom():
    """The gate's own sentence is a readout too.

    Section 10's units line is not only about plots. P6 printed its area per
    lipid as `A^2` and its thickness as `A` in the controls report, and
    `gate_six` wrote `A^2` into the reason a raised gate carries, which is the
    sentence somebody reads when the build stops.
    """
    from flexappeal.refstates.membrane import MembraneChunk, gate_six

    raised = gate_six([MembraneChunk(start_ns=0.0, end_ns=10.0,
                                     apl_a2=90.0, thickness_a=38.0)])
    assert raised["open"] is True, "90 A^2 should raise G6"
    assert "Å²" in raised["reason"], raised["reason"]
    assert "A^2" not in raised["reason"], raised["reason"]


def test_the_g6_reason_does_not_round_a_mean_onto_the_window():
    """The sentence that stops a build, at the precision that decides it.

    Same fault as 4bn, in the one readout whose only reader is a person
    deciding whether a bilayer invalidated a result. A mean of 55.04 Å² is
    outside the acceptable window and printed as "55.0" at one decimal, which
    reads as the boundary and so reads as inside, in a sentence whose whole
    job is to say the number is outside.
    """
    from flexappeal.refstates.membrane import MembraneChunk, gate_six

    raised = gate_six([MembraneChunk(start_ns=0.0, end_ns=10.0,
                                     apl_a2=54.96, thickness_a=38.0)])
    assert raised["open"] is True, "54.96 is below the acceptable 55"
    # The window's own bound prints as 55.0, so the measured mean is read off
    # the front of the sentence rather than by searching the whole of it.
    assert raised["reason"].startswith("mean area per lipid 54.96 "), \
        raised["reason"]
    assert raised["mean_apl_a2"] == 54.96


def test_2g1t_deciding_pair_is_what_the_about_page_prints():
    """The one number on a public page that nothing re-derived.

    The About page tells a reader that the published 6.0 Å spine threshold
    calls 2G1T's spine broken when it is intact, and prints 6.14 Å. That number
    was measured once and then written into DECISIONS, FINDINGS and the page,
    and until now nothing checked it again. It is the single reference that
    turns the recalibration in D015 from a preference into a demonstration, so
    it is the one most worth being able to re-derive.

    Skips rather than fails without the KLIFS cache, because it needs the
    deposited structure and the pocket alignment, and a checkout without the
    cache has nothing to measure.
    """
    from pathlib import Path as _Path

    root = _Path(__file__).resolve().parents[1]
    cache = root / "data" / "cache" / "klifs"
    if not (cache / "2g1t.cif.gz").is_file():
        pytest.skip("no cached 2G1T in this checkout")

    import gemmi

    from flexappeal.refstates import (Catalogue, kinase_id, kinase_references,
                                      pocket, refs)
    from flexappeal.refstates.pocket import (R_SPINE_CONTACT_A,
                                             R_SPINE_CONTACT_CRYSTAL_A,
                                             r_spine, r_spine_state)

    catalogue = Catalogue(cache)
    ensembles = kinase_references(catalogue, kinase_id(catalogue))
    entry = next((e for which in ("active", "inactive", "src_like")
                  if ensembles.get(which) is not None
                  for e in ensembles[which].entries
                  if str(e["pdb"]).upper() == "2G1T"), None)
    if entry is None:
        pytest.skip("2G1T is not in this checkout's ensembles")

    structure = gemmi.read_structure(str(refs.fetch_mmcif("2G1T", cache)))
    structure.setup_entities()
    structure.remove_alternative_conformations()
    chain = next(c for c in structure[0]
                 if c.name == str(entry.get("chain") or "A"))
    state = r_spine_state(chain, r_spine(pocket.pocket_map(
        catalogue, entry.get("structure_ID")), chain))

    widest = max(p["distance_a"] for p in state["pairs"] if "distance_a" in p)
    assert widest == pytest.approx(6.14, abs=0.01), state["pairs"]
    # The point of the number: the two thresholds disagree about this one
    # structure, and that disagreement is the argument for having moved it.
    assert widest < R_SPINE_CONTACT_A, "should be assembled at the build's 7.7"
    assert widest > R_SPINE_CONTACT_CRYSTAL_A, "should be broken at 6.0"
    assert state["state"] == "assembled"


# ===========================================================================
#  The documents that quote P4 agree with what P4 computes
# ===========================================================================


def _p4_now() -> dict[str, float]:
    """P4's four scores, from the artefact `controls` writes.

    Not recomputed here: P4 re-analyses each start model against its whole
    reference ensemble, which is tens of seconds and a network cache. The
    artefact is the same numbers the command printed, written every run.
    """
    import json

    root = Path(__file__).resolve().parents[1]
    artefact = root / "data" / "analysis" / "controls.json"
    if not artefact.is_file():
        pytest.skip("no controls.json in this checkout")
    # Stale against what P4 actually depends on, which is not the database.
    # P4 scores the start structures at frame 0 against the reference
    # ensembles; it never reads a trajectory. Comparing against
    # `refstates.db` made this skip every time the streaming loop re-ingested,
    # which is every half hour, so the check was almost always absent exactly
    # while production was running.
    age = artefact.stat().st_mtime
    inputs = [root / "data" / "systems", root / "tools" / "analyse_run.py"]
    inputs += [root / "flexappeal" / "refstates" / f"{name}.py" for name in
               # `setup` is where the family wiring lives since it moved out of
               # `analyse_run.py`, so it determines P4 as much as the script
               # that used to hold it.
               ("traj", "state", "ensembles", "kinase", "gpcr", "pocket",
                "setup")]
    newer = []
    for item in inputs:
        if item.is_dir():
            newer += [f for f in item.rglob("*.pdb") if f.stat().st_mtime > age]
        elif item.is_file() and item.stat().st_mtime > age:
            newer.append(item)
    if newer:
        pytest.skip("controls.json predates "
                    f"{newer[0].name}, which determines P4; re-run "
                    "`cli controls`")
    rows = json.loads(artefact.read_text())["controls"]["P4"]
    return {r["name"]: round(float(r["whole_s"]), 3)
            for r in rows if r.get("whole_s") is not None}


def test_the_documents_quoting_p4_agree_with_it():
    """FINDINGS 1.4 and GATE_OPEN both print P4's four numbers.

    They drifted on 11 October: the superposition core and the alignment target
    changed, `abl1_apo` moved +0.305 to +0.304 and `b2ar_carazolol` -0.934 to
    -0.923, and **both** documents kept the old pair, as did `controls.json`,
    which `controls` was not writing unless asked. The same two numbers wrong
    in three places is what a copied figure looks like, and none of the three
    is generated from the others.

    It reads the scores where they are printed as a signed three-decimal number
    and asks whether any is one P4 does not currently produce, rather than
    checking the documents mention the run at all, which would fire on every
    rewording.
    """
    import re

    current = _p4_now()
    if len(current) < 4:
        pytest.skip("not all four rows are scored in this checkout")
    allowed = {f"{v:+.3f}" for v in current.values()}

    root = Path(__file__).resolve().parents[1]
    checked = 0
    for name, pattern in (
            ("FINDINGS.md",
             r"\|\s*(?:abl1|b2ar)_\w+[^|]*\|[^|]*\|\s*\*\*([+\u2212-]0\.\d{3})\*\*"),
            # Exactly the four comma-separated scores, not the rest of the
            # cell: it also names the +0.300 pass band, and a looser capture
            # read that as a fifth score and failed on it.
            ("GATE_OPEN.md",
             r"P4 each start scores as its own state \|[^|]*?4 / 4:\s*"
             r"((?:[+\u2212-]0\.\d{3}(?:,\s*)?){4})")):
        text = (root / name).read_text()
        found = re.findall(pattern, text)
        if not found:
            continue
        quoted = []
        for chunk in found:
            quoted += re.findall(r"[+\u2212-]0\.\d{3}", chunk)
        quoted = [q.replace("\u2212", "-") for q in quoted]
        if not quoted:
            continue
        checked += 1
        stale = sorted({q for q in quoted if q not in allowed})
        assert not stale, (
            f"{name} quotes {stale} as P4 scores; P4 currently produces "
            f"{sorted(allowed)}. Re-derive the numbers in that document.")
    if not checked:
        pytest.skip("neither document prints P4's scores in a known shape")


def test_a_real_ligand_gets_its_symmetry_mapping():
    """Section 5.7's symmetry correction, on a real ligand rather than a fixture.

    `pose_and_contacts` had only ever run on synthetic fixtures: both bound
    rows are last in the queue. Rehearsed on the built systems before they
    produce, because an uncorrected pose RMSD does not fail, it just overstates
    drift, and the row would have been finished and written up before anyone
    looked at the `symmetry_corrected` flag.

    It also pins the distinction the rehearsal turned up. The **start models**
    carry heavy atoms only, so the graph inference has nothing to work with and
    the mapping comes back empty with a note. The **solvated** systems carry
    hydrogens, in a second residue of the same name, and that is what the
    production path analyses. Grouping by name rather than taking the largest
    residue is what makes the second case work, and carazolol splits the same
    way imatinib does: 22 heavy atoms and 22 hydrogens.

    Skips without a built system, which a clean checkout has: the systems are
    gigabytes and are not committed.
    """
    import tempfile

    root = Path(__file__).resolve().parents[1]
    solvated = root / "data" / "md" / "abl1_imatinib" / "solvated.pdb"
    if not solvated.is_file():
        pytest.skip("no built abl1_imatinib system in this checkout")

    import mdtraj as md

    from flexappeal.refstates import ligand as ligand_mod, traj

    with tempfile.TemporaryDirectory() as work:
        frame = Path(work) / "one.dcd"
        md.load(str(solvated)).save_dcd(str(frame))
        loaded = traj.load(solvated, frame)

    found = ligand_mod.ligand_residues(loaded.topology)
    names = {loaded.topology.residue(i).name for i in found}
    assert names == {"STI"}, f"expected imatinib alone, got {names}"
    assert len(found) == 2, ("the solvated system splits the ligand into a "
                             "heavy residue and a hydrogen residue; grouping "
                             "by name is what rejoins them")

    protein = [r.index for r in loaded.topology.residues if r.is_protein]
    out = ligand_mod.pose_and_contacts(loaded, protein[:200])
    assert out, "a bound row produced no ligand panel"
    assert out["resname"] == "STI"
    assert out["n_heavy_atoms"] == 37, out["n_heavy_atoms"]
    assert out["symmetry_corrected"] is True, out.get("note")
    # The pocket, not just a count. Asserted on residue names and distances
    # rather than numbers: this loads the solvated topology directly, and
    # OpenMM renumbers from 1 on solvation, so the catalytic lysine that is
    # Lys271 in the deposit is Lys16 here. Restoring the deposited numbering is
    # what `--start` does in the real pipeline, and it is not this test's
    # subject.
    contacts = out["persistent_contacts"]
    assert len(contacts) == 22, len(contacts)
    names = {row["name"] for row in contacts}
    assert {"LYS", "GLU"} <= names, sorted(names)
    assert min(row["closest_a"] for row in contacts) < 3.5


def _wrapped_complex(wrap: bool):
    """A two-residue protein with a ligand beside it, optionally wrapped.

    Generated here rather than committed. The ligand sits 3 Å from the protein
    in reality; when `wrap` is set it is moved one box length in y, which is
    what an integrator does to a molecule that drifts through the edge and is
    exactly what the first real `b2ar_carazolol` frames looked like.
    """
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    protein_atoms = []
    for index in range(6):
        residue = topology.add_residue("ALA", chain, resSeq=index + 1)
        for atom_name in ("N", "CA", "C", "O"):
            protein_atoms.append(
                topology.add_atom(atom_name, md.element.carbon
                                  if atom_name != "O" else md.element.oxygen,
                                  residue))
    ligand_chain = topology.add_chain()
    ligand_residue = topology.add_residue("LIG", ligand_chain, resSeq=900)
    for index in range(4):
        topology.add_atom(f"C{index}", md.element.carbon, ligand_residue)

    box = np.array([4.0, 4.0, 4.0], dtype=np.float32)
    protein = np.zeros((len(protein_atoms), 3), dtype=np.float32)
    protein[:, 0] = np.arange(len(protein_atoms)) * 0.15
    protein[:, 1] = 2.0
    protein[:, 2] = 2.0
    ligand = np.zeros((4, 3), dtype=np.float32)
    ligand[:, 0] = np.arange(4) * 0.15 + 0.2
    ligand[:, 1] = 2.0 + 0.30          # 3 A from the protein in y
    ligand[:, 2] = 2.0
    if wrap:
        ligand[:, 1] -= box[1]          # through the edge, as an integrator does

    xyz = np.stack([np.vstack([protein, ligand])] * 3)
    traj = md.Trajectory(xyz.astype(np.float32), topology)
    traj.unitcell_lengths = np.tile(box, (3, 1))
    traj.unitcell_angles = np.full((3, 3), 90.0, dtype=np.float32)
    return traj


def test_a_wrapped_ligand_is_imaged_back_onto_its_pocket():
    """Section 5.7's readings both depend on the ligand being in the right image.

    The first real `b2ar_carazolol` frames reported carazolol 110 Å from its
    own receptor, holding **zero** contacts, with a pose RMSD of 5.4 Å. The
    ligand had been wrapped through the edge of the box: it sat at y = 0.18 nm
    against a receptor at 12.83 in a box 13.11 tall, so the two were touching
    and the arithmetic was not. Imaged about the receptor the same frames give
    a pose RMSD of 0.45 Å and 15 contacts.

    The failure is silent and reads as a result, which is why it has a test: a
    lost pose and no contacts is exactly what a ligand that had genuinely come
    out of its pocket would look like.
    """
    from flexappeal.refstates import ligand as ligand_mod

    core = list(range(6))
    here = ligand_mod.pose_and_contacts(_wrapped_complex(wrap=False), core)
    there = ligand_mod.pose_and_contacts(_wrapped_complex(wrap=True), core)

    assert here["persistent_contacts"], "the unwrapped case should hold contacts"
    assert len(there["persistent_contacts"]) == len(here["persistent_contacts"]), (
        "a wrapped ligand must read the same as the same ligand unwrapped; "
        f"got {len(there['persistent_contacts'])} against "
        f"{len(here['persistent_contacts'])}")
    assert there["max_pose_rmsd_a"] == pytest.approx(here["max_pose_rmsd_a"],
                                                     abs=1e-6)
    # The note this fixture does carry is about symmetry: four carbons in a
    # line have no graph to infer, which is a property of the fixture. What
    # must not appear is a note about imaging, which would mean the correction
    # declined to run.
    assert "image" not in (there.get("note") or ""), there.get("note")
    assert "box" not in (there.get("note") or ""), there.get("note")


def test_a_trajectory_with_no_box_says_so_rather_than_imaging():
    """Omission over a stand-in: no box means no imaging, and the note says it."""
    from flexappeal.refstates import ligand as ligand_mod

    traj = _wrapped_complex(wrap=False)
    traj.unitcell_lengths = None
    traj.unitcell_angles = None
    out = ligand_mod.pose_and_contacts(traj, list(range(6)))
    assert "no box" in (out.get("note") or ""), out.get("note")


def test_a_wrapped_sodium_is_imaged_onto_the_receptor(tmp_path: Path):
    """Section 5.5 asks whether a Na+ sits within 3.5 A of D2.50.

    The ions were loaded as raw coordinates. An ion is free to wrap through
    the edge of the box, and a wrapped sodium reads as tens of angstrom from a
    pocket it is sitting in, so the question could only ever be answered one
    way. On `b2ar_apo` the raw coordinates put the nearest sodium to one
    aspartate at 42.78 A and the minimum image puts it at 2.18.

    Per atom, unlike the ligand: an ion is a single atom with no geometry to
    preserve, so each belongs in whichever copy sits nearest the receptor. A
    molecule imaged that way would come apart, which is why
    `pose_and_contacts` shifts the ligand as a whole.
    """
    import mdtraj as md
    import numpy as np

    from flexappeal.refstates import traj as traj_mod

    topology = md.Topology()
    chain = topology.add_chain()
    for index in range(4):
        residue = topology.add_residue("ALA", chain, resSeq=index + 1)
        topology.add_atom("CA", md.element.carbon, residue)
    ion_chain = topology.add_chain()
    ion_residue = topology.add_residue("NA", ion_chain, resSeq=500)
    topology.add_atom("NA", md.element.sodium, ion_residue)

    box = np.array([4.0, 4.0, 4.0], dtype=np.float32)
    protein = np.zeros((4, 3), dtype=np.float32)
    protein[:, 0] = np.arange(4) * 0.15
    protein[:, 1] = 2.0
    protein[:, 2] = 2.0
    # 2 A from the first residue in reality, then wrapped through the edge.
    ion = np.array([[0.0, 2.0 - box[1] + 0.20, 2.0]], dtype=np.float32)

    xyz = np.stack([np.vstack([protein, ion])] * 2)
    built = md.Trajectory(xyz.astype(np.float32), topology)
    built.unitcell_lengths = np.tile(box, (2, 1))
    built.unitcell_angles = np.full((2, 3), 90.0, dtype=np.float32)

    structure = tmp_path / "system.pdb"
    trajectory = tmp_path / "run.dcd"
    built[0].save_pdb(str(structure))
    built.save_dcd(str(trajectory))

    positions = traj_mod.ion_positions(structure, trajectory)
    assert positions is not None, "the sodium should have been found"
    # In angstrom, and beside the protein rather than a box away from it.
    nearest = float(np.linalg.norm(positions[0][0] - protein[0] * 10.0))
    assert nearest < 5.0, (
        f"a wrapped sodium 2 A from the protein read as {nearest:.1f} A away")


def _chain_trajectory(break_it: bool):
    """A short chain, optionally with one residue wrapped through the box."""
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    for index in range(8):
        residue = topology.add_residue("ALA", chain, resSeq=index + 1)
        topology.add_atom("CA", md.element.carbon, residue)
    coords = np.zeros((8, 3), dtype=np.float32)
    coords[:, 0] = np.arange(8) * 0.38          # a peptide bond, in nm
    if break_it:
        coords[4:, 0] -= 4.0                     # the tail wrapped one box
    xyz = np.stack([coords] * 2)
    traj = md.Trajectory(xyz.astype(np.float32), topology)
    traj.unitcell_lengths = np.tile(np.array([4.0, 4.0, 4.0]), (2, 1))
    traj.unitcell_angles = np.full((2, 3), 90.0, dtype=np.float32)
    return traj


def test_a_whole_protein_is_not_reported_as_broken():
    module = _analyse_run_module()
    assert module.report_broken_chain(_chain_trajectory(break_it=False)) is False


def test_a_protein_split_across_the_box_is_reported(capsys):
    """The one geometry nothing images, so it is checked instead.

    Four cross-molecule distances in this build needed the minimum image and
    none of them announced itself; each was found by a number that looked
    wrong. The protein is the case where nothing images anything, because
    OpenMM keeps a molecule whole, and every state score depends on that. A
    split chain does not fail: it measures the long way round the box and
    returns a plausible number.
    """
    module = _analyse_run_module()
    assert module.report_broken_chain(_chain_trajectory(break_it=True)) is True
    said = capsys.readouterr().out
    assert "broken" in said and "periodic boundary" in said, said


def _p5_database(tmp_path: Path, *, status: str, length_ns: float,
                 scores: list[tuple[float, float]]):
    """A one-motif inactive-start run with the given per-frame scores."""
    from flexappeal.refstates import store

    database = tmp_path / "refstates.db"
    store.ingest_analysis(
        database,
        {"name": "f", "kind": "gpcr",
         "core": {"n_residues": 10, "source": "fixture"},
         "motif_residues": {"TM6": list(range(260, 291))},
         "ensembles": {"active": {"codes": ["3P0G"]},
                       "inactive": {"codes": ["2RH1"]}},
         "frames": [{"time_ns": t, "n_conflicts": 0, "motifs": [
             {"motif": "TM6", "d_act_a": 3.0, "d_inact_a": 1.0, "s": s,
              "label": "inactive"}]} for t, s in scores],
         "notes": []},
        system_id="b2ar", run_id="b2ar_carazolol", family="gpcr",
        environment="membrane",
        run_fields={"start_state": "inactive", "status": status,
                    "length_ns": length_ns})
    return database


def test_p5_does_not_pass_a_run_that_is_still_going(tmp_path: Path):
    """The final 20 ns of a running row is not final.

    P5 reported PASS on `b2ar_carazolol` 8.4 ns into a 50 ns run. At that
    length the window is the whole trajectory and every frame of it is still
    near the structure the run started from, so the control was measuring the
    start state, which is what P4 measures. Section 9 puts the window at the
    end precisely to make the answer independent of the start.
    """
    from flexappeal.refstates.cli import control_p5

    database = _p5_database(
        tmp_path, status="running", length_ns=8.4,
        scores=[(t / 10.0, -0.34) for t in range(84)])
    rows = control_p5(database)
    assert len(rows) == 1, rows
    assert rows[0]["status"] == "pending", rows[0]
    assert "passed" not in rows[0], rows[0]
    assert "not final until it stops" in rows[0]["reason"], rows[0]


def test_p5_judges_a_finished_run(tmp_path: Path):
    """And still answers once the run has stopped."""
    from flexappeal.refstates.cli import control_p5

    database = _p5_database(
        tmp_path, status="complete", length_ns=50.0,
        scores=[(t / 10.0, -0.34) for t in range(500)])
    rows = control_p5(database)
    assert rows[0]["status"] == "ok", rows[0]
    assert rows[0]["passed"] is True, rows[0]
    # The window is the final 20 ns of 50, not the whole run.
    assert rows[0]["frames_in_window"] < 500, rows[0]["frames_in_window"]




# ---------------------------------------------------------------------------
#  Scoring inside a generated bundle
# ---------------------------------------------------------------------------


def test_the_scoring_path_needs_nothing_a_bundle_does_not_have():
    """The architectural guard for the whole arrangement.

    Resolving which residues a motif is needs KLIFS or GPCRdb, their deposited
    structures and gemmi. Scoring a frame against the answer needs none of
    them, and that asymmetry is why the work can be split: the resolving half
    runs where the bundle is built and the scoring half runs on a machine with
    an OpenMM environment and nothing else.

    If an import creeps into the scoring half, the symptom is a bundle that
    builds and then fails on somebody's laptop, so it is checked here rather
    than discovered there.
    """
    import subprocess
    import sys

    script = """
import sys
class Block:
    BANNED = {"gemmi", "requests", "tenacity", "sse_comparison"}
    def find_module(self, name, path=None):
        if name.split(".")[0] in self.BANNED:
            raise ImportError(name + " is not in a generated bundle")
        return None
sys.meta_path.insert(0, Block())
from flexappeal.refstates import traj, state, sse
assert hasattr(traj, "analyse") and hasattr(traj, "deposited_numbering")
assert hasattr(state, "smooth") and hasattr(sse, "secondary_structure")
assert hasattr(sse, "profiles_from_json")
print("ok")
"""
    result = subprocess.run([sys.executable, "-c", script],
                            capture_output=True, text=True,
                            cwd=str(Path(__file__).resolve().parents[1]))
    assert result.returncode == 0, result.stderr[-800:]
    assert "ok" in result.stdout


def test_the_bundle_carries_the_ensemble_and_the_scoring_only_when_asked():
    """Off by default, and a run that did not ask pays nothing for it."""
    from flexappeal import bundle, schema

    structure = (b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  "
                 b"1.00  0.00           N\nEND\n")

    plain = schema.validate({"job_name": "plain"}).config
    files = bundle.unpack(bundle.build(plain, structure).content)
    assert not [n for n in files if "refstates" in n]

    asked = schema.validate({"job_name": "asked", "use_refstates": "on",
                             "refstates_family": "kinase"}).config
    files = bundle.unpack(bundle.build(asked, structure,
                                       refstates=b'{"version":1}').content)
    assert files[bundle.REFSTATES_MEMBER] == b'{"version":1}'
    for name in bundle.VENDORED_REFSTATES:
        vendored = files[f"fxalib/refstates/{name}"]
        assert vendored == (Path(__file__).resolve().parents[1] / "flexappeal"
                            / "refstates" / name).read_bytes(), (
            f"{name} should travel byte for byte")


def test_reference_profiles_survive_a_round_trip_through_json():
    """A Counter is a dict after JSON and an int key is a string.

    The ribbon simply does not draw when they come back as plain objects, with
    every score around it correct, which is the quiet kind of failure.
    """
    import json

    from collections import Counter

    from flexappeal.refstates import sse

    original = {"active": {"A-loop": {3: Counter({"H": 4, "E": 1})}},
                "inactive": {"A-loop": {3: Counter({"E": 5})}}}
    restored = sse.profiles_from_json(json.loads(json.dumps(original)))

    assert restored["active"]["A-loop"][3].most_common(1) == [("H", 4)]
    assert isinstance(restored["inactive"]["A-loop"][3], Counter)
    assert sse.profiles_from_json({}) == {}
    assert sse.profiles_from_json(None) == {}


def test_a_structure_no_catalogue_knows_is_refused_at_build_time():
    """With a sentence, rather than hours later on somebody's machine."""
    from flexappeal.refstates import prepare

    # Generated here: a short peptide is neither a kinase nor a receptor.
    peptide = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "peptide12.pdb"
    if not peptide.is_file():
        pytest.skip("tests/fixtures/peptide12.pdb is missing")

    with pytest.raises(prepare.PrepareError) as caught:
        prepare.resolve(peptide, {"refstates_family": "auto"})
    assert "kinase or a receptor" in str(caught.value)


def test_a_reference_with_no_resolution_is_left_without_one():
    """KLIFS returns it as a string and the table formats with `.2f`.

    A reference whose resolution cannot be read is a reference with no
    resolution, not one at zero angstrom, and the first version of this crashed
    the whole references table rather than omitting a cell.
    """
    from flexappeal.refstates import prepare

    assert prepare._resolution({"resolution": "1.91"}) == {"resolution": 1.91}
    assert prepare._resolution({"resolution": 2.4}) == {"resolution": 2.4}
    assert prepare._resolution({"resolution": None}) == {}
    assert prepare._resolution({"resolution": "NOT DETERMINED"}) == {}
    assert prepare._resolution({}) == {}


def test_the_state_map_cloud_is_carried_rather_than_recomputed():
    """It is a property of the references, not of the trajectory.

    Measuring each reference on its published criterion plane needs the
    catalogue's own structures, and for a kinase the vendored annotator's atom
    lookups. Both are build-time things, so the cloud travels with the ensemble
    and `analyse.py` copies it into the results file.
    """
    source = (Path(__file__).resolve().parents[1] / "flexappeal" / "runtime"
              / "analyse.py.j2").read_text()
    block = source[source.index("def reference_states("):]
    block = block[: block.index("\ndef ")]
    assert 'payload.get("reference_structures")' in block
    assert "dfg_d_alphac" not in block, (
        "the bundle must not measure the references itself: it has neither "
        "their structures nor the annotator")


def test_the_tm7_motif_is_named_for_the_receptor_class():
    """7.49 to 7.53 is NPxxY in class A and PxxG in class B.

    The span is the same structural place in either, because a generic number
    is a position rather than a sequence, so only the name changes. Everything
    this build scored until GLP-1R was class A and the class A name was a
    constant in the motif list, which handed a class B receptor a motif it does
    not have.
    """
    from flexappeal.refstates.gpcr import tm7_motif_name

    assert tm7_motif_name("Class A (Rhodopsin)") == "NPxxY"
    assert tm7_motif_name("Class B1 (Secretin)") == "PxxG"
    assert tm7_motif_name("Class B2 (Adhesion)") == "PxxG"
    # A class this build has not met is named by the numbers themselves, which
    # are true of any of them, rather than by class A's motif.
    assert tm7_motif_name("Class C (Glutamate)") == "TM7 7.49-7.53"
    assert tm7_motif_name("") == "TM7 7.49-7.53"
    assert tm7_motif_name(None) == "TM7 7.49-7.53"


def test_the_state_map_labels_its_axis_from_the_payload():
    """A file carries the name its own class uses, and the panel follows it.

    A results file written before the names were class-aware has no such key,
    and every one of those is class A, so the fallback is NPxxY.
    """
    from flexappeal.refstates import panels

    base = {"kind": "gpcr", "frames": [], "reference_structures": [
        {"pdb": "6XOX", "source_state": "active",
         "metrics": {"tm3_tm6": 14.0, "npxxy_to_inactive": 2.1}}]}

    class_b = dict(base, report={"map_axis": {"motif": "PxxG"}})
    assert "PxxG RMSD" in panels.state_map(class_b)["layout"]["yaxis"]["title"]["text"]

    class_a = dict(base, report={"map_axis": {"motif": "NPxxY"}})
    assert "NPxxY RMSD" in panels.state_map(class_a)["layout"]["yaxis"]["title"]["text"]

    # No key at all: an older file, and older files are class A.
    assert "NPxxY RMSD" in panels.state_map(base)["layout"]["yaxis"]["title"]["text"]


def test_the_camera_is_measured_the_same_way_in_both_places():
    """A run of your own should open like the library's runs do.

    The packer measures two things for the four library rows: the bilayer
    normal, from the axis the phosphorus atoms are least spread along, and
    which end of it is extracellular, from helix 8 lying on the intracellular
    face. A bundle measures the same two, and the first version of it did not:
    it took the protein's centroid against the phosphorus mean, which is nearly
    the same point, so the sign was arbitrary and a receptor could open upside
    down with the camera rotated anywhere about the normal.

    Built here rather than committed: a slab of phosphorus either side of a
    protein whose helix 8 is on the low side of it.
    """
    import numpy as np
    import mdtraj as md

    from flexappeal.refstates.membrane import bilayer_orientation

    rng = np.random.default_rng(0)
    # A bilayer across z: the leaflets are far wider in the plane than the
    # pair is thick, which is the shape the normal is found from. The real
    # receptor runs measure about 38 and 40 angstrom of spread in the plane
    # against 19 across, and the test below is the same ratio.
    plane = rng.uniform(-60, 60, size=(120, 2))
    z = np.where(np.arange(120) < 60, -18.0, 18.0)
    phosphorus = np.column_stack([plane, z])

    lines = []
    for serial, xyz in enumerate(phosphorus, start=1):
        lines.append(
            f"HETATM{serial:>5d} P    POP A{serial:>4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           P")
    for offset in range(4):
        serial = 121 + offset
        lines.append(
            f"ATOM  {serial:>5d} CA   ALA B{offset + 1:>4d}    "
            f"{0.0:8.3f}{0.0:8.3f}{offset * 2.0:8.3f}  1.00  0.00           C")
    import tempfile
    from pathlib import Path as _Path

    path = _Path(tempfile.mkdtemp()) / "slab.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    traj = md.load(str(path))

    # The protein spans z = 0 to 6 and is widest along x; helix 8 is its lowest
    # part, so up is +z and the camera looks across x.
    protein = np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 2.0],
                        [0.0, 5.0, 4.0], [0.0, 0.0, 6.0]])
    helix8 = np.array([[0.0, 0.0, -4.0], [1.0, 0.0, -4.0]])

    found = bilayer_orientation(traj, 0, protein_xyz=protein, helix8_xyz=helix8)
    assert found["normal_axis"] == 2, found
    assert found["up"][2] > 0, "helix 8 is below, so up is the positive normal"
    # And the camera looks across the membrane rather than down it.
    assert abs(found["along"][2]) < 1e-6, found["along"]

    # Flip helix 8 to the other side and the view flips with it.
    flipped = bilayer_orientation(traj, 0, protein_xyz=protein,
                                  helix8_xyz=helix8 + np.array([0, 0, 14.0]))
    assert flipped["up"][2] < 0

    # Without helix 8 there is no way to tell one end from the other, and a
    # guess is refused rather than made.
    assert bilayer_orientation(traj, 0, protein_xyz=protein) == {}
