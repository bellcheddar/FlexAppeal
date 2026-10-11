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


def test_receptor_tilt_is_folded_into_the_first_quadrant():
    along_z = np.array([[0.0, 0.0, float(i)] for i in range(10)])
    assert membrane.receptor_tilt(along_z) == pytest.approx(0.0, abs=1e-6)
    flipped = np.array([[0.0, 0.0, -float(i)] for i in range(10)])
    assert membrane.receptor_tilt(flipped) == pytest.approx(0.0, abs=1e-6)


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


def test_the_pages_render_without_a_database(tmp_path: Path):
    """The feature is additive: a deployment without the file must not 500."""
    from flexappeal.webapp import create_app

    app = create_app()
    app.config["REFSTATES_DB"] = str(tmp_path / "absent.db")
    client = app.test_client()
    assert client.get("/refstates/").status_code == 200
    assert client.get("/refstates/run/nothing").status_code == 404


def test_the_finding_is_generated_from_the_data(tmp_path: Path):
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    with app.test_request_context():
        from flexappeal.refstates.web import finding

        text = finding("r")
    assert "xDFG" in text
    assert "towards active" in text


def test_the_run_page_shows_every_motif(tmp_path: Path):
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    client = app.test_client()
    page = client.get("/refstates/run/r").get_data(as_text=True)
    assert "xDFG" in page and "hinge" in page
    series = client.get("/refstates/run/r/frames.json").get_json()
    assert set(series["motifs"]) == {"s:xDFG", "s:hinge"}
    assert series["motifs"]["s:xDFG"]["s"] == [0.6, 0.8]


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


def test_the_dssp_ribbon_is_decimated_and_keeps_its_gaps(tmp_path: Path):
    """An undecimated per-frame heatmap is what made one page 24 MB before.

    The gap check matters as much: an unassigned residue has to arrive as null
    so the ribbon leaves it blank, rather than as coil, which would read as a
    measurement that was never taken.
    """
    import json as _json

    from flexappeal.refstates.web import MAX_RIBBON_COLUMNS
    from flexappeal.webapp import create_app

    analysis = tmp_path / "analysis"
    analysis.mkdir()
    frames = ["HHH-E"] * (MAX_RIBBON_COLUMNS * 3)
    (analysis / "r.json").write_text(_json.dumps({
        "dssp": {"n_frames": len(frames), "residues": [10, 11, 12, 13, 14],
                 "frames": frames, "motifs": {}}}))

    app = create_app()
    client = app.test_client()
    with app.app_context():
        import flexappeal.refstates.web as web

        original = web.analysis_file
        web.analysis_file = lambda run_id: analysis / f"{run_id}.json"
        try:
            payload = client.get("/refstates/run/r/dssp.json").get_json()
        finally:
            web.analysis_file = original

    assert len(payload["frames"]) <= MAX_RIBBON_COLUMNS
    assert payload["n_frames_total"] == len(frames)
    assert payload["stride"] >= 3
    # The unassigned residue survives as the blank marker, not as coil.
    assert payload["frames"][0][3] == "-"


def test_the_viewer_route_serves_only_its_two_suffixes(tmp_path: Path):
    """The route must not be talked into returning anything else in there.

    The analysis directory also holds each run's JSON, and a route that served
    by filename would hand it over. Served by suffix from a fixed allowlist
    instead, with the run id checked against the app's own pattern.
    """
    from flexappeal.webapp import create_app

    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (analysis / "r.pdb").write_text("ATOM\nEND\n")
    (analysis / "r.xtc").write_bytes(b"\x00\x01")
    (analysis / "r.json").write_text('{"secret": 1}')

    app = create_app()
    client = app.test_client()
    import flexappeal.refstates.web as web

    original_root = app.root_path
    app.root_path = str(tmp_path / "flexappeal")
    (tmp_path / "flexappeal").mkdir(exist_ok=True)
    (tmp_path / "data").mkdir(exist_ok=True)
    try:
        # The route resolves data/analysis relative to the app root's parent.
        import shutil

        shutil.copytree(analysis, tmp_path / "data" / "analysis",
                        dirs_exist_ok=True)
        assert client.get("/refstates/run/r/viewer.pdb").status_code == 200
        assert client.get("/refstates/run/r/viewer.xtc").status_code == 200
        # Not in the allowlist, so not served even though the file is there.
        assert client.get("/refstates/run/r/viewer.json").status_code == 404
        # And an id that is not an id.
        assert client.get("/refstates/run/..%2F..%2Fetc/viewer.pdb").status_code == 404
    finally:
        app.root_path = original_root


def test_the_viewer_panel_is_hidden_without_its_files(tmp_path: Path):
    """A viewer that fails to load is worse than one that is not offered."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    page = app.test_client().get("/refstates/run/r").get_data(as_text=True)
    assert "molstar-viewer" not in page
    # The rest of the page still renders.
    assert "xDFG" in page


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

    `receptor_tilt` measures the principal axis against z, so a helix built
    along z and tilted by a known angle in the xz plane must report it.
    """
    import numpy as np

    straight = np.array([[0.0, 0.0, 1.5 * i] for i in range(20)])
    for degrees in (0.0, 15.0, 30.0, 45.0):
        angle = np.radians(degrees)
        rotation = np.array([[np.cos(angle), 0.0, np.sin(angle)],
                             [0.0, 1.0, 0.0],
                             [-np.sin(angle), 0.0, np.cos(angle)]])
        tilted = straight @ rotation.T
        assert membrane.receptor_tilt(tilted) == pytest.approx(degrees, abs=1e-4)


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
    resume = source.index("if checkpoint.is_file():")
    heating = source.index("heating to")
    assert resume < heating, "equilibration is reached before the resume check"
    assert "loadCheckpoint" in source
    assert "already at" in source, "a completed run must report and stop"


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

REFSTATES_PAGES = ["/refstates/", "/refstates/about"]


def _page_client(tmp_path: Path):
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    return app.test_client()


def test_every_referenced_asset_resolves(tmp_path: Path):
    """A 404 on a script or stylesheet is the usual cause of a console error."""
    import re

    client = _page_client(tmp_path)
    missing = []
    for url in REFSTATES_PAGES + ["/refstates/run/r", "/refstates/references/s"]:
        response = client.get(url)
        if response.status_code != 200:
            continue
        html = response.get_data(as_text=True)
        for ref in set(re.findall(r'(?:src|href)="(/static/[^"]+)"', html)):
            if client.get(ref).status_code != 200:
                missing.append(f"{url} -> {ref}")
    assert not missing, f"referenced but not served: {missing}"


def test_the_pages_meet_the_markup_rules_axe_would_check(tmp_path: Path):
    """The subset of axe-core's serious rules that live in the markup.

    Not a substitute for axe-core, and FINDINGS says so. These are the ones a
    template can get wrong and a server can see: a table without header scope,
    an image without alt text, a page without a language.
    """
    import re

    client = _page_client(tmp_path)
    problems = []
    for url in REFSTATES_PAGES + ["/refstates/run/r", "/refstates/references/s"]:
        response = client.get(url)
        if response.status_code != 200:
            continue
        html = response.get_data(as_text=True)
        if not re.search(r'<html[^>]+lang=', html):
            problems.append(f"{url}: <html> has no lang")
        for image in re.findall(r"<img[^>]*>", html):
            if "alt=" not in image:
                problems.append(f"{url}: an img has no alt")
        # Every th in these pages is a header for a row or a column, and
        # axe-core flags one without a scope in a table that has both.
        # `<th(?=[\s>])` so the pattern does not also match `<thead>`, which
        # it did on the first run.
        for header in re.findall(r"<th(?=[\s>])[^>]*>", html):
            if "scope=" not in header:
                problems.append(f"{url}: a th has no scope: {header}")
        if "<table" in html and "<thead" not in html and "scope=\"row\"" not in html:
            problems.append(f"{url}: a table has neither a thead nor row headers")
    assert not problems, problems


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


def test_the_state_map_never_sends_a_point_it_cannot_plot(tmp_path: Path):
    """A point with no y renders as nothing while the endpoint returns 200.

    The receptor map had 100 points with `y: None`, so the panel would have
    appeared empty on a page that reported success. Points missing either axis
    are left out instead, and a map with none left is a 404 so the panel is
    hidden rather than blank.
    """
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="b2ar", run_id="r",
                          family="gpcr", environment="POPC")
    with store.connect(path) as connection:
        # One complete, one missing the second axis.
        store.put_reference(connection, "b2ar:3SN6", "b2ar", pdb="3SN6",
                            source_state="active", our_state="active",
                            metrics={"tm3_tm6": 15.4, "npxxy_to_inactive": 2.8})
        store.put_reference(connection, "b2ar:2RH1", "b2ar", pdb="2RH1",
                            source_state="inactive", our_state="inactive",
                            metrics={"tm3_tm6": 8.4})

    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    payload = app.test_client().get("/refstates/map/b2ar.json").get_json()

    assert [p["pdb"] for p in payload["points"]] == ["3SN6"]
    assert all(p["y"] is not None for p in payload["points"])
    assert payload["axes"][1], "the receptor map must name its second axis"


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


def test_a_two_hundred_never_carries_arrays_nothing_can_draw(tmp_path: Path):
    """Every series an endpoint returns must be plottable or absent."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _lopsided_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    client = app.test_client()

    payload = client.get("/refstates/run/r/frames.json").get_json()
    assert payload["motifs"], "a 200 with no motifs at all"
    for motif_id, series in payload["motifs"].items():
        length = len(series["time_ns"])
        assert length > 0, f"{motif_id} has metric rows but no time axis"
        assert len(series["s"]) == length, f"{motif_id} s is a different length"
        assert len(series["label"]) == length, f"{motif_id} label length"
        for name, values in series["metrics"].items():
            assert len(values) == length, (
                f"{motif_id} metric {name} is {len(values)} against a "
                f"{length}-point time axis")


def test_a_motif_that_resolved_one_side_still_has_a_time_axis(tmp_path: Path):
    """The specific case: eight metric rows, and before the fix, no points."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _lopsided_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    series = client_series = app.test_client().get(
        "/refstates/run/r/frames.json").get_json()["motifs"]["s:one_sided"]
    assert len(series["time_ns"]) == 8
    assert series["metrics"]["d_inact_a"][0] == pytest.approx(5.0)
    # It has no s score, and that is reported as a gap rather than as a zero.
    assert all(v is None for v in series["s"])


def test_an_endpoint_with_nothing_to_show_returns_404_not_an_empty_body(tmp_path: Path):
    """A panel may be absent. It may not be present and blank."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.create(path)
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    client = app.test_client()
    for url in ("/refstates/run/nothing/frames.json",
                "/refstates/run/nothing/dssp.json",
                "/refstates/map/nothing.json",
                "/refstates/references/nothing",
                "/refstates/run/nothing/viewer.pdb"):
        assert client.get(url).status_code == 404, url


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


@pytest.mark.parametrize("active_s,inactive_s,held,ordered", [
    (0.60, -0.85, True, True),      # both hold their state
    (0.10, -0.85, False, True),     # one drifted into the band
    (-0.40, 0.30, False, False),    # the two crossed
])
def test_the_pair_panel_says_which_of_the_three_things_happened(
        tmp_path: Path, active_s, inactive_s, held, ordered):
    """Each outcome is reported as itself, not rounded towards the hoped-for one."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    _pair_db(path, active_s, inactive_s)
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    page = app.test_client().get("/refstates/").get_data(as_text=True)

    assert "The two start states" in page
    if held:
        assert "stayed on the side they started" in page
    elif ordered:
        assert "drifted into the intermediate band" in page
    else:
        assert "have crossed" in page
    # The separation is always stated, whichever way it went.
    assert f"{active_s - inactive_s:.3f}" in page


def test_no_pair_panel_without_both_halves(tmp_path: Path):
    """One run is not a comparison, and the page must not imply it is."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water",
                          run_fields={"start_state": "DFG-in"})
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    page = app.test_client().get("/refstates/").get_data(as_text=True)
    assert "The two start states" not in page


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
        root / "flexappeal" / "refstates" / "web.py",
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
            for name in (n.strip() for n in names.split(",")):
                if not name or name.startswith("#") or name == "as":
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
