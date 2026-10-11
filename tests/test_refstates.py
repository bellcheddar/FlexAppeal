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
    # The element, not the string: the click handler references the id by name
    # in JavaScript whether or not the panel renders, and a substring check
    # matched that instead of the panel.
    assert 'id="molstar-viewer"' not in page
    assert "<h3>The structure</h3>" not in page
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


def test_a_failed_row_still_renders(tmp_path: Path):
    """A failed row has to show up as failed, not 500 and not vanish."""
    from flexappeal.webapp import create_app

    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC",
                          run_fields={"status": "failed:Particle coordinate is NaN",
                                      "start_state": "active",
                                      "start_model": "4LDE"})
    app = create_app()
    app.config["REFSTATES_DB"] = str(database)
    client = app.test_client()
    for url in ("/refstates/", "/refstates/run/b2ar_apo"):
        response = client.get(url)
        assert response.status_code == 200, f"{url} returned {response.status_code}"
    page = client.get("/refstates/run/b2ar_apo").get_data(as_text=True)
    assert "failed" in page, "the failure is not shown on the run page"


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


def test_the_viewer_route_serves_the_lipid_pair_and_nothing_else():
    """The allowlist is what stops this route serving the database.

    The lipid suffixes are added explicitly rather than by relaxing the match,
    because the protection is the list.
    """
    from flexappeal.refstates import web

    assert "lipids.pdb" in web.VIEWER_SUFFIXES
    assert "lipids.xtc" in web.VIEWER_SUFFIXES
    for forbidden in ("db", "json", "sqlite", "lipids.json", "pdb.db"):
        assert forbidden not in web.VIEWER_SUFFIXES


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


def test_the_thumbnail_route_serves_only_the_two_systems():
    """Section 6.2 asks the cards for a thumbnail; the figures are not static.

    They are generated output, so they live in `data/figures` and need a route,
    and the route uses a fixed set of names for the same reason `viewer_file`
    does: the protection is the list, not the pattern.
    """
    from flexappeal.refstates import web

    assert web.THUMBNAILS == {"abl1", "b2ar"}


def test_a_card_offers_no_thumbnail_when_the_file_is_absent(tmp_path: Path):
    """An <img> whose source 404s is a console error on a real page."""
    from flexappeal.webapp import create_app

    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr", environment="POPC")
    app = create_app()
    app.config["REFSTATES_DB"] = str(database)
    # The real figures directory is beside the app, so this asserts on the
    # route rather than on the card: an unknown system is refused outright.
    client = app.test_client()
    assert client.get("/refstates/thumbnail/nope.png").status_code == 404
    assert client.get("/refstates/thumbnail/b2ar.png").status_code in (200, 404)


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


def test_the_dssp_ribbon_carries_the_reference_consensus(tmp_path: Path):
    """Section 6.2 asks for "active/inactive reference consensus rows".

    The data was already in the analysis payload as a count of DSSP codes per
    state, per motif, per position within that motif. What was missing was
    picking the modal code and resolving the motif position to a residue
    number, which is what puts it on the ribbon's own axis.
    """
    from flexappeal.refstates import web

    payload = {
        "motif_residues": {"TM6": [270, 271, 272]},
        "report": {"sse_profiles": {
            "active": {"TM6": {"0": {"H": 70, "T": 2},
                               "1": {"H": 36, "T": 36},
                               "2": {}}},
            "inactive": {"TM6": {"0": {"T": 20, "H": 8}}},
        }},
    }
    found = web.sse_consensus(payload, [269, 270, 271, 272])

    # Residue 269 is in no motif, so it has no consensus and says so.
    assert found["active"][0] is None
    assert found["active"][1] == {"code": "H", "agreement": 0.972}
    # A residue the references disagree on still reports, with the number that
    # says they disagree: a 50/50 split is a result, not a missing value.
    assert found["active"][2]["agreement"] == 0.5
    # An empty histogram is not a consensus of anything.
    assert found["active"][3] is None
    assert found["inactive"][1] == {"code": "T", "agreement": 0.714}
    # The rows line up with the ribbon's residues, so the browser does not
    # have to join them.
    assert len(found["active"]) == 4


def test_the_dssp_endpoint_serves_the_consensus_beside_the_frames(tmp_path: Path):
    """And it reaches the page, not just the helper."""
    from flexappeal.webapp import create_app

    app = create_app()
    client = app.test_client()
    response = client.get("/refstates/run/b2ar_apo/dssp.json")
    if response.status_code == 404:
        pytest.skip("no analysed b2ar_apo on this machine")
    body = response.get_json()
    assert "consensus" in body, sorted(body)
    assert set(body["consensus"]) == {"active", "inactive"}
    for state, row in body["consensus"].items():
        assert len(row) == len(body["residues"]), state
        assert any(cell for cell in row), f"{state} consensus is entirely empty"


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


def test_every_state_the_strip_can_show_has_a_shape(tmp_path: Path):
    """Section 6.1: never by colour alone.

    The label always carried that, which is why nothing ever broke, but there
    was no `.md-state` rule in the stylesheet at all: the `data-shape`
    attribute the template wrote was read by nothing. The shape follows from
    the state, so it is drawn from the state class and the separate map that
    produced the attribute is gone.
    """
    root = Path(__file__).resolve().parents[1]
    css = (root / "flexappeal" / "static" / "brand.css").read_text()
    for state in ("active", "inactive", "intermediate", "conflict"):
        assert f".md-state-{state}::before" in css, state
    assert ".md-state::before" in css

    # And the fact lives in one place now.
    web_py = (root / "flexappeal" / "refstates" / "web.py").read_text()
    assert "STATE_SHAPES" not in web_py, (
        "the shapes map is back, and nothing reads the attribute it produces")
    run_html = (root / "flexappeal" / "templates"
                / "refstates_run.html").read_text()
    assert "data-shape" not in run_html


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


def test_the_run_view_gives_each_motif_its_residue_endpoints(tmp_path: Path):
    """The ribbon's motif bands need numbers, not the span sentence."""
    from flexappeal.webapp import create_app

    database = tmp_path / "t.db"
    store.ingest_analysis(database, _analysis_payload(), system_id="abl1",
                          run_id="abl1_apo", family="kinase",
                          environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(database)
    page = app.test_client().get("/refstates/run/abl1_apo").get_data(as_text=True)
    assert 'id="motif-meta"' in page
    found = json.loads(page.split('id="motif-meta">')[1].split("</script>")[0])
    assert found, "the motif metadata is empty, so no band can be drawn"
    for entry in found:
        assert entry["first"] <= entry["last"]
        assert entry["token"], f"{entry['name']} has no colour token"


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


def test_a_uniform_ensemble_keeps_the_first_reference_as_its_target():
    """Choosing by modal length must not move a set that is already uniform.

    Every ABL1 and beta2AR ensemble in this build is uniform, so this is what
    makes the change safe for the production rows: same target, same transform,
    same scores.
    """
    found = _ensemble_of({"active": [(f"B{i:03d}", 43) for i in range(5)],
                          "inactive": [(f"C{i:03d}", 43) for i in range(14)]})
    report = ensembles.align_to_common_frame(found)
    assert report["target"] == "B000", "the target moved on a uniform set"
    assert report["n_aligned"] == 19
    assert report["n_dropped"] == 0


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


def test_the_measurement_table_shows_what_moved_and_omits_what_was_not_measured(tmp_path: Path):
    """Driven by a table, so adding a measurement does not mean editing a page.

    And a reading absent from both ends is left out rather than rendered blank,
    so the panel is what was measured rather than a list of what might have
    been.
    """
    import json as _json

    from flexappeal.refstates import web
    from flexappeal.webapp import create_app

    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (analysis / "r.json").write_text(_json.dumps({"frames": [
        {"time_ns": 0.0, "measurements": {
            "alphac_salt_bridge": {"distance_a": 3.40, "state": "in"},
            "p_loop": {"tip_to_hinge_a": 16.0}}},
        {"time_ns": 5.0, "measurements": {
            "alphac_salt_bridge": {"distance_a": 2.60, "state": "in"},
            "p_loop": {"tip_to_hinge_a": 15.0}}},
    ]}))

    app = create_app()
    original = web.analysis_file
    web.analysis_file = lambda run_id: analysis / f"{run_id}.json"
    try:
        with app.app_context():
            rows = web.frame_measurements("r")
    finally:
        web.analysis_file = original

    by_label = {r["label"]: r for r in rows}
    assert "alphaC salt bridge" in by_label
    bridge = by_label["alphaC salt bridge"]
    assert bridge["start"] == pytest.approx(3.40)
    assert bridge["end"] == pytest.approx(2.60)
    assert bridge["change"] == pytest.approx(-0.80)
    assert bridge["state"] == "in"

    # Nothing in this run measured TM6, so it is absent rather than blank.
    assert "TM6 axis rotation" not in by_label
    assert "Ionic lock" not in by_label


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


def test_the_sse_panel_shows_the_margin_and_withholds_an_undecided_lean(tmp_path: Path):
    """A lean on a small margin must not appear as a finding on the page."""
    import json as _json

    from flexappeal.refstates import web
    from flexappeal.webapp import create_app

    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (analysis / "r.json").write_text(_json.dumps({"dssp": {"motifs": {
        "decides": {"helix_fraction_last": 0.9, "reference_agreement": {
            "n_residues": 20, "mean_active_agreement": 0.90,
            "mean_inactive_agreement": 0.70, "margin": 0.20,
            "decisive": True, "leans": "active"}},
        "undecided": {"helix_fraction_last": 0.8, "reference_agreement": {
            "n_residues": 12, "mean_active_agreement": 0.864,
            "mean_inactive_agreement": 0.888, "margin": 0.024,
            "decisive": False, "leans": "neither"}},
        "unscored": {"helix_fraction_last": 0.5, "reference_agreement": {
            "n_residues": 0}},
    }}}))

    app = create_app()
    original = web.analysis_file
    web.analysis_file = lambda run_id: analysis / f"{run_id}.json"
    try:
        with app.app_context():
            rows = web.sse_agreement("r")
    finally:
        web.analysis_file = original

    by_motif = {r["motif"]: r for r in rows}
    assert set(by_motif) == {"decides", "undecided"}, "an unscored motif leaked in"
    assert by_motif["decides"]["decisive"] is True
    assert by_motif["decides"]["leans"] == "active"
    assert by_motif["undecided"]["decisive"] is False
    assert by_motif["undecided"]["margin"] == pytest.approx(0.024)


def test_the_whole_protein_series_reaches_the_plot_with_both_versions(tmp_path: Path):
    """Section 5.6 wants the smoothed line for display and the raw stored.

    Both are served, and the time axis has to be the same length as the series
    or the plot silently drops points off one end.
    """
    import json as _json

    from flexappeal.refstates import web
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    analysis = tmp_path / "analysis"
    analysis.mkdir()
    (analysis / "r.json").write_text(_json.dumps({
        "frames": [{"time_ns": t / 10.0} for t in range(40)],
        "whole_s_raw": [0.1 + t * 0.001 for t in range(40)],
        "whole_s_smoothed": [0.12 + t * 0.001 for t in range(40)],
        "smoothing_window_ns": 2.0,
    }))

    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    original = web.analysis_file
    web.analysis_file = lambda run_id: analysis / f"{run_id}.json"
    try:
        payload = app.test_client().get("/refstates/run/r/frames.json").get_json()
    finally:
        web.analysis_file = original

    whole = payload["whole"]
    assert len(whole["raw"]) == 40
    assert len(whole["smoothed"]) == 40
    assert len(whole["time_ns"]) == len(whole["raw"]), (
        "the time axis and the series must be the same length")
    assert whole["window_ns"] == pytest.approx(2.0)


def test_a_run_with_no_analysis_file_still_serves_its_motifs(tmp_path: Path):
    """The whole series is additive: without it the motif lines still plot."""
    from flexappeal.webapp import create_app

    path = tmp_path / "t.db"
    store.ingest_analysis(path, _analysis_payload(), system_id="s", run_id="r",
                          family="kinase", environment="water")
    app = create_app()
    app.config["REFSTATES_DB"] = str(path)
    payload = app.test_client().get("/refstates/run/r/frames.json").get_json()
    assert payload["motifs"], "the motif series must survive a missing analysis file"
    assert payload["whole"] == {}


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
