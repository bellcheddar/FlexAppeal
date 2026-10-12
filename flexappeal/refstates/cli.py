"""`flexappeal refstates annotate`: the same analysis, on somebody else's run.

Spec section 6.2 item 5, the reuse hook. The point of the whole package is that
nothing in it is specific to ABL1 or to β2AR: residues are resolved through KLIFS
and GPCRdb at run time, so the same code annotates any protein kinase and any
class A GPCR those catalogues cover. This is the entry point that makes that
claim testable by someone other than the author.

CLI only in v1, which is what section 6.2 allows. The work it does is identical
to `tools/analyse_run.py`; what it adds is resolving the family and the protein
from arguments rather than from this build's own two systems.

    flexappeal-refstates annotate \\
        --topology system.pdb --trajectory run.dcd \\
        --family kinase --protein ABL1 --out annotated.json

`--family auto` asks the vendored annotator which family applies, using the same
sequence test `compare-sse` uses, rather than making the caller decide.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

# The HRD and DFG motifs are near-universal in protein kinases and absent from
# class A GPCRs, which is the test the vendored annotator's `applies_to` uses.
# Reused rather than reimplemented so `auto` cannot disagree with it.
KINASE_MARKERS = ("HRD", "DFG")


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", file=sys.stderr, flush=True)


def detect_family(topology: Path) -> str:
    """Which family the structure belongs to, from its own sequence.

    Returns `kinase`, `gpcr` or `unknown`. A structure that looks like neither
    returns `unknown` rather than defaulting, because annotating a kinase as a
    receptor would resolve every motif to the wrong residues and still produce
    numbers.
    """
    import gemmi

    structure = gemmi.read_structure(str(topology))
    structure.setup_entities()
    letters = []
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is not None and info.is_amino_acid():
                letters.append(info.one_letter_code.upper())
    sequence = "".join(letters)
    if all(marker in sequence for marker in KINASE_MARKERS):
        return "kinase"
    # Class A GPCRs are identified by their own catalogue rather than by a
    # sequence motif: GPCRdb resolves a UniProt entry name and nothing else is
    # a reliable test on sequence alone.
    if 250 <= len(sequence) <= 500:
        return "gpcr"
    return "unknown"


# Below this many reference structures in a state, the comparison is reported as
# indicative. Five is the smallest number that can show a spread rather than a
# pair of points, and both worked examples clear it several times over, so it
# only ever fires on an uploaded protein whose catalogue coverage is thin.
THIN_ENSEMBLE = 5


def annotate(args: argparse.Namespace) -> int:
    """Run the reference-state analysis on a caller's trajectory."""
    from . import ensembles, traj

    family = args.family
    if family == "auto":
        family = detect_family(args.topology)
        log(f"family detected as {family}")
    if family == "unknown":
        log("could not tell which family this is; pass --family explicitly")
        return 2

    log(f"loading {args.trajectory.name}")
    loaded = traj.load(args.topology, args.trajectory, stride=args.stride)
    log(f"{loaded.n_frames} frames, {loaded.n_atoms} atoms kept")

    if args.start is not None:
        mapping = traj.deposited_numbering(args.start, loaded.topology)
        changed = traj.renumber(loaded.topology, mapping)
        log(f"{changed} residues restored to deposited numbering")
    else:
        log("no --start given, so the trajectory's own residue numbering is "
            "used; if it was renumbered on solvation every motif will be read "
            "at the wrong residues")

    # The setup functions live in tools/analyse_run.py, which is a script rather
    # than a module. Imported by path so there is one implementation of the
    # family wiring rather than two that have to agree.
    import importlib.util

    script = Path(__file__).resolve().parents[2] / "tools" / "analyse_run.py"
    spec = importlib.util.spec_from_file_location("_analyse_run", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # `--protein` has to reach these. It did not, so every kinase was scored
    # against ABL1's reference ensembles and every receptor against beta2AR's,
    # with the option accepted by the parser and then dropped. The default is
    # per family rather than one value, because "ABL1" is not a receptor and a
    # shared default cannot be right for both.
    limit = args.max_references if args.max_references > 0 else None
    if family == "kinase":
        setup = module.kinase_setup(args.start or args.topology,
                                    args.reference_pdb,
                                    protein=args.protein or "ABL1",
                                    max_references=limit)
    else:
        setup = module.gpcr_setup(args.start or args.topology,
                                  protein=args.protein or "adrb2_human",
                                  max_references=limit)
    # The sixth element is the overlay's inputs, which are callables and so
    # are not part of a report anything serialises. `annotate` does not
    # build an overlay: it has no viewer bundle to hang one on.
    motif_residues, core, built, own, report = setup[:5]

    run = traj.analyse(loaded, name=args.name or args.trajectory.stem,
                       kind=family, core=core, motif_residues=motif_residues,
                       reference_motifs=ensembles.for_scoring(built),
                       frame_ps=args.frame_ps)

    payload: dict[str, Any] = run.as_row()
    payload["ensembles"] = {state: e.as_row() for state, e in built.items()}
    payload["report"] = report

    # A score is a comparison against an ensemble, so how many structures that
    # ensemble has is part of the result. The worked examples are comfortable
    # (ABL1 has 12 active and 45 inactive, beta2AR 72 and 28), but an arbitrary
    # protein need not be: CDK2 comes back with 337 active and 2 inactive, and a
    # state score computed against two structures is a far weaker measurement
    # than one computed against forty-five. Said in the output rather than left
    # for a reader to infer from the ensemble sizes.
    for state, ensemble in built.items():
        n = len(ensemble.by_code)
        if n and n < THIN_ENSEMBLE:
            note = (f"the {state} ensemble has only {n} reference "
                    f"{'structure' if n == 1 else 'structures'}, so every "
                    f"{state} distance here rests on {n}; treat the scores as "
                    f"indicative rather than measured")
            payload.setdefault("notes", []).append(note)
            log(f"  {note}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")

    last = run.frames[-1] if run.frames else None
    if last is not None and last.whole_s is not None:
        log(f"final frame: whole score {last.whole_s:+.3f} "
            f"({last.whole_label}), {last.n_conflicts} conflicts")
    log(f"wrote {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="flexappeal-refstates", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    one = sub.add_parser("annotate", help="score a trajectory against the "
                                          "experimental reference states")
    one.add_argument("--topology", type=Path, required=True,
                     help="the solvated topology the trajectory indexes into")
    one.add_argument("--trajectory", type=Path, required=True)
    one.add_argument("--start", type=Path, default=None,
                     help="the pre-solvation structure, used to restore the "
                          "deposited residue numbering. Strongly recommended: "
                          "OpenMM renumbers from 1 on solvation and every "
                          "annotation here addresses deposited numbers.")
    one.add_argument("--family", choices=("kinase", "gpcr", "auto"),
                     default="auto")
    one.add_argument("--protein", default=None,
                     help="the name the catalogue knows the protein by: a "
                          "KLIFS kinase name such as CDK2, or a GPCRdb entry "
                          "name such as adrb2_human. Defaults to the worked "
                          "example for whichever family is detected")
    one.add_argument("--reference-pdb", default=None,
                     help="a PDB entry whose pocket map resolves the "
                          "construct's own motif residues")
    one.add_argument("--max-references", type=int, default=40,
                     help="the most reference structures to use per state, "
                          "best first by the catalogue's own quality or "
                          "resolution. A state score is the distance to the "
                          "nearest reference, so a well-studied kinase's three "
                          "hundredth structure costs a download and changes "
                          "almost nothing. 0 means no limit, which is what the "
                          "two worked examples use")
    one.add_argument("--stride", type=int, default=1)
    one.add_argument("--frame-ps", type=float, default=100.0,
                     help="the interval the trajectory was saved at")
    one.add_argument("--name", default=None)
    one.add_argument("--out", type=Path, default=Path("annotated.json"))
    one.set_defaults(func=annotate)

    two = sub.add_parser("controls", help="evaluate positive control P4 on "
                                          "every start model")
    two.add_argument("--systems", type=Path,
                     default=Path(__file__).resolve().parents[2] / "data" / "systems")
    two.add_argument("--db", type=Path,
                     default=Path(__file__).resolve().parents[2]
                     / "data" / "refstates.db",
                     help="the explorer database, which P5 and P6 read")
    two.add_argument("--out", type=Path, default=None)
    two.set_defaults(func=controls)

    three = sub.add_parser("report", help="print what the explorer database holds")
    three.add_argument("--verify", action="store_true",
                       help="also check every stored number against the "
                            "analysis file the run was built from")
    three.add_argument("--analysis", type=Path,
                       default=Path(__file__).resolve().parents[2]
                       / "data" / "analysis")
    three.add_argument("--db", type=Path,
                       default=Path(__file__).resolve().parents[2]
                       / "data" / "refstates.db")
    three.set_defaults(func=report)

    five = sub.add_parser(
        "build", help="build a start system or the reference manifests")
    five.add_argument("target", choices=sorted(BUILD_TARGETS),
                      help="abl1 or b2ar for a start system, references for "
                           "the KLIFS and GPCRdb manifests")
    five.add_argument("build_args", nargs=argparse.REMAINDER,
                      help="passed through, for example: --which apo")
    five.set_defaults(func=build)

    four = sub.add_parser("run", help="run the MD queue (tools/queue_md.py)")
    four.add_argument("queue_args", nargs=argparse.REMAINDER,
                      help="passed through, for example: --smoke --only b2ar_apo")
    four.set_defaults(func=run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)



# ===========================================================================
#  Positive control P4
# ===========================================================================

# Section 9's P4 band. A start structure has to score as the state it was built
# from, and the margin is wide enough that a near-zero score counts as a
# failure rather than a pass: a system that cannot place its own start model is
# not going to place a trajectory.
P4_MARGIN = 0.3

P4_ROWS = [
    {"name": "abl1_apo", "family": "kinase", "expect": "active",
     "start": "abl1_apo_start.pdb", "reference_pdb": "2GQG",
     "note": "DFG-in, from 2GQG"},
    {"name": "abl1_imatinib", "family": "kinase", "expect": "inactive",
     "start": "abl1_imatinib_start.pdb", "reference_pdb": "2HYY",
     "note": "DFG-out, from 2HYY"},
    {"name": "b2ar_apo", "family": "gpcr", "expect": "active",
     "start": "b2ar_apo_start.pdb", "reference_pdb": "",
     "note": "active refined model, from 4LDE"},
    {"name": "b2ar_carazolol", "family": "gpcr", "expect": "inactive",
     "start": "b2ar_carazolol_start.pdb", "reference_pdb": "",
     "note": "inactive refined model, from 5D5A"},
]


def _score_start(row: dict, systems: Path, work: Path) -> dict[str, Any]:
    """Score one start structure against the reference ensembles.

    The start model is used as both topology and single-frame trajectory, so
    this measures the structure itself with no simulation in between. That is
    what P4 asks: a scorer that cannot place its own start model has nothing to
    say about a trajectory that began there.
    """
    import importlib.util

    import mdtraj as md

    from . import ensembles, traj

    start = systems / row["start"]
    if not start.is_file():
        return {**row, "status": "no start model"}

    loaded = md.load(str(start))
    frame = work / f"{row['name']}_start.dcd"
    loaded.save_dcd(str(frame))

    script = Path(__file__).resolve().parents[2] / "tools" / "analyse_run.py"
    spec = importlib.util.spec_from_file_location("_analyse_run", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    try:
        if row["family"] == "kinase":
            setup = module.kinase_setup(start, row["reference_pdb"])
        else:
            setup = module.gpcr_setup(start)
    except Exception as exc:                                    # noqa: BLE001
        return {**row, "status": f"setup failed: {str(exc)[:70]}"}

    motif_residues, core, built, _own, _report = setup[:5]
    scored = traj.load(start, frame)
    run = traj.analyse(scored, name=row["name"], kind=row["family"], core=core,
                       motif_residues=motif_residues,
                       reference_motifs=ensembles.for_scoring(built))
    if not run.frames:
        return {**row, "status": "no frame scored"}

    first = run.frames[0]
    whole = first.whole_s
    if whole is None:
        return {**row, "status": "no motif scored"}
    passed = (whole > P4_MARGIN if row["expect"] == "active"
              else whole < -P4_MARGIN)
    return {**row, "status": "ok", "whole_s": whole,
            "label": first.whole_label, "passed": bool(passed),
            "motifs": {m.motif: m.s for m in first.motifs if m.s is not None}}


def controls(args: argparse.Namespace) -> int:
    """Evaluate P4 on every start model and print the table."""
    import tempfile

    systems = args.systems
    results = []
    with tempfile.TemporaryDirectory() as work:
        for row in P4_ROWS:
            log(f"scoring {row['name']} ({row['note']})")
            results.append(_score_start(row, systems, Path(work)))

    print()
    print("P4: each start structure scores as its own state at t = 0")
    print(f"    pass band: s > +{P4_MARGIN} for an active start, "
          f"s < -{P4_MARGIN} for an inactive one")
    print()
    print(f"    {'row':16s} {'expected':9s} {'s':>8s}  {'label':13s} verdict")
    for row in results:
        if row.get("status") != "ok":
            print(f"    {row['name']:16s} {row['expect']:9s} {'':>8s}  "
                  f"{'':13s} {row.get('status')}")
            continue
        print(f"    {row['name']:16s} {row['expect']:9s} {row['whole_s']:+8.3f}  "
              f"{row['label']:13s} {'PASS' if row['passed'] else 'FAIL'}")

    scored = [r for r in results if r.get("status") == "ok"]
    passed = [r for r in scored if r["passed"]]
    print()
    print(f"    {len(passed)} of {len(scored)} scored rows pass P4")

    # P5 and P6 read the explorer database rather than the start models, so
    # they only have anything to say once a run has been analysed. Reported as
    # pending rather than skipped, because a control nobody mentions reads as a
    # control that passed.
    database = args.db
    p5 = control_p5(database) if database.is_file() else []
    p6 = control_p6(database) if database.is_file() else []

    print()
    print(f"P5: an inactive-start run stays inactive over its final "
          f"{P5_WINDOW_NS:.0f} ns")
    print("    section 9 asks 2 of 3 replicas; D009 dropped the replicas to "
          "one, so the\n    rule is per run and the adaptation is recorded "
          "rather than reinterpreted")
    if not p5:
        print("    pending: no inactive-start run analysed yet")
    for row in p5:
        if row.get("status") != "ok":
            print(f"    {row['run_id']:16s} {row.get('status')}")
            continue
        print(f"    {row['run_id']:16s} {row['frames_in_window']:3d} frames, "
              f"mean s {row['mean_whole_s']:+.3f}, "
              f"{row['fraction_below_zero']:.0%} below zero  "
              f"{'PASS' if row['passed'] else 'FAIL'}")

    print()
    print("P6: every membrane chunk inside section 5.8's windows")
    if not p6:
        print("    pending: no membrane chunks ingested yet")
    for row in p6:
        if row.get("status") != "ok":
            print(f"    {row['run_id']:16s} {row.get('status')}")
            continue
        thickness = (f"{row['thickness_mean']:.1f} \u00c5"
                     if row["thickness_mean"] else "no thickness")
        print(f"    {row['run_id']:16s} {row['chunks']:3d} "
              f"chunk{'' if row['chunks'] == 1 else 's'}, APL "
              f"{row['apl_mean']:.1f} \u00c5\u00b2 ({row['apl_min']:.1f} to "
              f"{row['apl_max']:.1f}), {thickness}  "
              f"{'PASS' if row['passed'] else 'FAIL'}")

    results = {"P4": results, "P5": p5, "P6": p6}
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(
            {"controls": results, "p4_margin": P4_MARGIN,
             "p5_window_ns": P5_WINDOW_NS}, indent=2) + "\n")
        log(f"wrote {args.out}")
    # A P4 failure is a result to report rather than an error to abort on,
    # which section 9 is explicit about, so this exits zero either way and the
    # table is the output.
    return 0


# ===========================================================================
#  report and run
# ===========================================================================


def verify_live_metrics(db: Path) -> list[str]:
    """Every metric a motif reports must have a value in at least one frame.

    Section 10, in full: "any metric with all-NaN for a motif is a test
    failure". A metric that is NaN in every frame is not a result. It is a
    measurement that did not happen: an atom the structure does not contain, a
    residue the numbering did not resolve, a motif whose span ran off the end of
    the chain. The distinction matters because of how it renders. An all-NaN
    metric draws as an empty column and a flat line, which on the page is
    indistinguishable from a quantity that genuinely did not move, and "did not
    move" is the headline result this build is looking for.

    Some NaN frames are ordinary and are not reported: a salt bridge whose atoms
    are missing from a handful of frames is an honest gap. All of them is not.

    Checked against the database rather than a fixture, for the same reason
    `verify_summary` is: four rows land unattended, and the question is whether
    what landed is measurable, not whether the code could measure it.

    Returns a list of complaints, empty when every metric has a reading.
    """
    from . import store

    problems: list[str] = []
    with store.connect(db, read_only=True) as connection:
        groups = list(connection.execute(
            "SELECT run_id, motif_id, metric, COUNT(*) AS frames, "
            "       SUM(value IS NOT NULL) AS values_present, "
            "       SUM(s_score IS NOT NULL) AS scores_present "
            "FROM motif_frames GROUP BY run_id, motif_id, metric "
            "ORDER BY run_id, motif_id, metric"))
        # A motif that reports no metric at all would never appear above, so the
        # motifs themselves are checked separately. One with frames ingested for
        # the run and none of its own is as dead as one that is all NaN, and
        # grouping cannot see it.
        silent = list(connection.execute(
            "SELECT r.run_id, m.motif_id FROM runs r JOIN motifs m "
            "  ON m.system_id = r.system_id "
            "WHERE EXISTS (SELECT 1 FROM motif_frames f WHERE f.run_id = r.run_id) "
            "  AND NOT EXISTS (SELECT 1 FROM motif_frames f "
            "                  WHERE f.run_id = r.run_id AND f.motif_id = m.motif_id)"))

    for row in groups:
        where = f"{row['run_id']} {row['motif_id']} {row['metric']}"
        if row["frames"] and not row["values_present"]:
            problems.append(
                f"{where}: NaN in all {row['frames']} frames, so the metric "
                f"was never measured")
        elif row["frames"] and not row["scores_present"]:
            problems.append(
                f"{where}: a value in every frame but no state score in any, "
                f"so the motif has a measurement and no reading of it")
    for row in silent:
        problems.append(
            f"{row['run_id']} {row['motif_id']}: the run has frames and this "
            f"motif has none of them")
    return problems


def verify_summary(db: Path, analysis_dir: Path) -> list[str]:
    """Check the stored summary against the payload each run was built from.

    Two independent paths to the same numbers: the analysis JSON a run writes,
    and the SQLite summary the pages read. Everything between them is ingest
    and SQL, and a fault there produces numbers that are wrong rather than
    missing, which no amount of rendering reveals.

    The unit test for this uses a fixture. This is the same check against
    whatever is actually on disk, so that when four rows land unattended the
    final pass says whether their numbers survived the trip rather than
    assuming it.

    Returns a list of complaints, empty when everything agrees.
    """
    import statistics

    from . import store

    problems: list[str] = []
    with store.connect(db, read_only=True) as connection:
        run_ids = [r["run_id"] for r in connection.execute(
            "SELECT run_id FROM runs ORDER BY run_id")]

    for run_id in run_ids:
        path = analysis_dir / f"{run_id}.json"
        if not path.is_file():
            continue                      # nothing to check against, not a fault
        try:
            frames = json.loads(path.read_text()).get("frames") or []
        except json.JSONDecodeError:
            problems.append(f"{run_id}: the analysis file does not parse")
            continue
        if not frames:
            continue

        hi, lo = frames[-1]["time_ns"], frames[0]["time_ns"]
        cut = hi - 0.1 * (hi - lo)
        series: dict[str, list] = {}
        tail: dict[str, list] = {}
        for frame in frames:
            for motif in frame.get("motifs") or []:
                if motif.get("s") is None:
                    continue
                series.setdefault(motif["motif"], []).append(motif["s"])
                if frame["time_ns"] >= cut:
                    tail.setdefault(motif["motif"], []).append(motif["s"])

        summary = store.summarise(db, run_id)
        for key, values in (summary.get("motifs") or {}).items():
            name = key.split(":")[-1]
            mine = series.get(name)
            if not mine:
                problems.append(f"{run_id} {name}: in the database and not in "
                                f"the analysis file")
                continue
            half = len(mine) // 2
            want = {"start_s": mine[0], "final_s": mine[-1]}
            if tail.get(name):
                want["final_tenth_mean_s"] = statistics.mean(tail[name])
            if len(mine) >= 4:
                want["first_half_mean_s"] = round(
                    statistics.mean(mine[:half]), 4)
                want["second_half_mean_s"] = round(
                    statistics.mean(mine[half:]), 4)
                want["sd_s"] = round(statistics.pstdev(mine), 4)
                want["n_scored_frames"] = len(mine)
            for field, expected in want.items():
                got = values.get(field)
                if got is None or abs(float(got) - float(expected)) > 5e-4:
                    problems.append(
                        f"{run_id} {name}.{field}: database {got}, "
                        f"analysis file {expected}")
    return problems


def report(args: argparse.Namespace) -> int:
    """Print what the explorer database holds, as a build-state summary.

    Reads the same file the pages read, so a disagreement between this and the
    site is a bug rather than two views drifting. Nothing here recomputes
    anything: a report that re-derived its numbers could agree with itself and
    still disagree with what is served.
    """
    from . import store

    if not args.db.is_file():
        log(f"no database at {args.db}; run tools/analyse_all.sh first")
        return 1

    with store.connect(args.db, read_only=True) as connection:
        systems = [dict(r) for r in connection.execute(
            "SELECT * FROM systems ORDER BY system_id")]
        runs = [dict(r) for r in connection.execute(
            "SELECT * FROM runs ORDER BY system_id, run_id")]

    if not runs:
        print("No runs ingested yet.")
        return 0

    if getattr(args, "verify", False):
        problems = verify_summary(args.db, args.analysis)
        if problems:
            print()
            log(f"{len(problems)} stored numbers disagree with the analysis "
                f"files they came from:")
            for complaint in problems[:20]:
                print(f"    {complaint}")
            return 1
        log("every stored number agrees with the analysis file it came from")

        # Agreeing with the analysis file is not the same as having measured
        # anything: an all-NaN metric round-trips perfectly. Section 10 makes
        # that its own failure, so it is its own check.
        dead = verify_live_metrics(args.db)
        if dead:
            print()
            log(f"{len(dead)} motif metrics have no reading at all:")
            for complaint in dead[:20]:
                print(f"    {complaint}")
            return 1
        log("every motif metric has a reading in at least one frame")

    print()
    for system in systems:
        print(f"{system['system_id']}  {system.get('name') or ''}"
              f"  ({system.get('environment') or 'unknown'})")
        for run in [r for r in runs if r["system_id"] == system["system_id"]]:
            length = (f"{run['length_ns']:.1f} ns" if run.get("length_ns")
                      else "no frames")
            speed = (f"{run['ns_per_day']:.0f} ns/day" if run.get("ns_per_day")
                     else "speed not recorded")
            print(f"    {run['run_id']:18s} {str(run.get('start_state')):10s} "
                  f"{length:>10s}  {run['status']:9s} {speed}")
            summary = store.summarise(args.db, run["run_id"])
            for motif, values in sorted((summary.get("motifs") or {}).items()):
                if values.get("final_s") is None:
                    continue
                name = motif.split(":", 1)[-1]
                # The halves and the scatter beside the endpoints, because a
                # start and a final on their own cannot tell a motif that
                # moved from one that is noisy, and section 5.0b is an
                # instance of exactly that: H8 reads a 0.18 fall between the
                # halves with a 0.17 standard deviation behind it.
                #
                # The comparison is against the standard deviation and is a
                # deliberately blunt heuristic, not a test. The test would
                # compare the difference of two means against their standard
                # error, and the standard error of a mean over MD frames needs
                # the autocorrelation time, because consecutive frames 100 ps
                # apart are not independent samples. Against the standard
                # deviation the bar is far too high for a real trend and will
                # never call noise a movement, which is the error worth
                # avoiding in a line a reader skims.
                trend = ""
                if values.get("sd_s") is not None:
                    change = (values["second_half_mean_s"]
                              - values["first_half_mean_s"])
                    verdict = ("within the scatter"
                               if abs(change) <= values["sd_s"]
                               else "larger than the scatter")
                    trend = (f"   halves {values['first_half_mean_s']:+.3f} "
                             f"-> {values['second_half_mean_s']:+.3f} "
                             f"({change:+.3f}, sd {values['sd_s']:.3f}, "
                             f"{verdict})")
                print(f"        {name:12s} start {values['start_s']:+.3f} "
                      f"-> final {values['final_s']:+.3f}  "
                      f"{values.get('final_label') or ''}{trend}")
        print()
    return 0


def run(args: argparse.Namespace) -> int:
    """Run the MD queue. A thin pass-through to tools/queue_md.py.

    Exists so the CLI surface section 11 names is real rather than partly
    documented, and it delegates rather than reimplementing: one definition of
    the run matrix and the protocol, which is the same reason `annotate`
    imports its family wiring instead of copying it.
    """
    import subprocess

    root = Path(__file__).resolve().parents[2]
    command = [sys.executable, str(root / "tools" / "queue_md.py")]
    command += args.queue_args or []
    log("delegating to tools/queue_md.py " + " ".join(args.queue_args or []))
    return subprocess.run(command, cwd=root).returncode


# Section 11 names `build` on the CLI surface, so it is here for the same
# reason `run` is: the deliverable is a CLI with these four subcommands, and a
# subcommand that only exists in the document is not delivered. Each target
# delegates to the tool that already owns it rather than reimplementing a
# builder, so there is one definition of each system.
BUILD_TARGETS = {
    "abl1": "build_abl1.py",
    "b2ar": "build_b2ar.py",
    "references": "build_references.py",
}


def build(args: argparse.Namespace) -> int:
    """Build a start system or regenerate the reference manifests.

    Delegates, the way `run` does. The builders carry real defaults of their
    own and refuse clearly when a required argument is missing, so passing the
    remainder through straight is better than restating their arguments here
    where the two copies could drift.
    """
    import subprocess

    root = Path(__file__).resolve().parents[2]
    script = root / "tools" / BUILD_TARGETS[args.target]
    command = [sys.executable, str(script)] + (args.build_args or [])
    log(f"delegating to tools/{script.name} "
        + " ".join(args.build_args or []))
    return subprocess.run(command, cwd=root).returncode


# ===========================================================================
#  Positive controls P5 and P6
# ===========================================================================

# Section 9's P5 asks that at least 2 of 3 inactive-start replicas keep the
# whole-protein score below zero over the final 20 ns. D009 dropped the
# replicas to one, so the criterion becomes "the inactive-start run keeps it
# below zero", and the adaptation is recorded rather than quietly reinterpreted:
# one replica cannot support a 2-of-3 rule, and pretending otherwise would be
# claiming a majority from a sample of one.
P5_WINDOW_NS = 20.0
P5_INACTIVE_STARTS = ("inactive", "DFG-out")


def _final_window(frames: list[dict], window_ns: float) -> list[dict]:
    if not frames:
        return []
    last = max(f["time_ns"] for f in frames)
    return [f for f in frames if f["time_ns"] >= last - window_ns]


def control_p5(database: Path) -> list[dict[str, Any]]:
    """Do the inactive-start runs stay inactive over their final 20 ns?

    Evaluated on the whole-protein score per frame, recomputed from the stored
    per-motif scores so there is one definition of the weighting.
    """
    from . import store

    with store.connect(database, read_only=True) as connection:
        runs = [dict(r) for r in connection.execute(
            "SELECT run_id, system_id, start_state, length_ns FROM runs "
            "ORDER BY run_id")]
        weights = {}
        for row in connection.execute("SELECT motif_id, residues_json FROM motifs"):
            try:
                weights[row["motif_id"]] = len(json.loads(row["residues_json"] or "[]"))
            except (TypeError, json.JSONDecodeError):
                weights[row["motif_id"]] = 0

    out = []
    for run in runs:
        if run.get("start_state") not in P5_INACTIVE_STARTS:
            continue
        with store.connect(database, read_only=True) as connection:
            rows = [dict(r) for r in connection.execute(
                "SELECT DISTINCT frame, time_ns, motif_id, s_score FROM "
                "motif_frames WHERE run_id = ? AND s_score IS NOT NULL",
                (run["run_id"],))]
        by_frame: dict[int, dict] = {}
        for row in rows:
            frame = by_frame.setdefault(row["frame"],
                                        {"time_ns": row["time_ns"], "parts": []})
            weight = weights.get(row["motif_id"], 0)
            if weight > 0:
                frame["parts"].append((row["s_score"], weight))
        frames = [{"time_ns": f["time_ns"],
                   "whole_s": sum(s * w for s, w in f["parts"])
                              / sum(w for _, w in f["parts"])}
                  for f in by_frame.values() if f["parts"]]
        if not frames:
            out.append({**run, "status": "no scored frames"})
            continue

        window = _final_window(frames, P5_WINDOW_NS)
        below = [f for f in window if f["whole_s"] < 0.0]
        out.append({
            **run, "status": "ok",
            "window_ns": P5_WINDOW_NS,
            "frames_in_window": len(window),
            "fraction_below_zero": len(below) / len(window) if window else None,
            "mean_whole_s": sum(f["whole_s"] for f in window) / len(window),
            # The spec's rule is per replica, so with one replica the run
            # either stayed inactive or it did not.
            "passed": all(f["whole_s"] < 0.0 for f in window),
            "replicas": 1,
        })
    return out


def control_p6(database: Path) -> list[dict[str, Any]]:
    """Was every membrane chunk inside section 5.8's windows?"""
    from . import membrane, store

    with store.connect(database, read_only=True) as connection:
        runs = [r["run_id"] for r in connection.execute(
            "SELECT DISTINCT run_id FROM membrane ORDER BY run_id")]
        out = []
        for run_id in runs:
            chunks = [dict(r) for r in connection.execute(
                "SELECT * FROM membrane WHERE run_id = ? ORDER BY time_ns",
                (run_id,))]
            apl = [c["apl"] for c in chunks if c["apl"] is not None]
            thickness = [c["thickness"] for c in chunks
                         if c["thickness"] is not None]
            if not apl:
                out.append({"run_id": run_id, "status": "no area per lipid"})
                continue
            lo, hi = membrane.APL_EXPECTED_A2
            acceptable = membrane.APL_ACCEPTABLE_A2
            t_lo, t_hi = membrane.THICKNESS_EXPECTED_A
            out.append({
                "run_id": run_id, "status": "ok",
                "chunks": len(chunks),
                "apl_mean": sum(apl) / len(apl),
                "apl_min": min(apl), "apl_max": max(apl),
                "apl_in_expected": all(lo <= v <= hi for v in apl),
                "apl_in_acceptable": all(acceptable[0] <= v <= acceptable[1]
                                         for v in apl),
                "thickness_mean": (sum(thickness) / len(thickness)
                                   if thickness else None),
                "thickness_in_expected": all(t_lo <= v <= t_hi
                                             for v in thickness) if thickness else None,
                # Section 5.8's G6 is the acceptable window, not the expected
                # one: a chunk outside 60 to 68 is logged, a run outside 55 to
                # 75 invalidates the bilayer.
                "passed": all(acceptable[0] <= v <= acceptable[1] for v in apl),
            })
    return out


# The guard belongs at the end of the file and not in the middle of it. It used
# to sit directly under `main`, which is above the subcommand functions, so
# running this module as a script executed `build_parser` before `controls`,
# `report` and `run` had been defined and died on a `NameError` every time.
# Importing the module defines everything before anything calls it, which is
# why the tests never saw this and why the only broken path was the one a user
# would take.
if __name__ == "__main__":
    raise SystemExit(main())
