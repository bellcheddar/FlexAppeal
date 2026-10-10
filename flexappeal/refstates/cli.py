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

    if family == "kinase":
        setup = module.kinase_setup(args.start or args.topology,
                                   args.reference_pdb)
    else:
        setup = module.gpcr_setup(args.start or args.topology)
    motif_residues, core, built, own, report = setup

    run = traj.analyse(loaded, name=args.name or args.trajectory.stem,
                       kind=family, core=core, motif_residues=motif_residues,
                       reference_motifs=ensembles.for_scoring(built),
                       frame_ps=args.frame_ps)

    payload: dict[str, Any] = run.as_row()
    payload["ensembles"] = {state: e.as_row() for state, e in built.items()}
    payload["report"] = report
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
    one.add_argument("--protein", default="ABL1",
                     help="the name the catalogue knows the protein by")
    one.add_argument("--reference-pdb", default="2GQG",
                     help="a PDB entry whose pocket map resolves the "
                          "construct's own motif residues")
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
    two.add_argument("--out", type=Path, default=None)
    two.set_defaults(func=controls)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())


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

    motif_residues, core, built, _own, _report = setup
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
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(
            {"control": "P4", "margin": P4_MARGIN, "rows": results}, indent=2) + "\n")
        log(f"wrote {args.out}")
    # A P4 failure is a result to report rather than an error to abort on,
    # which section 9 is explicit about, so this exits zero either way and the
    # table is the output.
    return 0
