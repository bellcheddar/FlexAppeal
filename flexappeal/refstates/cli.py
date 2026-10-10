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
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
