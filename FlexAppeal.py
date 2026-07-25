#!/usr/bin/env python3
"""FlexAppeal -- prepare OpenMM molecular dynamics runs, analyse what comes back.

    ./FlexAppeal.py serve                     run the web app locally
    ./FlexAppeal.py inspect 1aki.pdb          report what is in a structure
    ./FlexAppeal.py fetch 1AKI -o 1aki.pdb    retrieve a structure
    ./FlexAppeal.py validate config.json      check a configuration
    ./FlexAppeal.py docs                      regenerate docs/options.md
    ./FlexAppeal.py sweep                     remove expired scratch sessions

The web app is the primary interface; this CLI exists so every stage can be
driven and tested without a browser, and so the option registry can generate its
own documentation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from flexappeal import options as opts  # noqa: E402
from flexappeal import schema, sources, structure  # noqa: E402

# Brand palette as truecolor ANSI, matching install.sh and the generated
# bundle's console output so the whole toolchain looks like one thing.
BLUE, GREEN, AMBER, RED, DIM, BOLD, RESET = (
    "\033[38;2;30;115;190m", "\033[38;2;0;208;132m", "\033[38;2;252;185;0m",
    "\033[38;2;214;54;56m", "\033[2m", "\033[1m", "\033[0m",
)


def _info(msg): print(f"{BLUE}ℹ{RESET} {msg}")
def _ok(msg): print(f"{GREEN}✓{RESET} {msg}")
def _warn(msg): print(f"{AMBER}⚠{RESET} {msg}")
def _err(msg): print(f"{RED}✗{RESET} {msg}", file=sys.stderr)
def _step(msg): print(f"{BLUE}→{RESET} {msg}")


# ===========================================================================


def cmd_serve(args: argparse.Namespace) -> int:
    from flexappeal.webapp import create_app

    app = create_app()
    _info(f"FlexAppeal v{opts.FLEXAPPEAL_VERSION} on http://{args.host}:{args.port}")
    _info("simulations run on the user's machine -- this server only prepares and analyses")
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


def cmd_inspect(args: argparse.Namespace) -> int:
    data = Path(args.path).read_bytes()
    report = structure.analyse(data, Path(args.path).name)

    if args.json:
        print(json.dumps(report.to_dict(), indent=2))
        return 0

    print(f"\n{BOLD}{report.title or report.name}{RESET}")
    meta = [report.method] if report.method else []
    if report.resolution:
        meta.append(f"{report.resolution} Å")
    meta.append(f"{report.models} model(s)")
    print(f"{DIM}{' · '.join(meta)}{RESET}\n")

    print(f"{BOLD}Chains{RESET}")
    for c in report.chains:
        if c.kind == "water":
            continue
        missing = f"  {AMBER}{c.missing_residues} missing{RESET}" if c.missing_residues else ""
        print(f"  {c.id:3s} {c.kind:8s} {c.observed_residues:5d} residues{missing}")
        for gap in c.gaps:
            if gap["after"] is not None:
                print(f"      {DIM}gap of {gap['length']} after residue {gap['after']}{RESET}")

    if report.disulfides:
        print(f"\n{BOLD}Disulfides{RESET}")
        for d in report.disulfides:
            print(f"  {d.chain_a}{d.resid_a}–{d.chain_b}{d.resid_b}  {DIM}{d.distance} Å{RESET}")

    hetero = [h for h in report.heteroatoms if h.category != "water"]
    if hetero:
        print(f"\n{BOLD}Heteroatoms{RESET}")
        for h in hetero:
            desc = f"  {DIM}{h.description}{RESET}" if h.description else ""
            print(f"  {h.name:5s} ×{h.count:<4d} {h.category:10s}{desc}")

    cfg = opts.defaults()
    est = structure.estimate_system_size(report, cfg)
    cfg["_estimated_atoms"] = est["total_atoms"]
    cfg["_solute_atoms"] = report.solute_atoms
    derived = schema.derive(cfg)
    wall = schema.estimate_wall_time(cfg)

    print(f"\n{BOLD}With the default settings{RESET}")
    print(f"  {est['total_atoms']:,} atoms  {DIM}({est['water_molecules']:,} waters, {est['basis']}){RESET}")
    print(f"  {derived['total_ns']:.0f} ns  →  {derived['traj_frames']:,} frames, {derived['traj_size_human']}")
    print(f"  roughly {wall['human']} {DIM}at an estimated {wall['ns_per_day']} ns/day{RESET}")

    if report.warnings:
        print()
        for w in report.warnings:
            _warn(w)
    print()
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    cfg = {"input_source": args.source, "pdb_id": args.accession,
           "uniprot_id": args.accession, "assembly": args.assembly}
    try:
        result = sources.fetch(cfg)
    except sources.SourceError as exc:
        _err(str(exc))
        return 1
    out = Path(args.output or result.filename)
    out.write_bytes(result.data)
    _ok(f"{result.citation} → {out} ({len(result.data):,} bytes)")
    if warning := sources.SOURCE_WARNINGS.get(result.source):
        _warn(warning)
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    raw = json.loads(Path(args.path).read_text())
    if isinstance(raw, dict) and "config" in raw:
        raw = raw["config"]
    if not isinstance(raw, dict):
        _err("expected a JSON object of option ids to values")
        return 1

    result = schema.validate(raw, strict_unknown=args.strict)
    for issue in result.errors:
        _err(f"{issue.option_id or 'config'}: {issue.message}")
    for issue in result.warnings:
        _warn(f"{issue.option_id or 'config'}: {issue.message}")

    if not result.ok:
        return 1

    derived = schema.derive(result.config)
    _ok(f"valid — {derived['total_ns']:.0f} ns, {derived['total_steps']:,} steps, "
        f"{derived['traj_frames']:,} frames, {derived['traj_size_human']}")
    return 0


def cmd_bundle(args: argparse.Namespace) -> int:
    """Build a run bundle from a structure and an optional config, without a browser."""
    from flexappeal import bundle as bundle_mod

    data = Path(args.structure).read_bytes()
    report = structure.analyse(data, Path(args.structure).name)

    cfg = opts.defaults()
    if args.config:
        raw = json.loads(Path(args.config).read_text())
        cfg.update(raw.get("config", raw))
    if args.job_name:
        cfg["job_name"] = args.job_name
    elif cfg["job_name"] == "flexappeal_run":
        cfg["job_name"] = Path(args.structure).stem.replace(".", "_")[:64]

    # The size and wall-time estimates in the generated README come from these.
    estimate = structure.estimate_system_size(report, cfg)
    cfg["_estimated_atoms"] = estimate["total_atoms"]
    cfg["_solute_atoms"] = report.solute_atoms
    cfg.setdefault("chains", [c.id for c in report.chains if c.kind != "water"])

    result = schema.validate(cfg)
    for issue in result.warnings:
        _warn(f"{issue.option_id or 'config'}: {issue.message}")
    if not result.ok:
        for issue in result.errors:
            _err(f"{issue.option_id or 'config'}: {issue.message}")
        return 1

    try:
        built = bundle_mod.build(cfg, data, Path(args.structure).name)
    except bundle_mod.BundleError as exc:
        _err(str(exc))
        return 1

    out = Path(args.output or built.filename)
    out.write_bytes(built.content)
    out.chmod(0o755)

    derived = schema.derive(cfg)
    wall = schema.estimate_wall_time(cfg)
    _ok(f"wrote {out} ({built.size_human})")
    _info(f"{estimate['total_atoms']:,} atoms · {derived['total_ns']:.0f} ns · "
          f"{derived['traj_frames']:,} frames · roughly {wall['human']}")
    _info(f"run it with:  ./{out.name}")
    return 0


def cmd_unpack(args: argparse.Namespace) -> int:
    """Inspect a bundle's contents without running it."""
    from flexappeal import bundle as bundle_mod

    content = Path(args.bundle).read_bytes()
    try:
        files = bundle_mod.unpack(content)
    except bundle_mod.BundleError as exc:
        _err(str(exc))
        return 1

    if args.output:
        target = Path(args.output)
        for name, payload in files.items():
            path = target / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        _ok(f"extracted {len(files)} file(s) to {target}")
    else:
        for name in sorted(files):
            print(f"  {name:20s} {len(files[name]):>10,d} bytes")
    return 0


def cmd_docs(args: argparse.Namespace) -> int:
    """Generate docs/options.md from the registry, so it cannot go stale."""
    lines = [
        "# Every option FlexAppeal exposes",
        "",
        "Generated from `flexappeal/options.py` by `./FlexAppeal.py docs`.",
        "Do not edit this file by hand; edit the registry and regenerate.",
        "",
        f"{len(opts.OPTIONS)} options across {len(opts.GROUPS)} groups, "
        f"targeting OpenMM {opts.OPENMM_VERSION}.",
        "",
    ]
    for group in opts.GROUPS:
        lines += [f"## {group.icon} {group.title}", "", group.blurb, "",
                  "| Option | Type | Default | Units | Applies when | Description |",
                  "|---|---|---|---|---|---|"]
        for opt in opts.BY_GROUP[group.id]:
            default = opt.default
            if isinstance(default, list):
                default = ", ".join(str(v) for v in default) or "—"
            elif isinstance(default, bool):
                default = "on" if default else "off"
            elif default in (None, ""):
                default = "—"
            requires = f"`{opt.requires}`" if opt.requires else "always"
            help_text = opt.help.replace("|", "\\|")
            if opt.experimental:
                help_text = "**Experimental.** " + help_text
            lines.append(
                f"| `{opt.id}`{' ⚙️' if opt.advanced else ''} | {opt.widget} | "
                f"{default} | {opt.units or '—'} | {requires} | {help_text} |"
            )
            if opt.choices:
                for c in opt.choices:
                    # Escaped outside the f-string: Python 3.11 rejects a
                    # backslash inside an f-string expression, and the droplet
                    # and the bundle both target 3.11.
                    choice_help = c.help.replace("|", "\\|")
                    choice_when = f"`{c.requires}`" if c.requires else ""
                    lines.append(
                        f"| {'&nbsp;' * 4}↳ `{c.value}` | choice | | | "
                        f"{choice_when} | *{c.label}* — {choice_help} |"
                    )
        lines.append("")

    lines += ["## Legend", "",
              "- ⚙️ collapsed behind **Advanced** in the form",
              "- *Applies when* is the `requires` predicate: the field is hidden, "
              "skipped by the validator and omitted from the generated script "
              "whenever it evaluates false.", ""]

    out = SCRIPT_DIR / "docs" / "options.md"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(lines))
    _ok(f"wrote {out} ({len(opts.OPTIONS)} options, {len(lines)} lines)")
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    from flexappeal.webapp import DEFAULT_SCRATCH_ROOT, sweep_scratch

    root = Path(args.root) if args.root else DEFAULT_SCRATCH_ROOT
    removed = sweep_scratch(root, ttl=args.ttl)
    _ok(f"removed {removed} expired session(s) from {root}")
    return 0


# ===========================================================================


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="FlexAppeal.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("serve", help="run the web app locally")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8004,
                   help="8004 by default, matching the droplet's port allocation "
                        "(8000 AlphaFraud, 8001 chem_sage, 8002 chatPDB, 8003 BoltzMaker)")
    p.add_argument("--debug", action="store_true")
    p.set_defaults(func=cmd_serve)

    p = sub.add_parser("inspect", help="report what is in a structure file")
    p.add_argument("path")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("fetch", help="retrieve a structure from a public database")
    p.add_argument("accession", help="PDB ID or UniProt accession")
    p.add_argument("--source", choices=["rcsb", "opm", "alphafold"], default="rcsb")
    p.add_argument("--assembly", choices=["asymmetric", "biological"], default="asymmetric",
                   help="the deposited asymmetric unit is often not the biological "
                        "molecule; simulating one chain of a physiological dimer is "
                        "simulating the wrong thing")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("validate", help="check a configuration against the registry")
    p.add_argument("path")
    p.add_argument("--strict", action="store_true",
                   help="treat unknown options as errors rather than warnings")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("bundle", help="build a run bundle without using the web app")
    p.add_argument("structure", help="a PDB or mmCIF file")
    p.add_argument("-c", "--config", help="a config.json from a previous bundle")
    p.add_argument("-n", "--job-name")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_bundle)

    p = sub.add_parser("unpack", help="list or extract a bundle's contents")
    p.add_argument("bundle")
    p.add_argument("-o", "--output", help="extract here instead of listing")
    p.set_defaults(func=cmd_unpack)

    p = sub.add_parser("docs", help="regenerate docs/options.md from the registry")
    p.set_defaults(func=cmd_docs)

    p = sub.add_parser("sweep", help="remove expired scratch sessions")
    p.add_argument("--root")
    p.add_argument("--ttl", type=int, default=4 * 60 * 60, help="seconds")
    p.set_defaults(func=cmd_sweep)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    try:
        return args.func(args)
    except (structure.StructureError, sources.SourceError) as exc:
        _err(str(exc))
        return 1
    except KeyboardInterrupt:
        print()
        return 130


if __name__ == "__main__":
    sys.exit(main())
