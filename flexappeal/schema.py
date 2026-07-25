"""Validation, coercion and derived quantities for a FlexAppeal config.

Three jobs, in the order they happen to a submitted form:

1. ``normalise``  -- fill in values that are implied rather than typed
2. ``validate``   -- coerce types, check ranges and choices, apply physics rules
3. ``derive``     -- compute the numbers the UI and the run script both need
                     (step counts, frame counts, size and time estimates)

Two layers of strictness, following BoltzMaker's parser convention: an unknown
key is a *warning* and is dropped (the user's typo does not stop their run), but
a value that is out of range, not an offered choice, or physically inconsistent
with another value is an *error* and stops the build. The distinction matters
because a silently accepted bad value here becomes a wrong trajectory eight
hours later on someone else's laptop.

``evaluate_predicate`` is the shared conditional engine. The same predicate
string decides whether the browser shows a field, whether the validator demands
it, and whether the script template emits the corresponding lines -- so those
three can never disagree.
"""

from __future__ import annotations

import ast
import math
import re
from dataclasses import dataclass, field
from typing import Any

from . import options as opts
from .options import Option

# ===========================================================================
#  Predicate evaluation
# ===========================================================================

# Only these node types are ever evaluated. Anything else in a predicate is a
# programming error in options.py, not user input -- predicates are never
# user-supplied -- but the whitelist keeps that guarantee explicit and cheap to
# verify rather than resting on the claim.
_ALLOWED_NODES = (
    ast.Expression, ast.BoolOp, ast.UnaryOp, ast.Compare, ast.Name, ast.Load,
    ast.Constant, ast.Tuple, ast.List, ast.And, ast.Or, ast.Not,
    ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE, ast.In, ast.NotIn,
)

_predicate_cache: dict[str, ast.Expression] = {}


class PredicateError(ValueError):
    """A predicate in options.py is malformed or references an unknown option."""


def _parse_predicate(expr: str) -> ast.Expression:
    cached = _predicate_cache.get(expr)
    if cached is not None:
        return cached
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        raise PredicateError(f"cannot parse predicate {expr!r}: {exc}") from None
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise PredicateError(
                f"predicate {expr!r} uses unsupported syntax {type(node).__name__}; "
                f"predicates may only combine option ids, constants, comparisons "
                f"and and/or/not"
            )
        if isinstance(node, ast.Name) and node.id not in opts.BY_ID:
            raise PredicateError(
                f"predicate {expr!r} references unknown option {node.id!r}"
            )
    _predicate_cache[expr] = tree
    return tree


def _eval_node(node: ast.AST, config: dict[str, Any]) -> Any:
    if isinstance(node, ast.Expression):
        return _eval_node(node.body, config)
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        return config.get(node.id, opts.BY_ID[node.id].default)
    if isinstance(node, (ast.Tuple, ast.List)):
        return [_eval_node(e, config) for e in node.elts]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return not _eval_node(node.operand, config)
    if isinstance(node, ast.BoolOp):
        values = (_eval_node(v, config) for v in node.values)
        if isinstance(node.op, ast.And):
            return all(values)
        return any(values)
    if isinstance(node, ast.Compare):
        left = _eval_node(node.left, config)
        for op, comparator in zip(node.ops, node.comparators):
            right = _eval_node(comparator, config)
            if isinstance(op, ast.Eq):
                result = left == right
            elif isinstance(op, ast.NotEq):
                result = left != right
            elif isinstance(op, ast.Lt):
                result = left < right
            elif isinstance(op, ast.LtE):
                result = left <= right
            elif isinstance(op, ast.Gt):
                result = left > right
            elif isinstance(op, ast.GtE):
                result = left >= right
            elif isinstance(op, ast.In):
                result = left in right
            else:  # NotIn
                result = left not in right
            if not result:
                return False
            left = right
        return True
    raise PredicateError(f"unsupported node {type(node).__name__}")


def evaluate_predicate(expr: str | None, config: dict[str, Any]) -> bool:
    """Evaluate a ``requires`` predicate against a config. ``None`` means always."""
    if not expr:
        return True
    return bool(_eval_node(_parse_predicate(expr), config))


def is_active(opt: Option, config: dict[str, Any]) -> bool:
    """Whether an option applies given the rest of the config."""
    return evaluate_predicate(opt.requires, config)


def active_choices(opt: Option, config: dict[str, Any]) -> tuple[opts.Choice, ...]:
    """The choices offered for an option given the rest of the config."""
    return tuple(c for c in opt.choices if evaluate_predicate(c.requires, config))


def active_options(config: dict[str, Any]) -> tuple[Option, ...]:
    return tuple(o for o in opts.OPTIONS if is_active(o, config))


# ===========================================================================
#  Issues
# ===========================================================================


@dataclass(frozen=True)
class Issue:
    level: str  # "error" | "warning"
    message: str
    option_id: str | None = None

    def __str__(self) -> str:
        where = f"{self.option_id}: " if self.option_id else ""
        return f"{self.level}: {where}{self.message}"


@dataclass
class ValidationResult:
    config: dict[str, Any]
    issues: list[Issue] = field(default_factory=list)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "error"]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors


# ===========================================================================
#  Coercion
# ===========================================================================

_TRUTHY = {"1", "true", "yes", "on", "checked"}
_JOB_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")


def _coerce(opt: Option, raw: Any) -> Any:
    """Turn one raw form value into the option's declared type."""
    if opt.widget == "checkbox":
        if isinstance(raw, bool):
            return raw
        return str(raw).strip().lower() in _TRUTHY

    if opt.widget == "multiselect":
        if raw is None:
            return []
        if isinstance(raw, str):
            raw = [p.strip() for p in raw.split(",") if p.strip()]
        return [str(v).strip() for v in raw]

    if opt.widget == "int":
        if raw in (None, ""):
            return opt.default
        # Accept "10.0" from a browser number input that decided to add a decimal.
        return int(float(str(raw).strip()))

    if opt.widget == "number":
        if raw in (None, ""):
            return opt.default
        return float(str(raw).strip())

    if opt.widget == "file":
        return raw

    return str(raw).strip() if raw is not None else ""


# ===========================================================================
#  Normalisation -- values implied by other values
# ===========================================================================


def normalise(config: dict[str, Any]) -> dict[str, Any]:
    """Fill in values the user should not have to set by hand.

    Done before validation so that the rules below see a coherent config, and
    so the form can round-trip: a user who keeps a heteroatom should not also
    have to tick "system contains ligands".
    """
    cfg = dict(config)

    # Keeping any heteroatom means the system needs small-molecule parameters.
    if cfg.get("keep_heteroatoms"):
        cfg["has_ligands"] = True

    # ff19SB was fitted with OPC and is not meant to be used with anything else;
    # CHARMM36 requires its own LJ-on-hydrogen water. Rather than let the user
    # pick an invalid pair and then reject it, set the only valid partner.
    ff = cfg.get("protein_ff")
    if ff == "amber/protein.ff19SB.xml" and cfg.get("solvent_mode") == "explicit":
        cfg["water_model"] = "opc"
    elif ff == "charmm36.xml" and cfg.get("solvent_mode") == "explicit":
        cfg["water_model"] = "charmm_tip3p"
    elif cfg.get("water_model") == "charmm_tip3p" and ff != "charmm36.xml":
        cfg["water_model"] = "tip3p"

    # A membrane system must use the membrane barostat; an isotropic one would
    # squeeze the bilayer plane and the normal together, which is not physical.
    if cfg.get("use_membrane") and cfg.get("barostat") not in ("MonteCarloMembrane", "none"):
        cfg["barostat"] = "MonteCarloMembrane"
    elif not cfg.get("use_membrane") and cfg.get("barostat") == "MonteCarloMembrane":
        cfg["barostat"] = "MonteCarlo"

    # There is no box in implicit solvent or vacuum, so there is nothing for a
    # barostat to scale.
    if cfg.get("solvent_mode") in ("implicit", "vacuum"):
        cfg["barostat"] = "none"
        if cfg.get("nonbonded_method") in ("PME", "LJPME", "Ewald"):
            cfg["nonbonded_method"] = (
                "CutoffNonPeriodic" if cfg.get("solvent_mode") == "implicit" else "NoCutoff"
            )

    return cfg


# ===========================================================================
#  Physics rules -- the checks a single `requires` cannot express
# ===========================================================================


def _rule_timestep(cfg: dict[str, Any]) -> list[Issue]:
    """The timestep has to be justified by what the constraints actually remove."""
    dt = float(cfg.get("timestep") or 0.0)
    constraints = cfg.get("constraints")
    hmr = bool(cfg.get("use_hmr"))
    issues: list[Issue] = []

    if constraints == "None":
        limit = 1.0
    elif constraints == "HBonds":
        limit = 5.0 if hmr else 2.5
    else:  # AllBonds / HAngles
        limit = 5.0 if hmr else 3.5

    if dt > limit:
        issues.append(Issue(
            "error",
            f"a {dt:g} fs timestep is not stable with constraints={constraints} and "
            f"hydrogen mass repartitioning {'on' if hmr else 'off'}; the ceiling here "
            f"is about {limit:g} fs. Either shorten the step, constrain bonds to "
            f"hydrogen, or enable hydrogen mass repartitioning.",
            "timestep",
        ))
    if hmr and constraints == "None":
        issues.append(Issue(
            "error",
            "hydrogen mass repartitioning does nothing useful without constraints on "
            "bonds to hydrogen -- the bond vibration it is meant to slow is still the "
            "fastest motion in the system.",
            "use_hmr",
        ))
    if dt >= 4.0 and not hmr:
        issues.append(Issue(
            "warning",
            f"{dt:g} fs without hydrogen mass repartitioning is aggressive; expect "
            f"visible drift in total energy. Turning HMR on is nearly free.",
            "timestep",
        ))
    return issues


def _rule_box(cfg: dict[str, Any]) -> list[Issue]:
    """Padding must exceed the cutoff or a molecule can see its own image."""
    if cfg.get("solvent_mode") != "explicit" or cfg.get("use_membrane"):
        return []
    if cfg.get("box_sizing") != "padding":
        return []
    padding = float(cfg.get("padding") or 0.0)
    cutoff = float(cfg.get("nonbonded_cutoff") or 0.0)
    if padding < cutoff:
        return [Issue(
            "error",
            f"solvent padding ({padding:g} nm) is smaller than the non-bonded cutoff "
            f"({cutoff:g} nm). The solute would interact with its own periodic image, "
            f"which is a physically meaningless artefact. Use at least {cutoff:g} nm.",
            "padding",
        )]
    if padding < cutoff + 0.2:
        return [Issue(
            "warning",
            f"solvent padding ({padding:g} nm) barely clears the cutoff ({cutoff:g} nm). "
            f"The box shrinks under the barostat, so leave some headroom -- "
            f"{cutoff + 0.3:g} nm would be safer.",
            "padding",
        )]
    return []


def _rule_switching(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    if cfg.get("use_switching"):
        if float(cfg.get("switch_distance") or 0) >= float(cfg.get("nonbonded_cutoff") or 0):
            issues.append(Issue(
                "error",
                "the switching distance must be shorter than the cutoff -- the taper "
                "has to finish before the interaction is truncated.",
                "switch_distance",
            ))
        if cfg.get("dispersion_correction"):
            issues.append(Issue(
                "warning",
                "using a switching function and the analytic dispersion correction "
                "together double-counts the long-range Lennard-Jones energy. CHARMM "
                "force fields expect the switch and no correction.",
                "dispersion_correction",
            ))
    return issues


def _rule_forcefield(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    ff = cfg.get("protein_ff")
    if ff == "charmm36.xml" and float(cfg.get("nonbonded_cutoff") or 0) < 1.15:
        issues.append(Issue(
            "warning",
            "CHARMM36 was parameterised with a 1.2 nm cutoff and a force switch from "
            "1.0 nm. Using the AMBER-style 1.0 nm cutoff with it changes the balance "
            "the force field was fitted to.",
            "nonbonded_cutoff",
        ))
    if ff == "amoeba2018.xml":
        issues.append(Issue(
            "warning",
            "AMOEBA is polarisable and roughly two orders of magnitude slower than a "
            "fixed-charge force field. A 100 ns run that takes a day with ff14SB will "
            "not finish this year.",
            "protein_ff",
        ))
    if cfg.get("use_membrane") and ff != "charmm36.xml":
        issues.append(Issue(
            "warning",
            "the membrane builder uses CHARMM36 lipid parameters. Mixing them with an "
            "AMBER protein force field works mechanically but is not a combination "
            "either force field was validated for; CHARMM36m is the consistent choice.",
            "protein_ff",
        ))
    return issues


def _rule_solvation(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    if cfg.get("solvent_mode") == "vacuum":
        issues.append(Issue(
            "warning",
            "vacuum simulation of a protein is almost always a mistake: without solvent "
            "screening, surface salt bridges collapse and the structure compacts within "
            "a nanosecond.",
            "solvent_mode",
        ))
    if cfg.get("solvent_mode") == "explicit" and not cfg.get("neutralize"):
        issues.append(Issue(
            "warning",
            "a net-charged system under PME is neutralised by an implicit uniform "
            "background charge, which distorts electrostatics near the solute.",
            "neutralize",
        ))
    return issues


def _rule_protocol(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    if not cfg.get("minimize") and cfg.get("solvent_mode") == "explicit":
        issues.append(Issue(
            "error",
            "skipping minimisation on a freshly solvated system will almost certainly "
            "blow the run up in the first few steps -- addSolvent leaves clashes that "
            "dynamics cannot resolve.",
            "minimize",
        ))
    if float(cfg.get("equilibration_duration") or 0) <= 0 and cfg.get("solvent_mode") == "explicit":
        issues.append(Issue(
            "warning",
            "no equilibration means production starts before the density and solvation "
            "shell have settled; the first part of your trajectory will not be at the "
            "state point you asked for.",
            "equilibration_duration",
        ))
    if cfg.get("use_membrane") and float(cfg.get("equilibration_duration") or 0) < 5000:
        issues.append(Issue(
            "warning",
            "membrane systems normally need at least 5 ns of equilibration for the "
            "lipids to pack around the protein and the area per lipid to plateau.",
            "equilibration_duration",
        ))
    return issues


def _rule_output(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    if cfg.get("traj_selection") == "custom" and not str(cfg.get("traj_custom_selection") or "").strip():
        issues.append(Issue(
            "error",
            "a custom atom selection was chosen but no selection expression was given.",
            "traj_custom_selection",
        ))
    if float(cfg.get("traj_interval") or 0) <= 0:
        issues.append(Issue("error", "the trajectory interval must be greater than zero.", "traj_interval"))
    derived = derive(cfg)
    if derived["traj_frames"] > 200000:
        issues.append(Issue(
            "warning",
            f"this writes {derived['traj_frames']:,} frames. Analysis will be slow and "
            f"the file large; a longer interval loses nothing for the metrics on the "
            f"analysis tab.",
            "traj_interval",
        ))
    return issues


def _rule_job(cfg: dict[str, Any]) -> list[Issue]:
    name = str(cfg.get("job_name") or "")
    if not _JOB_NAME_RE.match(name):
        return [Issue(
            "error",
            "the job name becomes a directory and a filename on your machine, so it "
            "must start with a letter or digit and contain only letters, digits, "
            "hyphens and underscores (64 characters maximum).",
            "job_name",
        )]
    return []


def _rule_input(cfg: dict[str, Any]) -> list[Issue]:
    issues: list[Issue] = []
    source = cfg.get("input_source")
    if source in ("rcsb", "opm"):
        pdb_id = str(cfg.get("pdb_id") or "").strip()
        if not re.match(r"^[0-9][A-Za-z0-9]{3}$", pdb_id):
            issues.append(Issue(
                "error",
                "a PDB ID is four characters beginning with a digit, for example 1AKI.",
                "pdb_id",
            ))
    if source == "alphafold" and not str(cfg.get("uniprot_id") or "").strip():
        issues.append(Issue("error", "a UniProt accession is required.", "uniprot_id"))
    if cfg.get("use_membrane") and cfg.get("membrane_orientation") == "preoriented" \
            and source not in ("upload", "opm"):
        issues.append(Issue(
            "warning",
            "'already oriented' assumes the structure's z axis is the membrane normal. "
            "Structures straight from the RCSB are in a crystal frame and almost never "
            "are; fetching from OPM instead gets you the correct alignment.",
            "membrane_orientation",
        ))
    return issues


RULES = (
    _rule_job, _rule_input, _rule_forcefield, _rule_solvation, _rule_box,
    _rule_switching, _rule_timestep, _rule_protocol, _rule_output,
)


# ===========================================================================
#  Validation
# ===========================================================================


def validate(
    raw: dict[str, Any],
    *,
    strict_unknown: bool = False,
    form_post: bool = False,
) -> ValidationResult:
    """Coerce and check a raw submission. Returns a complete config either way.

    The returned config is always complete and usable for rendering the form
    back to the user, even when there are errors -- so a failed submission does
    not lose what they typed.

    ``form_post`` must be set when ``raw`` came from an HTML form. An unticked
    checkbox is simply absent from a form POST, so absence has to mean False
    there -- whereas for a partial config from a file or the CLI, absence means
    "leave it at the default". Guessing which case we are in from the shape of
    the input gets this wrong exactly when it matters, so the caller says.
    """
    issues: list[Issue] = []
    cfg = opts.defaults()

    # Unknown keys warn and are dropped. A typo should not stop a run, but it
    # should not be silently honoured either.
    for key in raw:
        if key not in opts.BY_ID and not key.startswith("_"):
            issues.append(Issue(
                "error" if strict_unknown else "warning",
                f"unknown option {key!r} -- ignored. Options are declared only in "
                f"flexappeal/options.py.",
            ))

    for opt in opts.OPTIONS:
        if opt.id not in raw:
            if opt.widget == "checkbox" and form_post:
                cfg[opt.id] = False
            continue
        try:
            cfg[opt.id] = _coerce(opt, raw[opt.id])
        except (TypeError, ValueError):
            issues.append(Issue(
                "error",
                f"{raw[opt.id]!r} is not a valid {opt.widget} value.",
                opt.id,
            ))

    cfg = normalise(cfg)

    # Per-option checks, skipping anything the config has made irrelevant.
    for opt in opts.OPTIONS:
        if not is_active(opt, cfg):
            continue
        value = cfg[opt.id]

        if opt.choices:
            allowed = {c.value for c in active_choices(opt, cfg)}
            values = value if opt.widget == "multiselect" else [value]
            for v in values:
                if v not in allowed and not (opt.widget == "multiselect" and v == "*"):
                    issues.append(Issue(
                        "error",
                        f"{v!r} is not an available choice here. Available: "
                        f"{', '.join(sorted(allowed))}.",
                        opt.id,
                    ))

        if opt.widget in ("number", "int") and isinstance(value, (int, float)):
            if opt.minimum is not None and value < opt.minimum:
                issues.append(Issue(
                    "error",
                    f"{value:g} is below the minimum of {opt.minimum:g}"
                    f"{' ' + opt.units if opt.units else ''}.",
                    opt.id,
                ))
            if opt.maximum is not None and value > opt.maximum:
                issues.append(Issue(
                    "error",
                    f"{value:g} is above the maximum of {opt.maximum:g}"
                    f"{' ' + opt.units if opt.units else ''}.",
                    opt.id,
                ))

    for rule in RULES:
        issues.extend(rule(cfg))

    return ValidationResult(cfg, issues)


# ===========================================================================
#  Derived quantities
# ===========================================================================

# Rough atoms-per-nm3 for water at 1 g/mL: 33.4 molecules × 3 atoms.
_ATOMS_PER_NM3_WATER = 100.0
# Bytes per atom per frame. XTC compresses to roughly 3 bytes/coordinate;
# DCD stores three 4-byte floats.
_BYTES_PER_ATOM = {"xtc": 9.0, "dcd": 12.0, "hdf5": 12.0}


def derive(cfg: dict[str, Any]) -> dict[str, Any]:
    """Numbers both the Prepare form and the run script need.

    Deliberately arithmetic only -- no OpenMM import, no structure parsing --
    so the browser can call this on every keystroke and the droplet never has
    to do real work to answer "how big will this be?".

    Atom counts depend on the actual structure, so anything that needs them
    reads ``cfg['_estimated_atoms']``, which ``inspect`` fills in once a
    structure is loaded. Until then the estimates are order-of-magnitude.
    """
    dt_fs = float(cfg.get("timestep") or 4.0)
    dt_ps = dt_fs / 1000.0

    prod_ps = float(cfg.get("production_duration") or 0.0) * 1000.0
    equil_ps = float(cfg.get("equilibration_duration") or 0.0)
    heat_ps = float(cfg.get("heat_duration") or 0.0)

    def steps(ps: float) -> int:
        return int(round(ps / dt_ps)) if dt_ps > 0 else 0

    traj_interval_ps = float(cfg.get("traj_interval") or 10.0)
    state_interval_ps = float(cfg.get("state_interval") or 10.0)
    checkpoint_ps = float(cfg.get("checkpoint_interval") or 100.0)

    atoms_total = int(cfg.get("_estimated_atoms") or 0)
    atoms_solute = int(cfg.get("_solute_atoms") or 0)

    selection = cfg.get("traj_selection")
    if selection == "all":
        atoms_saved = atoms_total
    elif selection in ("protein_ligand", "protein"):
        atoms_saved = atoms_solute
    elif selection == "protein_heavy":
        atoms_saved = int(atoms_solute * 0.5)
    elif selection == "backbone":
        atoms_saved = int(atoms_solute * 0.13)
    elif selection == "ca":
        atoms_saved = int(atoms_solute * 0.033)
    else:
        atoms_saved = atoms_solute

    traj_frames = int(prod_ps / traj_interval_ps) if traj_interval_ps > 0 else 0
    bytes_per_atom = _BYTES_PER_ATOM.get(str(cfg.get("traj_format")), 12.0)
    traj_bytes = int(traj_frames * max(atoms_saved, 1) * bytes_per_atom)

    replicates = max(1, int(cfg.get("replicates") or 1))

    return {
        "timestep_ps": dt_ps,
        "heat_steps": steps(heat_ps),
        "equilibration_steps": steps(equil_ps),
        "production_steps": steps(prod_ps),
        "total_steps": steps(heat_ps + equil_ps + prod_ps) * replicates,
        "total_ns": (heat_ps + equil_ps + prod_ps) / 1000.0 * replicates,
        "traj_interval_steps": max(1, steps(traj_interval_ps)),
        "state_interval_steps": max(1, steps(state_interval_ps)),
        "checkpoint_interval_steps": max(1, steps(checkpoint_ps)),
        "traj_frames": traj_frames * replicates,
        "atoms_saved": atoms_saved,
        "traj_bytes": traj_bytes * replicates,
        "traj_size_human": _human_bytes(traj_bytes * replicates),
        "replicates": replicates,
    }


def estimate_wall_time(cfg: dict[str, Any], ns_per_day: float | None = None) -> dict[str, Any]:
    """Estimate how long the run will take, given a throughput.

    ``ns_per_day`` comes from the bundle's own platform benchmark when one has
    been run. Before that we fall back to a crude scaling fitted to Apple
    Silicon CPU-platform throughput, which is honest about being approximate --
    the form labels it "rough" rather than presenting it as a prediction.
    """
    derived = derive(cfg)
    atoms = int(cfg.get("_estimated_atoms") or 0)

    if ns_per_day is None:
        if atoms <= 0:
            return {"ns_per_day": None, "hours": None, "basis": "unknown",
                    "human": "load a structure to estimate"}
        # ~250 ns/day at 25k atoms with a 4 fs step, scaling inversely with atom
        # count and linearly with timestep. Order-of-magnitude only.
        ns_per_day = 250.0 * (25000.0 / max(atoms, 1000)) * (float(cfg.get("timestep") or 4.0) / 4.0)
        basis = "estimated"
    else:
        basis = "benchmarked"

    hours = (derived["total_ns"] / ns_per_day) * 24.0 if ns_per_day > 0 else None
    return {
        "ns_per_day": round(ns_per_day, 1),
        "hours": round(hours, 1) if hours is not None else None,
        "basis": basis,
        "human": _human_duration(hours) if hours is not None else "unknown",
    }


def _human_bytes(n: int) -> str:
    if n <= 0:
        return "0 B"
    units = ("B", "KB", "MB", "GB", "TB")
    i = min(int(math.log(n, 1024)), len(units) - 1)
    return f"{n / (1024 ** i):.1f} {units[i]}"


def _human_duration(hours: float) -> str:
    if hours < 1:
        return f"{hours * 60:.0f} minutes"
    if hours < 48:
        return f"{hours:.1f} hours"
    return f"{hours / 24:.1f} days"
