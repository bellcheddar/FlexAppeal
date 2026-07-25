# FlexAppeal: project plan

Goals, non-goals and the reasoning behind the architecture. The checkable To Do
list lives in `README.md`; dated history lives in `CHANGELOG.md`. This document
explains *why*, and is updated when a decision changes, not when work completes.

---

## The problem

OpenMM is an excellent simulation engine with a Python API rather than a user
interface. Setting up a correct protein MD run means writing a script that
stitches together PDBFixer repair, force-field selection, solvation, system
construction, integrator and barostat choice, a staged minimise/heat/equilibrate/
produce protocol, and a set of reporters. That is roughly 120 interacting
parameters, and a large fraction of the ways to get it wrong produce a
trajectory that runs to completion and looks plausible:

- solvent padding smaller than the non-bonded cutoff, so the protein interacts
  with its own periodic image
- a 4 fs timestep without hydrogen mass repartitioning, so total energy drifts
  upward all run
- CHARMM36 with AMBER's cutoff scheme, so the force field is used outside the
  conditions it was fitted to
- a membrane protein built into a bilayer sideways because the input was in a
  crystal frame rather than an OPM frame

Analysing the result is then a second, separate pile of MDTraj and MDAnalysis
scripting, repeated from scratch for every project.

## What FlexAppeal does

Two tabs and a file that travels between them.

**Prepare** turns the option surface into a guided form with real defaults and
real explanations, validates the combination against the physics rules above,
and emits a single self-contained `.command` file.

**The bundle** runs on the user's own machine. It bootstraps its environment,
benchmarks the available OpenMM platforms, runs the simulation, then runs the
full analysis against the complete trajectory and packs a compact `.fxa` results
file.

**Analysis** consumes that `.fxa` and renders interactive Plotly panels and a
Mol\* structure viewer, with bounded on-demand re-analysis for the richer bundles.

## Non-goals

Stated explicitly, because each of these is a reasonable thing to want and the
answer is still no:

- **FlexAppeal does not run simulations.** The hosted app never invokes OpenMM.
  This is not a limitation to be lifted later; it is the design. MD needs hours
  of dedicated compute, and the droplet has 3.8 GB of RAM shared with four other
  applications.
- **Not a trajectory archive.** The `.fxa` is an analysis payload, not storage.
  The full trajectory stays on the user's machine, and the app never asks for it.
- **Not free-energy methods.** No FEP, no umbrella sampling, no metadynamics, no
  replica exchange. Each is a project in itself and none is "standard protein MD".
- **Not a structure predictor or a docking tool.** Coordinates come in from
  somewhere else. BoltzMaker is the sibling project for that.
- **Not a general workflow engine.** One protocol shape -- minimise, heat,
  equilibrate, produce -- with every parameter exposed, rather than an arbitrary
  stage graph.
- **No mixed-lipid membranes.** `Modeller.addMembrane` builds a single-component
  bilayer. Anything more complex should be built in CHARMM-GUI and uploaded.

## Architecture decisions

### The option registry is the whole design

`flexappeal/options.py` declares every parameter once, as data. The form, the
validator, the script template and the documentation all read from it. The
alternative -- a hand-written form, a hand-written validator and a hand-written
template -- drifts within a week, and the drift is invisible until a user hits it.

The corollary is a hard rule: **adding an OpenMM option means editing exactly one
file.** If a change requires touching a second, the abstraction is wrong and the
abstraction should be fixed rather than worked around.

`tests/test_openmm_symbols.py` resolves every declared API symbol against the
real installed OpenMM, so upstream renames fail a test here rather than a user's
run there.

### Predicates unify three separate decisions

A `requires` string decides whether the browser shows a field, whether the
validator demands it, and whether the template emits it. Those three can then
never disagree. The evaluator is an AST whitelist rather than `eval`, which costs
about forty lines and removes an entire category of question.

### Two layers of strictness

An unknown key warns and is dropped; an out-of-range or physically inconsistent
value is fatal. A typo should not stop someone's run, but a wrong value must not
be silently honoured -- it becomes a wrong trajectory eight hours later, on
someone else's laptop, where nobody can see it.

### Analysis runs where the trajectory is

The bundle computes the full metric panel against the complete trajectory,
because that is the only place the complete trajectory exists. What travels back
is metrics plus a decimated trajectory, sized by a user-chosen tier. This keeps
uploads in the tens of megabytes rather than the gigabytes, and means no metric
is computed on decimated data.

The hybrid `/reanalyse` route exists for the case the standard panel did not
anticipate -- a custom selection, an extra metric -- and is bounded by an
atom-by-frame budget and a concurrency cap of one.

### Honesty about Apple Silicon

OpenMM has no official Metal platform. The `openmm-metal` plugin is third-party
and experimental, and Apple has deprecated OpenCL, though it still works on
M-series hardware. Rather than promise a back end, the bundle benchmarks whatever
is actually present and picks the fastest, recording the result in the manifest.

The reliable performance levers on this hardware are hydrogen mass
repartitioning with a 4 fs timestep, mixed precision, and a thread count matched
to the performance cores. The defaults reflect that, and the help text says so.

## Build phases

| Phase | Deliverable | State |
|---|---|---|
| 0 | Repo skeleton, option registry, validation, pytest scaffolding | done |
| 1 | Structure introspection and the Prepare tab | done |
| 2 | Bundle generator and runtime templates | done |
| 3 | Local analysis pass and the `.fxa` contract | done |
| 4 | Analysis tab: Plotly panels and the Mol\* viewer | done |
| 5 | Ligand parameterisation | done |
| 6 | Membrane construction | done |
| 7 | Hybrid server-side re-analysis | done |
| 8 | Deployment to flexappeal.mdeller.com | kit ready, blocked on DNS |

All eight phases are built and tested. Deployment is blocked on one external
step: `flexappeal.mdeller.com` has no DNS record, and certbot proves domain
control over HTTP, so the TLS step cannot succeed until an A record points at
the droplet. `deploy/provision.sh` detects this and says so rather than leaving
the site serving another application's certificate.

## Deployment shape

Port 8004 on the shared droplet (8000 AlphaFraud, 8001 chem_sage, 8002 chatPDB,
8003 BoltzMaker). Flask app factory, gunicorn, systemd, nginx, certbot -- the
deploy kit is forked from AlphaFraud's, which is the most evolved of the four.
The service carries a `MemoryMax` cgroup limit, because this app runs MDTraj and
must not be able to take the box down with it.

---

## How this document is maintained

Update it when a decision changes or a non-goal is added or removed. Do not use
it as a progress log -- that is what the To Do list in `README.md` and the dated
entries in `CHANGELOG.md` are for.
