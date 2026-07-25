# Changelog

Dated history. Goals and non-goals live in `PROJECT_PLAN.md`; the roadmap is the
checkbox list in `README.md`.

## 2026-07-25

First build, phases 0 to 8.

### Added

- **Option registry** (`flexappeal/options.py`). 115 OpenMM parameters across 14
  groups, declared once as data. The form, the validator, the script template
  and `docs/options.md` all read from it.
- **Validation** (`flexappeal/schema.py`). An AST-whitelist predicate evaluator
  shared by the form, the validator and the templates, plus nine physics rules.
- **Structure introspection** (`flexappeal/structure.py`). gemmi-based: chains,
  gaps, disulfides, cofactors, metals, membrane markers, solvated-size estimates.
- **Prepare tab.** Upload or fetch from the RCSB, OPM or the AlphaFold DB, then
  configure against what was found, with live size and wall-time readouts.
- **Run bundles** (`flexappeal/bundle.py`). A self-extracting `.command` carrying
  the structure, a pinned environment, a readable OpenMM script, the analysis
  pass and its own installer.
- **Analysis tab.** Nine panels, six convergence small multiples, stat tiles, a
  vendored Mol\* viewer and a table view.
- **Ligands.** Chemistry from the RCSB Chemical Component Dictionary transferred
  onto deposited coordinates, parameterised with OpenFF Sage, GAFF2 or espaloma.
- **Membranes.** `Modeller.addMembrane` with seven lipids, OPM-oriented input,
  the membrane barostat, and area-per-lipid, thickness and order-parameter panels.
- **Server-side re-analysis** (`flexappeal/analysis.py`). Bounded, one job at a
  time, in a detached subprocess.
- **Deploy kit** (`deploy/`). Port 8004, gunicorn, systemd with a cgroup memory
  limit, nginx with rate limits, certbot with the HTTP/2 patch.

### Fixed during the build

These were all found by running the code rather than reading it, and each would
have shipped silently.

- **Periodic wrapping made distances wrong by a box length.** Trajectories are
  written with `enforcePeriodicBox`, which wraps every molecule into the primary
  cell independently, so a bound ligand measured 67.8 A from a contact that is
  really 7.5 A. `compute_contacts` applies the minimum image convention itself,
  so the contact map was entirely correct beside the wrong distances.
- **`Modeller.addSolvent(model=)` accepts only five geometry names.** OPC and the
  force-balance water models are not among them, and passing an unrecognised name
  does not raise: it silently produces the wrong water.
- **A PDB residue name field is three characters**, so POPC is stored as POP.
  Every `resname POPC` selection matched nothing and all three membrane metrics
  returned empty, with no error anywhere.
- **`Chem.AddHs` leaves new hydrogens without residue metadata**, so a ligand's
  hydrogens formed a separate UNK residue and MDTraj's contact analysis failed on
  a residue with zero heavy atoms.
- **Benchmarking before minimising sent every platform to NaN**, because
  `addSolvent` packs water without regard to overlap.
- **MDTraj's `image_molecules` anchor heuristic is unsatisfiable** for a protein
  plus a ligand: it looks for a molecule with more atoms than the largest one.
- **The view function `def analysis()` shadowed the imported `analysis` module.**
- **`tests/fixtures/*.fxa` was caught by the repo's own `*.fxa` ignore rule**, so
  the entire analysis suite would have skipped on a fresh clone.
- **Generated READMEs reported no size or wall-time estimate**, because
  `schema.validate` rebuilds from defaults and dropped the internal keys.
- **`deploy.sh` had an escaping error** that made it fail to parse at all.
- **Three brand accent colours failed the chart-palette checks** (lightness band
  and contrast floor) and were re-stepped for data visualisation.

### Known limitations

Metal cofactors are refused. Ligand protonation comes from the CCD's fixed state.
Area per lipid includes the protein cross-section. Single-component membranes
only. No free-energy methods. OpenMM has no official Metal back end.
