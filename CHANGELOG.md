# Changelog

Dated history. Goals and non-goals live in `PROJECT_PLAN.md`; the roadmap is the
checkbox list in `README.md`.

## 2026-07-27

### Added

- **A trajectory clip on every tab.** Each of the three tabs now opens with its
  panel at 75% and a new 25% panel beside it, playing the Example run's lysozyme
  as a cartoon ribbon. One partial (`templates/_clip.html`) renders it in all
  three places, so they cannot drift.
- **Deferred, cheap loading for it.** 396 KB of H.264 (25 s, 480x552, 20 fps,
  `-tune animation`, no audio), which beat both VP9 and AV1 at matched SSIM on
  this flat cel-style artwork -- the AV1 file was 15% larger at the same
  quality. The markup ships a 13 KB WebP poster and no source at all;
  `static/clip.js` attaches the video only after `window.load`, so it never
  competes with the page's own assets, and never downloads it at all for a
  visitor who has asked for reduced motion. Tested, including a size budget so a
  future re-export cannot quietly land the 15 MB original on every page.
  It is the first 25 s of the 50 s playback, so it covers half the run rather
  than all of it: doubling the speed to fit the whole trajectory into the same
  25 s cost 58% more bytes, and the caption says which it is. The trajectory is
  a real one, so there is no seamless loop point anywhere in it -- measured, not
  assumed, and mirroring it into a palindrome would have meant showing the
  molecule running backwards half the time.

## 2026-07-26

### Added

- **Example tab.** A real 10 ns run of hen egg-white lysozyme (PDB 1AKI, 19,433
  atoms with water), committed to `examples/lysozyme_10ns/` and rendered through
  the same code path as an upload: `plots.build_all()` for the panels, the same
  Mol\* viewer reading the topology and trajectory out of the same `.fxa`, and
  option tables generated from the registry. It cannot drift from the product,
  because a change that breaks the Analysis tab breaks this page in CI.
  `bash examples/lysozyme_10ns/rebuild.sh` regenerates every artefact.
- **Terminal captures** (`scripts/terminal_capture.py`). Turns a `script(1)`
  recording of a run into HTML for the Example tab -- the bytes the scripts
  actually wrote, replayed through the library that produced them, with the
  progress bars caught mid-stage rather than finished. Both phases are captured:
  `run.py` and `analyse.py`. Home directories are shortened to `~`, which is the
  only edit made and is disclosed on the page.
- **Rich terminal output** (`flexappeal/console.py`). The CLI and both generated
  scripts share one palette and print through rich: tables for structure
  reports, determinate bars for the integration stages, an indeterminate spinner
  for minimisation, and a live memory readout beside every bar. Off a terminal
  it degrades to plain milestone lines with no escape sequences.
- **Memory and swap readouts** during a run, with a warning the first time a run
  pushes the machine into swap -- the failure that costs an order of magnitude
  with no other symptom.
- **Growth projection** (`watch_growth`). Samples resident size during
  production, extrapolates to the end of the run, and warns if the projection
  will not fit, naming the remedy: the bundle resumes from its checkpoint, so
  restarting continues the run in a fresh process.
- **Trajectory playback on load.** The structure panel loops the trajectory,
  one pass every 50 seconds; skipped under `prefers-reduced-motion`.

### Fixed

- **OpenCL was never selected on Apple Silicon.** Apple's implementation is
  single-precision only; the default request for mixed precision made context
  creation fail with a message naming the platform rather than the property, and
  the benchmark fell back to the CPU with nothing in the output to say so.
  Measured: 29.5 ns/day on CPU against 206 on OpenCL for the same system.
  Unsupported properties are now dropped and the substitution announced.
- **Example session expiry.** The Example tab cached a scratch session token and
  validated it with a helper that aborts with 400 when the directory is missing
  -- correct for a token off a form, wrong when the sweeper is expected to have
  removed it. The page would have worked for four hours after each deploy and
  then started failing on a timer.
- **HTML was heuristically cached.** Flask sends no `Cache-Control` on a
  rendered template, so browsers invented their own freshness lifetime and held
  the page -- which pinned the previous `?v=` asset URL, itself served
  `immutable` for a year. Every CSS and JS change was invisible to returning
  visitors for days. HTML now sends `no-cache`.

### Known

- **OpenMM's OpenCL platform leaks** about 3 kB per integration step on an
  M1 Max, measured with no reporters attached. The 10 ns reference run ended
  8 GB up; a 100 ns run would need roughly 80 GB. Neither
  `Context.reinitialize(preserveState=True)` nor rebuilding the `Simulation`
  from a saved `State` reclaims it. Restarting the process does, and the bundle
  resumes from its checkpoint.

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

### Licensing

- **MIT**, with third-party notices in `LICENSE`. Compatible with everything
  FlexAppeal depends on: the two libraries this repository redistributes
  (Plotly.js and Mol\*) are both MIT and their notices are retained, while
  MDTraj (LGPL-2.1+) and gemmi (MPL-2.0) are imported unmodified through their
  published interfaces, which both licences permit.
- **MDAnalysis removed.** It was declared in `requirements.txt`, the development
  environment and every generated bundle, and nothing ever imported it. MDTraj
  covers every metric FlexAppeal computes and owns the trajectory formats, so
  carrying MDAnalysis cost the droplet an install and every bundle an extra
  dependency in its environment solve, for nothing.

### Known limitations

Metal cofactors are refused. Ligand protonation comes from the CCD's fixed state.
Area per lipid includes the protein cross-section. Single-component membranes
only. No free-energy methods. OpenMM has no official Metal back end.
