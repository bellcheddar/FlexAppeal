# 🧬 FlexAppeal

> **Turn the whole OpenMM option surface into a form, run the simulation on your own machine, and explore what comes back.**

[![live](https://img.shields.io/badge/live-flexappeal.mdeller.com-00897B?logo=icloud&logoColor=white)](https://flexappeal.mdeller.com) ![python](https://img.shields.io/badge/python-3.11-3776AB?logo=python&logoColor=white) ![OpenMM](https://img.shields.io/badge/OpenMM-8.5-467FF7) ![MDTraj](https://img.shields.io/badge/analysis-MDTraj-00897B) ![Mol*](https://img.shields.io/badge/viewer-Mol*%205.11-9b51e0) ![tests](https://img.shields.io/badge/tests-417-00897B) ![licence](https://img.shields.io/badge/licence-MIT-467FF7) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://flexappeal.mdeller.com" target="_blank" rel="noopener noreferrer">flexappeal.mdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/FlexAppeal" target="_blank" rel="noopener noreferrer">bellcheddar/FlexAppeal</a></td>
</tr>
</table>

---

FlexAppeal is a four-tab web application for protein molecular dynamics. The **Prepare** tab turns the OpenMM option surface into a guided form and emits a single self-contained file. You run that file on your own machine, where it installs its environment, benchmarks the available compute platforms, simulates, and analyses the full trajectory. It produces one results file, which the **Analysis** tab reads and renders as interactive plots and a three-dimensional structure viewer. The **Example** tab is a real 100 ns run of T4 lysozyme with benzene bound in its buried cavity, rendered through that same path, so what the documentation claims and what the software does cannot drift apart. The **Runs** tab is a library of kept results files: open any row for the same panels an upload gets, with each panel's description tucked under its own chart rather than above it.

**Why it matters:** setting up a correct MD run means stitching together structure repair, force-field selection, solvation, system construction, integrator and barostat choice, a staged equilibration protocol, and a set of reporters, across roughly 120 interacting parameters. A large fraction of the ways to get that wrong produce a trajectory that runs to completion and looks entirely plausible: solvent padding smaller than the non-bonded cutoff, so the protein interacts with its own periodic image; a 4 fs timestep without hydrogen mass repartitioning, so total energy climbs all run; CHARMM36 used with AMBER's cutoff scheme. FlexAppeal encodes those rules as validation rather than leaving them to be remembered. It is useful for: setting up a first simulation without a week of reading, reproducing someone else's protocol exactly, sweeping a parameter across replicates, and getting a standard analysis panel without writing MDTraj by hand.

The hosted application never simulates. That is a design decision, not a limitation to be lifted later: MD needs hours of dedicated compute, and your trajectory should not have to travel anywhere.

---

## ✨ Features

| | |
|---|---|
| **A worked example** | A real 100 ns lysozyme-benzene run, committed and rendered by the same code that renders an upload |
| **Every OpenMM option, explained** | 115 parameters across 14 groups, each with a full sentence saying why its default is what it is |
| **Validation that knows the physics** | Padding versus cutoff, timestep versus constraints, switching versus dispersion correction, membrane versus barostat |
| **Live estimates** | Atom count, trajectory size and wall time update as you type |
| **One file to run** | A self-extracting bundle carrying the structure, a pinned environment, a readable run script and its own installer |
| **Auditable output** | The generated `run.py` writes every value literally (`310.0 * unit.kelvin`), so it reads as a record of exactly what ran |
| **Platform benchmarking** | The bundle times the available OpenMM platforms on your real system and picks the fastest |
| **Ligands, by default** | A real ligand or cofactor in the file is kept without being asked for, with its chemistry fetched from the RCSB Chemical Component Dictionary and parameterised with OpenFF Sage, GAFF2 or espaloma. Crystallisation additives, ions and anything bearing a metal a small-molecule force field cannot type are offered rather than assumed, and the page says why |
| **Membranes** | POPC and six other lipids via `Modeller.addMembrane`, with pre-oriented structures fetched from OPM |
| **Analysis without scripting** | RMSD, RMSF, radius of gyration, SASA, DSSP, hydrogen bonds, contacts, PCA, clustering, and ligand or membrane panels |
| **Graded study packages** | `flexappeal package` collects finished runs into one file that DeepFlex renders as a graded report: pose, occupancy, interactions, energy, convergence and conformation, ranked within each receptor |
| **Interactive** | Plotly plots and a Mol\* viewer, both vendored: the page needs no third-party network access |
| **A trajectory beside the controls** | Every tab opens with its panel at 75% and the lysozyme-benzene run playing at 25%: 446 KB of H.264, fetched after the page has loaded and skipped entirely for reduced motion |
| **Reference states** | A motif's drift reported as a distance from the experimental active and inactive structures of the same protein, not as a number with no scale. Kinase residues resolved through KLIFS, receptor residues through GPCRdb generic numbering, so the same code annotates a protein it has never seen |

## 🏗️ How it works

```
  Prepare tab                    your machine                   Analysis tab
  ───────────                    ────────────                   ────────────
  structure  ─┐
  115 options ├─► flexappeal_run.command ─► install ─► simulate ─► .fxa ─► plots
  validation ─┘      (one file, ~70 KB)      benchmark   analyse    (~25 MB)  viewer
```

Three properties fall out of that shape:

- **No metric is computed on thinned data.** The bundle analyses the complete trajectory locally; only the results and a decimated copy for viewing travel back.
- **Uploads stay small.** Tens of megabytes rather than the gigabytes a raw trajectory would be.
- **Your trajectory stays yours.** Nothing in the bundle uploads anything.

The **Example** tab is that whole path already walked: a 100 ns run of T4 lysozyme L99A with benzene, whose bundle, terminal recording and results file are all in `examples/lysozyme_benzene_100ns/`, rendered by the code that renders an upload. The ligand is the point: an apo run draws the whole-system panels and not one of the eleven binding-site panels, so an apo example would leave most of the Analysis tab untested by the page whose job is to prove it works.

## 📋 Requirements

**To use the hosted app:** a browser.

**To run a bundle:** macOS on Apple Silicon, and about 4 GB of disk for the environment it installs on first run. Nothing needs to be installed in advance; the bundle bootstraps [pixi](https://pixi.sh) if you do not have it.

**To develop FlexAppeal:**

| Requirement | Notes |
|---|---|
| Python 3.11 | The ceiling for `openff-toolkit` and `openmmforcefields` on conda-forge |
| pixi | Installed by `install.sh` if absent |
| ~6 GB disk | OpenMM, AmberTools and the OpenFF toolkit |

## 🔧 Installation

```bash
git clone https://github.com/bellcheddar/FlexAppeal.git
cd FlexAppeal
./install.sh
```

The web application itself needs none of that. Its own dependencies are pip-installable with no compiled toolchain:

```bash
pip install -r requirements.txt
```

## 🚀 Usage

```bash
pixi run serve                    # http://127.0.0.1:8004
pixi run test                     # the full suite
```

The command-line interface exists so every stage can be driven without a browser:

| Command | What it does |
|---|---|
| `./FlexAppeal.py serve` | Run the web app locally |
| `./FlexAppeal.py inspect 1aki.pdb` | Report chains, gaps, disulfides, cofactors and the estimated system size |
| `./FlexAppeal.py fetch 1AKI` | Retrieve a structure from the RCSB, OPM or the AlphaFold DB |
| `./FlexAppeal.py bundle 1aki.pdb` | Build a run bundle without the browser |
| `./FlexAppeal.py unpack run.command` | List or extract a bundle's contents |
| `./FlexAppeal.py validate config.json` | Check a configuration against the registry |
| `./FlexAppeal.py package -o study.fxa run1 run2` | Collect finished runs into one `.fxa` version 2 study package |
| `./FlexAppeal.py docs` | Regenerate `docs/options.md` from the option registry |
| `./FlexAppeal.py sweep` | Remove expired scratch sessions |

### 📦 Packaging a study

A run writes an `.fxa` **version 1**: one system, one replicate, read by the
Analysis tab here. A study of several receptor-ligand complexes is an `.fxa`
**version 2**, which is a different kind of file rather than a bigger one, and it
is read by [DeepFlex](https://deepflex.mdeller.com) as a fully graded report with
a league table.

```bash
./FlexAppeal.py package -o study.fxa flexappeal_output/*/
```

It has to be a package because the central panels are comparative: the energy
sub-score is a percentile within a receptor's own cohort, and the two rank columns
are places in that cohort, so one complex has no cohort and those quantities are
undefined rather than zero. FlexAppeal runs one system per job, so packaging is a
separate, deliberate act: nothing about how a job runs changes, and a study can be
assembled from runs done weeks apart.

Sources can be finished run directories or version 1 `.fxa` files, mixed freely. A
run directory is the better source: it carries the complete trajectory rather than
a decimated payload. Two options worth knowing:

- `--energies <file.json>` carries in MM/GBSA and interaction energies. FlexAppeal
  cannot compute those, so without them the energy sub-score is left pending,
  which the command says out loud rather than writing a zero.
- `--annotate` resolves each receptor's UniProt accession over the network. Off by
  default, because it is the only part of the format that needs the internet.

The format is specified in `DEEPFLEX_FXA2_SPEC.md`. No chemistry travels in a
package: no structure strings and no connectivity, only pocket-cropped coordinates
and compound names in clear.

### ▶️ Running a bundle

```bash
chmod +x flexappeal_myrun.command
./flexappeal_myrun.command
```

The `chmod` is required: a file downloaded from a browser arrives without the executable bit. The first run installs the environment (several gigabytes, 5 to 15 minutes); later runs go straight to simulating. Re-running resumes from the last checkpoint.

## 📊 Output

Everything lands in one directory beside the bundle:

| File | What it is |
|---|---|
| `<job>.fxa` | **The results file.** Upload this to the Analysis tab |
| `prepared.pdb` | After repair, before solvation |
| `solvated.pdb` | The complete system as simulated |
| `trajectory.xtc` | The full trajectory, which stays on your machine |
| `state_data.csv` | Energies, temperature, pressure and density over time |
| `state.chk` | Checkpoint, for resuming |
| `system.xml`, `integrator.xml`, `final_state.xml` | Full serialisation, for exact reproduction or extension |
| `run_manifest.json` | Versions, platform, achieved ns/day, wall time |

## 🔬 Study panels

A version 1 `.fxa` is one system, and the Analysis tab plots how it moved. A
version 2 package is a study of several receptor-ligand complexes, and
[DeepFlex](https://deepflex.mdeller.com) renders it as a graded report. What
`flexappeal package` measures, and the panel each number feeds:

| Panel | What it carries |
|---|---|
| **Grade and six sub-scores** | pose, occupancy, interactions, energy, convergence and conformation. FlexAppeal supplies the measurements and nothing else: the weights and the bands live in DeepFlex, so a grade has one definition rather than one per writer |
| **League table** | every ligand ranked within its own receptor. A rank is a place in a cohort, which is why a package carries many complexes and not one |
| **Pose** | ligand RMSD to its starting pose (median, 90th centile, maximum and drift slope) after superposing on the pocket Calpha, and beside it the same RMSD fitted on the ligand itself. The two apart are what separate the ligand moving from the ligand changing shape |
| **Occupancy** | contacts retained per frame, the fraction of scored frames still in the pocket, exit events, and how far the ligand centroid wandered |
| **Interactions** | a per-frame fingerprint of polar, ionic and hydrophobic contacts: occupancy and mean lifetime per residue, similarity to the first frame, and a barcode of the busiest twenty-eight |
| **Convergence** | the scored window split four ways with the divergence between the blocks, plus pose clusters and whether the dominant pose changed and then held |
| **Conformation** | MMFF94 ligand strain in excess of the starting pose, the rotatable-bond count and the torsion flip rate |
| **Receptor flexibility** | Calpha RMSD, radius of gyration, per-residue RMSF with the pocket residues marked, a frame-against-frame RMSD matrix and a residue cross-correlation matrix |
| **3D viewer** | the pocket cropped to 10 A at the first and last frame, plus the trajectory as one multi-model PDB, superposed so the receptor holds still and what moves is the ligand |
| **Ligand diagram** | a 2D depiction with each atom shaded by how far it moved |

Three things hold across all of it. Every length is in angstrom and every time in
nanoseconds. Every scored number is measured from the first equilibrated frame
onward, and the file records which frame that was, so a chart and the number
printed beside it describe the same window. And a measurement FlexAppeal cannot
make is left out rather than written as zero, which DeepFlex reports as a pending
sub-score instead of reading a nought as a result.

The binding energies used to be the one gap. `flexappeal package --mmgbsa` now
computes them, and a ligand bundle does it unprompted at the end of its analysis.

It is AmberTools' `MMPBSA.py`, not a reimplementation, and that is the point. A
percentile across a cohort measured two different ways is not a percentile of
anything, so the protocol is pinned to the one that produced the numbers these
will be ranked against: single trajectory, `igb=2`, mbondi2 radii, salt at the
simulation's own ionic strength. There is no entropy term, which is the usual
choice and the usual caveat -- the result is systematically over-negative and
ranks ligands rather than predicting an affinity. Benzene in T4 lysozyme L99A
comes out at -9.8 kcal/mol here against roughly -5.2 from experiment.

How you can tell it worked: in a single-trajectory calculation the bonded terms
cancel exactly, because the receptor and the ligand are read out of the complex's
own coordinates. FlexAppeal refuses the result if they do not, because three
topologies that disagree about the molecule produce a plausible number rather than
an obvious failure.

`--energies` still carries numbers in from elsewhere, and takes precedence where
both are available. Where AmberTools is absent no energy is computed and the other
five sub-scores are unaffected.

A single complex is scored on all six. The energy sub-score is normally a percentile
among the other ligands run against the same receptor, and a lone complex has no
cohort, so it falls back to the absolute thresholds in DeepFlex's `grade.toml`: full
marks at -60 kcal/mol, zero at 0. Zero is physics, since no binding free energy means
no binding. -60 is the strong end of what that pipeline has measured.

The fallback runs only where there is no cohort, and the page says when it did. It is
a weaker instrument and the reason is worth knowing: within a single pocket these
energies do not separate good compounds from poor ones, so nothing anchors -60 beyond
the spread of the data. Collect runs into a study where you can, and the percentile
replaces it.

The Analysis tab draws the same binding-site measurements for a single run, and
grades it. Eleven panels, the 2D depiction shaded by atom mobility, a frame slider
that moves the structure viewer and every chart together, and a meter carrying the
letter, the composite score and each sub-score with the weight it actually bore.

The grade is DeepFlex's, not a second opinion. `flexappeal/grade.py` is a
transcription of its `grade.toml` and `grade.py`, and the test suite pins every
weight, ramp, band and gate against the published values so the two cannot drift.
A single complex can be graded because DeepFlex already grades one whose
sub-scores are incomplete: the missing weight is renormalised across the rest
rather than filled with a neutral value, which is exactly the case a lone complex
presents. The energy sub-score is the one that needs the cohort, being a
percentile among the other ligands run against the same receptor, so its 15% is
spread over the five that were measured and the meter says so. Both rank columns
need the cohort too, and a package is still what produces them.

The Example tab shows them on the committed 100 ns benzene run, where two of the
twelve stay empty for reasons worth stating. Benzene has no nitrogen, oxygen or
sulphur, so there are no polar contacts to occupy; and it has no rotatable bond,
so there is no torsion to flip and the conformation sub-score rests on the strain
alone. Both are properties of the molecule rather than of the measurement, and
the page says so where the panel would have been.

## 🎯 Reference states

A simulation tells you what moved. It does not tell you whether what moved went
anywhere meaningful. The **Reference states** tab answers the second question by
scoring every frame against the experimental structures of the same protein, so
a motif's drift arrives as a distance from the active and inactive ensembles.

Nothing is hand-annotated. Kinase residues are resolved through the KLIFS pocket
alignment and receptor residues through GPCRdb generic numbering, which is what
lets one implementation handle a protein it has never seen and refuse a
construct that renumbers rather than measure it at the wrong residues.

Each motif carries **two readings side by side**: how close it sits to the
nearest active and nearest inactive reference after superposing on the part of
the fold that does not move on activation, and what the published rule-based
criterion says. Where they disagree the motif is marked `conflict` rather than
averaged, because a disagreement between geometry and rule is the interesting
case rather than noise.

### Where the published thresholds did not hold

Two were recalibrated, both because the published number was fitted on crystal
structures and sits at the edge of a cluster, leaving no room for a side chain
at 300 K.

| Criterion | Published | Here | Why |
|---|---|---|---|
| DFG two-distance | 11.0 Å absolute | nearer anchor | The cut called six of ABL1's own DFG-out references neither state. The classes separate by 8.5 Å on the difference, so no threshold is needed where the data is dense. Agreement with KLIFS goes 11/17 → **17/17** |
| R-spine contact | 6.0 Å centroid | 7.7 Å | 6.0 sits a tenth of an angstrom above the top of the active range, so every trajectory frame read `broken` and the score carried a standing false conflict. The references still classify perfectly at 7.7 |

### Positive controls

| Control | Result |
|---|---|
| DFG label against KLIFS | 17 / 17 |
| αC label against KLIFS | 17 / 17 |
| TM3–TM6 state separation, leave-one-out | 100 / 100 at 11.39 Å |
| R-spine label after recalibration | 5 / 5 active, 14 / 14 inactive |
| Each start structure scores as its own state | 4 / 4 |
| Bilayer inside its published windows | 67.7 Å² per lipid, 38.6 Å thick, every chunk inside |
| An inactive start stays inactive | **pending**: no inactive-start run has finished yet |

The last row is listed rather than left out. It is the one control of the six
that has not been decided, because the run that would decide it is third in a
queue that is still draining, and a table of five passes with the undecided one
missing would read as six.

### Reuse

The same analysis runs on anybody's trajectory:

```bash
python FlexAppeal.py refstates annotate \
    --topology system.pdb --trajectory run.dcd \
    --start start.pdb --family auto --out annotated.json
```

`--family auto` reads the structure's own sequence and returns `unknown` rather
than defaulting, because annotating a kinase as a receptor would resolve every
motif to the wrong residues and still produce numbers. `--start` restores the
deposited numbering, which OpenMM discards on solvation.

## 🧪 Validation

The physics is checked against known values rather than against itself. Every fixture in the test suite is a real results file from a real run.

| System | Measured | Expected |
|---|---|---|
| Lysozyme (1AKI) | R<sub>g</sub> 14.13 Å | ~14.3 Å |
| Lysozyme | Helix content 39% | ~40% |
| Lysozyme | 4 disulfides at 6–127, 30–115, 64–80, 76–94 | Exactly those four |
| Trypsin + benzamidine (3PTB) | Asp189 contact 100% of frames | The defining S1 salt bridge |
| Trypsin + benzamidine | Ser190, Trp215, Gly216, Gly219, Gly226 | The canonical S1 pocket |
| Glycophorin A in POPC (1AFO) | Bilayer thickness 3.64 nm | 3.7–4.0 nm |
| Glycophorin A in POPC | Order parameters 0.20–0.25 plateau, falling to 0.12 | The fluid-bilayer profile |

A further 30 tests resolve every OpenMM symbol the option registry names against the installed package, and 8 more actually run `addSolvent` for each water model rather than merely checking that its file loads.

## 🧭 Caveats

Stated plainly, because each is a real limit rather than a rough edge:

- **Metal cofactors are refused.** Haem, zinc fingers and iron-sulfur clusters cannot go through GAFF, OpenFF or espaloma, which are organic force fields with no transition-metal parameters. FlexAppeal says so at build time instead of failing inside your run.
- **Ligand protonation comes from the RCSB chemical definition**, which is deposited in one fixed state. That state is frequently not the dominant one at the pH the protein is simulated at. The formal charge used is printed in the run log.
- **Area per lipid counts the protein's cross-section as membrane area.** It is a convergence diagnostic, not a number to compare against a pure-bilayer literature value.
- **OpenMM has no official Metal back end.** The `openmm-metal` plugin is third-party and experimental, and Apple has deprecated OpenCL. Rather than promise a back end, the bundle benchmarks what is actually present. On an M1 Max the reliable levers are hydrogen mass repartitioning with a 4 fs timestep, mixed precision, and a sensible thread count.
- **Single-component membranes only.** `Modeller.addMembrane` builds one lipid type. Anything more complex should be built in CHARMM-GUI and uploaded.
- **No free-energy methods.** No FEP, umbrella sampling, metadynamics or replica exchange. Each is a project in itself.
- **The reference-state timescale is the honest limit.** The β2AR active-to-inactive relaxation these runs look for was observed on the **microsecond** scale. Its absence at 50 ns is an expected limitation and is reported as one, not as a negative result. For kinases the whole-protein score understates a real state call, because motifs of unequal discriminating power carry equal residue weight: read the per-motif table rather than the headline number.
- **POPC only, and apo cores.** No cholesterol, no palmitoylation, and no ligand in the reference-state core runs. Each would change the answer and each is a separate build.

## 🧱 Stack

| Layer | Choice |
|---|---|
| Simulation | OpenMM 8.5, PDBFixer, openmmforcefields, OpenFF toolkit |
| Analysis | MDTraj, NumPy, SciPy |
| Web | Flask, Jinja, vanilla JavaScript (no framework, no build step) |
| Visualisation | Plotly 2.35, Mol\* 5.11 (both vendored) |
| Structure parsing | gemmi |
| Environments | pixi for development and for each generated bundle |

## 🛠️ Web deployment

FlexAppeal is deployed to a DigitalOcean droplet shared with four other applications, on port 8004.

```bash
cp deploy/.env.example .env      # fill in DROPLET_SSH and SERVER_NAME
# On the droplet, once, as root:
bash /opt/flexappeal/deploy/provision.sh
# From your Mac, thereafter:
bash deploy/deploy.sh
```

| File | Role |
|---|---|
| `deploy/provision.sh` | One-time and idempotent: user, venv, systemd units, nginx site, TLS |
| `deploy/deploy.sh` | Sync, install dependencies, restart, verify `/healthz` |
| `deploy/gunicorn.conf.py` | Two sync workers, 300 s timeout |
| `deploy/flexappeal-web.service` | Hardened unit with a `MemoryMax` cgroup limit |
| `deploy/flexappeal-scratch-clean.{service,timer}` | Sweeps abandoned sessions every 15 minutes |
| `deploy/nginx-flexappeal.conf` | Static assets from disk, rate limits, 250 MB body cap, and logging left on for `/static/` |
| `deploy/nginx-flexappeal-limits.conf` | The rate-limit zone, which must be declared at http scope |

The values duplicated across Python, nginx, gunicorn and systemd are annotated "MUST stay in sync" in both files, and `tests/test_deploy.py` enforces the agreement.

`/static/` deliberately keeps `access_log` on, and a test enforces that too. The [mdeller.com launcher](https://github.com/bellcheddar/mdeller-landing) counts a visit to this app by watching its access log for `/static/app.js`: a request only a browser that rendered the page ever makes, since a scanner fetches the HTML and stops. Turning logging off for that location is invisible in every other respect and simply pins the hit count at zero.

## ✅ To Do

Roadmap, newest ideas at the top. Suggestions welcome.

- [ ] Espaloma charges as a faster alternative to AM1-BCC, which takes minutes per ligand
- [ ] Replicate comparison in the Analysis tab: overlay several `.fxa` files on one panel
- [ ] Trajectory concatenation, so an extended run analyses as one piece
- [ ] Umbrella-sampling and steered-MD protocols (currently out of scope by design)
- [ ] A "reproduce this" button that rebuilds a bundle from an uploaded `.fxa`
- [x] **Trajectory playback on load.** The structure panel starts looping as soon as the frames are in, one pass every 50 seconds, skipped for anyone whose system asks for reduced motion.
- [ ] **Work around the OpenCL memory leak.** OpenMM's OpenCL platform grows by roughly 3 kB per integration step on an M1 Max, measured with no reporters attached, so it is the platform rather than anything FlexAppeal does. The 10 ns reference run ended 8 GB up. Neither `Context.reinitialize(preserveState=True)` nor destroying and rebuilding the `Simulation` reclaims it, and the run measures its own slope and warns in time to restart from the checkpoint, which does release it.

  The 100 ns benzene run since then says the warning's model is wrong, though, and in the direction that costs a user eight hours. Resident memory climbed linearly at 0.78 GB per percent to 46 GB at 64% of the run, then *fell* to 34.4 GB by the end, finishing with 8.4 MB of swap added on a 96 GB machine. The growth plateaus and the allocator gives it back. A linear fit over the first few percent does not see that, so the run projected 105 GB still to come and advised a restart it did not need. The slope needs refitting against this trace before the advice is trustworthy.
- [x] **A worked example.** A real 100 ns run of T4 lysozyme L99A with benzene bound in its buried cavity (181L, 25,463 atoms with water) is committed to the repository, and the Example tab renders it through the same code path as an upload: the panels come from `plots.build_all()`, the Mol\* viewer loads the topology and trajectory out of the same `.fxa`, and the option tables are read from the registry. So the page cannot drift from the product -- if the Analysis tab breaks, this breaks with it, in CI. The terminal captures are the bytes the run actually printed, replayed through the library that produced them. `bash examples/lysozyme_benzene_100ns/rebuild.sh` regenerates all of it, in about eight and a half hours.
- [x] **Apple Silicon runs on the GPU.** OpenMM's OpenCL platform was being disqualified on every Mac: Apple's implementation supports single precision only, the default configuration asks for mixed (correct on CUDA), and context creation fails with "No compatible OpenCL platform is available" -- which names the platform rather than the property, so the benchmark caught it and silently fell back to the CPU. Measured on solvated lysozyme, 19,433 atoms, M1 Max: **29.5 ns/day on CPU against 206 on OpenCL**, so a 10 ns run took eight hours instead of seventy minutes. Unsupported properties are now dropped rather than fatal, the substitution is announced, and the check probes the platform in front of it rather than hard-coding a list.
- [x] **Rich terminal output.** The CLI and both generated scripts share one palette and print through [rich](https://github.com/Textualize/rich): tables for structure reports, determinate bars for the integration stages, an indeterminate spinner for minimisation (whose iteration counter is not monotonic, so a percentage would run backwards), and a live memory readout beside every bar. Swap is the figure that matters -- a run that spills out of RAM slows by an order of magnitude with no other symptom -- so it is on screen throughout and warns once, loudly, if this run is what pushed the machine there. Off a terminal it degrades to plain milestone lines with no escape sequences, which is what makes a `nohup` log readable.
- [x] **Deployed to [flexappeal.mdeller.com](https://flexappeal.mdeller.com).** Port 8004 on the shared droplet, behind nginx with a Let's Encrypt certificate and HTTP/2. Verified in production rather than merely health-checked: a structure upload renders all 115 options and finds lysozyme's four disulfides, a bundle builds and its generated scripts compile, and a results file renders eleven panels. Each application on the box serves its own certificate, which is the failure the unconditional-certbot rule exists to prevent.
- [x] **Server-side re-analysis.** The Analysis tab can compute a metric the run did not, or restrict one to part of the molecule, against the trajectory packed in the results file. Bounded by an atom-by-frame budget, a separate ceiling on quadratic pairwise work, one job at a time via an expiring lock, and a cgroup memory limit. The work runs in a detached subprocess so MDTraj never blocks a worker.
- [x] **Membrane systems.** `Modeller.addMembrane` with seven lipids, pre-oriented structures from OPM, the membrane barostat, and area-per-lipid, thickness and deuterium order parameter panels.
- [x] **Ligands and cofactors.** Chemistry from the RCSB Chemical Component Dictionary transferred onto the deposited coordinates, parameterised with OpenFF Sage, GAFF2 or espaloma, with ligand RMSD and per-residue contact occupancy.
- [x] **Analysis tab.** Nine panels plus six convergence small multiples, a Mol\* viewer, stat tiles and a table view. The categorical palette is validated rather than chosen: three of the brand accents failed the contrast and lightness checks for chart use and were re-stepped.
- [x] **Self-contained run bundles.** One file carrying the structure, a pinned environment, a readable OpenMM script, the analysis pass and its own installer.
- [x] **The option registry.** 115 options declared once as data; the form, the validator, the script template and the documentation all read from it.

## 📝 Licence

MIT. See [`LICENSE`](LICENSE) for the text, and [`NOTICE.md`](NOTICE.md) for the third-party notices.

MIT is compatible with everything FlexAppeal depends on. The two libraries this
repository actually redistributes, Plotly.js and Mol\*, are both MIT, and their
copyright notices are retained as that licence requires. MDTraj (LGPL-2.1+) and
gemmi (MPL-2.0) are imported unmodified through their published interfaces,
which both licences permit from differently licensed code.

OpenMM, PDBFixer, openmmforcefields, the OpenFF toolkit and AmberTools are never
distributed as part of FlexAppeal. A generated bundle installs them into its own
environment on your machine, and the run script calls them there.

**Cite what you use.** Force-field parameters (AMBER, CHARMM36, OpenFF, GAFF),
chemical definitions from the RCSB Chemical Component Dictionary, and structures
from the RCSB PDB, OPM or the AlphaFold Database each carry their own terms and
citation expectations. FlexAppeal records the accession and the force field in
every results file's manifest so the provenance is there when you write it up.

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar/FlexAppeal" target="_blank" rel="noopener noreferrer">github.com/bellcheddar/FlexAppeal</a></td>
</tr>
</table>
