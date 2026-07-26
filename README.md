# 🧬 FlexAppeal

> **Turn the whole OpenMM option surface into a form, run the simulation on your own machine, and explore what comes back.**

[![live](https://img.shields.io/badge/live-flexappeal.mdeller.com-00897B?logo=icloud&logoColor=white)](https://flexappeal.mdeller.com) ![python](https://img.shields.io/badge/python-3.11-3776AB?logo=python&logoColor=white) ![OpenMM](https://img.shields.io/badge/OpenMM-8.5-467FF7) ![MDTraj](https://img.shields.io/badge/analysis-MDTraj-00897B) ![Mol*](https://img.shields.io/badge/viewer-Mol*%205.11-9b51e0) ![tests](https://img.shields.io/badge/tests-351-00897B) ![licence](https://img.shields.io/badge/licence-MIT-467FF7) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://flexappeal.mdeller.com" target="_blank" rel="noopener noreferrer">flexappeal.mdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/FlexAppeal" target="_blank" rel="noopener noreferrer">bellcheddar/FlexAppeal</a></td>
</tr>
</table>

---

FlexAppeal is a two-tab web application for protein molecular dynamics. The **Prepare** tab turns the OpenMM option surface into a guided form and emits a single self-contained file. You run that file on your own machine, where it installs its environment, benchmarks the available compute platforms, simulates, and analyses the full trajectory. It produces one results file, which the **Analysis** tab reads and renders as interactive plots and a three-dimensional structure viewer.

**Why it matters:** setting up a correct MD run means stitching together structure repair, force-field selection, solvation, system construction, integrator and barostat choice, a staged equilibration protocol, and a set of reporters, across roughly 120 interacting parameters. A large fraction of the ways to get that wrong produce a trajectory that runs to completion and looks entirely plausible: solvent padding smaller than the non-bonded cutoff, so the protein interacts with its own periodic image; a 4 fs timestep without hydrogen mass repartitioning, so total energy climbs all run; CHARMM36 used with AMBER's cutoff scheme. FlexAppeal encodes those rules as validation rather than leaving them to be remembered. It is useful for: setting up a first simulation without a week of reading, reproducing someone else's protocol exactly, sweeping a parameter across replicates, and getting a standard analysis panel without writing MDTraj by hand.

The hosted application never simulates. That is a design decision, not a limitation to be lifted later: MD needs hours of dedicated compute, and your trajectory should not have to travel anywhere.

---

## ✨ Features

| | |
|---|---|
| **Every OpenMM option, explained** | 115 parameters across 14 groups, each with a full sentence saying why its default is what it is |
| **Validation that knows the physics** | Padding versus cutoff, timestep versus constraints, switching versus dispersion correction, membrane versus barostat |
| **Live estimates** | Atom count, trajectory size and wall time update as you type |
| **One file to run** | A self-extracting bundle carrying the structure, a pinned environment, a readable run script and its own installer |
| **Auditable output** | The generated `run.py` writes every value literally (`310.0 * unit.kelvin`), so it reads as a record of exactly what ran |
| **Platform benchmarking** | The bundle times the available OpenMM platforms on your real system and picks the fastest |
| **Ligands** | Chemistry from the RCSB Chemical Component Dictionary, parameterised with OpenFF Sage, GAFF2 or espaloma |
| **Membranes** | POPC and six other lipids via `Modeller.addMembrane`, with pre-oriented structures fetched from OPM |
| **Analysis without scripting** | RMSD, RMSF, radius of gyration, SASA, DSSP, hydrogen bonds, contacts, PCA, clustering, and ligand or membrane panels |
| **Interactive** | Plotly plots and a Mol\* viewer, both vendored: the page needs no third-party network access |

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
| `./FlexAppeal.py docs` | Regenerate `docs/options.md` from the option registry |
| `./FlexAppeal.py sweep` | Remove expired scratch sessions |

### Running a bundle

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
| `deploy/nginx-flexappeal.conf` | Static assets from disk, rate limits, 250 MB body cap |
| `deploy/nginx-flexappeal-limits.conf` | The rate-limit zone, which must be declared at http scope |

The values duplicated across Python, nginx, gunicorn and systemd are annotated "MUST stay in sync" in both files, and `tests/test_deploy.py` enforces the agreement.

## ✅ To Do

Roadmap, newest ideas at the top. Suggestions welcome.

- [ ] Espaloma charges as a faster alternative to AM1-BCC, which takes minutes per ligand
- [ ] Replicate comparison in the Analysis tab: overlay several `.fxa` files on one panel
- [ ] Trajectory concatenation, so an extended run analyses as one piece
- [ ] Umbrella-sampling and steered-MD protocols (currently out of scope by design)
- [ ] A "reproduce this" button that rebuilds a bundle from an uploaded `.fxa`
- [x] **A worked example.** A real 10 ns run of hen egg-white lysozyme (1AKI, 19,433 atoms with water) is committed to the repository, and the Example tab renders it through the same code path as an upload: the panels come from `plots.build_all()`, the Mol\* viewer loads the topology and trajectory out of the same `.fxa`, and the option tables are read from the registry. So the page cannot drift from the product -- if the Analysis tab breaks, this breaks with it, in CI. The terminal captures are the bytes the run actually printed, replayed through the library that produced them. `bash Examples/lysozyme_10ns/rebuild.sh` regenerates all of it.
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
