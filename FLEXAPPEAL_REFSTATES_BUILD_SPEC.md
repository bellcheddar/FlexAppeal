# FLEXAPPEAL: Reference-State Examples Build Spec (ABL1 kinase + β2AR in membrane)

**Spec version:** 1.0 (2026-10-10)
**Project root:** `/Users/dellboy/Documents/Vibe_Coding/FlexAppeal/` (if the existing FlexAppeal code lives elsewhere under `Vibe_Coding/`, use that folder and log the path in `DECISIONS.md`)
**Mode:** Unattended (Claude Code on the Mac Studio). This spec is authoritative.

**One-paragraph summary.** Add two worked example systems to FlexAppeal: the human ABL1 kinase domain in water, and the human β2-adrenergic receptor (non-fusion model from GPCRdb) in an explicit POPC bilayer. Each system is simulated from two starting states (active-like and inactive). Every trajectory frame is annotated with family motifs (KLIFS for the kinase, GPCRdb generic numbering for the receptor) and secondary structure (DSSP). Each frame is then compared against ensembles of experimental active and inactive reference structures, so the user can see which state the protein is drifting towards, motif by motif. The motif annotation code is not rewritten: it is lifted from BoltzMaker's `compare-sse` (`sse_comparison/` package), which already does exactly this for single structures, and extended to trajectories.

---

## 0. Autonomy directive

### 0.1 Continuous Build Rules

1. **Run to completion without asking.** When a choice is ambiguous, pick the most defensible option, write it to `DECISIONS.md` (decision, alternatives, reason, how to reverse it) and carry on.
2. **Stop only at a Gate** (table below). A failing tool, missing package or slow step is a problem to solve, not a reason to stop.
3. **Fallback ladder for any failure:** retry once with a fix, then the documented alternative, then the next-best tool, then degrade gracefully (skip the item, flag it in the data as `status = "failed:<reason>"`, continue). One replica, one reference structure or one optional metric never blocks the pipeline.
4. **Resumable and idempotent.** Every stage writes a manifest (`data/manifests/<stage>.jsonl`). On restart, read `PROGRESS.md` and the manifests and skip completed items. MD runs write OpenMM checkpoints every 1 ns and resume from the last checkpoint. Assume the process can be killed at any moment.
5. **Log as you go.** `BUILD_LOG.md` (timestamped, one line per meaningful event) and `PROGRESS.md` (phase checklist with percentages, rewritten after every phase and at least every 30 minutes during long compute).
6. **Never idle.** While MD runs, build the analysis and app against the short smoke-test trajectories, mark missing fields "pending", and backfill automatically when the manifests show completion.
7. **Commit locally, never push.** `git commit` at the end of every phase with `phase N: <summary>`. No `git push` (Gate G4).
8. **Stay inside the project root.** No writes outside it except package caches (pixi, uv, pip, Homebrew). Never `sudo`.
9. **Keep the Mac awake:** wrap long compute in `caffeinate -dimsu <command>`.
10. **House style everywhere:** British English, no em dashes (colons or parentheses), none of groundbreaking, revolutionary, paradigm-shifting, game-changing, cutting-edge, unprecedented, seamless, leverage (verb), delve.

### 0.2 Gates (the only reasons to stop)

| Gate | Opens when | Needs from Marc |
|---|---|---|
| **G1 Compute** | The OpenMM platform smoke test fails on every local platform, or the Tier 2 set exceeds `MAX_COMPUTE_HOURS` (default 96) after automatic scope reduction | Accept reduced scope, raise budget, or approve Tier 3 on RunPod |
| **G2 Secrets** | Any credential is needed (RunPod API key, SSH keys) | The secret, via environment variable, never written into files |
| **G3 Deploy** | The feature passes QC locally and is ready for mdeller.com / ZeroGPU | Approval to rsync / push the Space |
| **G4 Publish** | Any `git push`, repo creation or public release | Approval |
| **G5 Disk** | Free space on the project volume falls below 50 GB | What to prune |
| **G6 Science** | A positive control (Section 9) still fails after the fallback ladder **and** that would invalidate a headline result (state classification, membrane integrity) | Decision to proceed with a flagged caveat |

**Gate protocol:** write `GATE_OPEN.md` (which gate, exactly what is needed, the command or answer that unblocks it, what is done, what happens next); send `osascript -e 'display notification "<one line>" with title "FlexAppeal gate <id>" sound name "Glass"'`; then keep working on everything the gate does not block, and stop only when nothing unblocked remains.

### 0.3 Launch recipe

```bash
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal
tmux new -s flexappeal
claude
# then: "Read FLEXAPPEAL_REFSTATES_BUILD_SPEC.md and execute it end to end under Section 0."
```

---

## 1. Goals

Every figure and feature maps to one of these.

1. **G-REF:** For each system, build an active-like and an inactive reference ensemble from experimental structures, annotated through KLIFS (ABL1) or GPCRdb (β2AR), with no hand-typed residue numbers anywhere.
2. **G-SIM:** Simulate each system from both starting states with replicas: ABL1 in water, β2AR in an explicit POPC bilayer.
3. **G-MOTIF:** Per frame, report the state of each named motif (kinase: DFG, αC helix, A-loop, P-loop, hinge, R-spine; GPCR: TM6, TM3–TM6 ionic lock, NPxxY, PIF, toggle switch, helix 8) and the DSSP secondary structure of every residue.
4. **G-STATE:** Per frame and per motif, report whether the conformation is closer to the active or the inactive reference ensemble, as a continuous score and a categorical label.
5. **G-SHOW:** Present all of this in FlexAppeal as two example cards that open into a trajectory explorer: Mol* with motif colouring and a per-frame state badge, a DSSP ribbon over time, 2D state maps with the reference cloud and the trajectory path, and a motif-versus-reference table.
6. **G-REUSE:** Package the annotation and comparison as a FlexAppeal module that works for any kinase or class A GPCR a user uploads, not only these two examples.

---

## 2. Hardware (Mac Studio Performance Rules, specialised)

### 2.1 Probe first (Phase 0)

`tools/hwprobe.py` writes `config/hardware.toml` from `sysctl -n machdep.cpu.brand_string`, `hw.perflevel0.physicalcpu`, `hw.perflevel1.physicalcpu`, `hw.memsize`, `system_profiler SPDisplaysDataType`, `df`, plus **every OpenMM platform available** (`openmm.Platform.getNumPlatforms()`) and the speed each reports.

### 2.2 Derive settings (`config/tuning.toml`)

| Setting | Rule |
|---|---|
| `md_platform` | Benchmark each available platform on a 20 ps run of the solvated ABL1 system; keep the fastest that passes an energy-conservation check. Expected order of preference on Apple Silicon: `OpenCL` (precision `mixed`), then `CPU`. If a `Metal` plugin platform is installed and passes, it competes like any other. `CUDA` only on RunPod |
| `md_precision` | `mixed` on GPU platforms; never `single` for production |
| `cpu_workers` | performance cores minus 2 (minimum 2); used for analysis, never alongside a CPU-platform MD run |
| `io_workers` | 8 for API fetches, with polite rate limits (KLIFS and GPCRdb: at most 2 requests/second each) and tenacity retries |
| `concurrent_md` | 1 on a single GPU (OpenCL contexts compete); analysis runs concurrently on CPU |
| `memory_headroom_gb` | 24 |

Record measured ns/day per system and platform in `config/benchmarks.toml`. Every compute budget in Section 5.6 is derived from these measurements, never assumed.

### 2.3 Run MD and analysis concurrently

```
fetch refs ─► annotate refs ─┐
                              ├─► MD queue (one GPU worker, checkpointed) ─► frames ─► analysis pool (CPU) ─► Parquet ─► SQLite + JSON
build systems ────────────────┘
```

The analysis pool polls `data/manifests/md.jsonl` and analyses each completed 10 ns chunk as it lands, so the explorer fills in while production runs.

---

## 3. Environment

**Compute:** pixi (conda-forge, `osx-arm64` + `linux-64` so the same lock runs on RunPod). **App:** whatever FlexAppeal already uses; if nothing exists yet, Flask + Jinja + gunicorn in a uv environment, vanilla JS, Plotly and Mol* vendored under `static/vendor/`.

| Package | Purpose | Fallback |
|---|---|---|
| `openmm` (>= 8.1) | MD engine (Marc's existing choice for FlexAppeal) | none: G1 if it will not run |
| `openmmforcefields` | Only if the optional ligand phase runs (Section 5.7) | skip the ligand phase |
| `pdbfixer` | Missing atoms, termini, protonation | `openmm.app.Modeller.addHydrogens` |
| `mdtraj` | DSSP (`compute_dssp`), distances, dihedrals, fast per-frame metrics | `MDAnalysis` |
| `MDAnalysis` | Membrane analysis (thickness, area per lipid), `HELANAL` helix axes | hand-rolled helix axis fit (Kahn-style local axis) |
| `gemmi` | Structure I/O, superposition (as in BoltzMaker) | `mdtraj.superpose` |
| `biopython` | Sequence alignment for numbering maps | `parasail` |
| `pandas`, `pyarrow`, `duckdb` | Analytical store | none needed |
| `scipy`, `numpy` | Geometry, clustering | none needed |
| `plotly` | Static report generation | none needed |
| `requests`, `tenacity` | API access with retries | `urllib` |
| `ambertools` (`packmol-memgen`) | Membrane builder fallback, and the only route if a cholesterol-containing bilayer is wanted later | OpenMM `Modeller.addMembrane` is primary |
| `dssp` (`mkdssp`) | Full DSSP for reference structures, matching BoltzMaker | `mdtraj.compute_dssp(simplified=False)` |
| Kincore standalone (optional) | Modi–Dunbrack dihedral labels (BLAminus etc.) | report spatial DFG labels only; dihedral label `N/A` |

---

## 4. Data acquisition

### 4.1 Reuse BoltzMaker's `sse_comparison/` (do this first)

BoltzMaker (`github.com/bellcheddar/BoltzMaker`, public, MIT) already contains a tested `sse_comparison/` package with:

- pluggable motif annotators: GPCRdb generic numbering (TM1–7, H8, ECL1–3, ICL1–3), KLIFS 85-residue pocket (hinge, gatekeeper, HRD, DFG, αC-Glu, catalytic Lys), PDBe SIFTS → Pfam fallback;
- an injectable API client with a disk cache (`--refresh-cache` bypasses it) and real recorded API responses as test fixtures;
- apo-onto-holo superposition on stable, non-binding-site residues via gemmi;
- the per-motif metrics: Cα RMSD, centroid shift, helix-axis rotation, helix kink, SSE boundary shift, flagged φ/ψ residues, DFG-in/out (Dunbrack two-distance criterion on the DFG-Phe ring, threshold 11.0 Å), αC-in/out (Lys NZ to Glu carboxylate below 4.0 Å);
- the **numbering-offset reconciliation** that was added precisely because a trimmed ABL1 kinase domain numbered 1–272 shared no residue numbers with references numbered 223–500 (2GQG): every candidate offset is scored by residue-type agreement and the best wins;
- the 70 % identity gate that refuses a reference of a different protein.

Clone BoltzMaker at the current `main` commit into `vendor/boltzmaker/` (read-only), record the commit SHA in `DECISIONS.md`, and import `sse_comparison` from there. Do **not** reimplement any of the above. FlexAppeal's new module (`flexappeal/refstates/`) wraps it and adds: trajectory iteration, reference ensembles (many structures per state, not one apo), the state-similarity score, membrane analysis and the GPCR micro-switch metrics that `compare-sse` does not yet have. If a wrapper needs a change inside `sse_comparison`, put it in a FlexAppeal subclass or adapter, never an edit to the vendored copy, and log it in `DECISIONS.md` as a candidate upstream change for BoltzMaker.

Also reuse BoltzMaker's existing ABL1 and β2AR knowledge: its `adrb2_gs_panel` example and test fixtures (a real apo β2AR structure and real KLIFS/GPCRdb responses) are the first smoke-test data for Phase 1.

### 4.2 ABL1 (kinase, water)

| Item | Rule |
|---|---|
| Identity | Human ABL1, UniProt P00519, KLIFS `kinase_ID` resolved at run time via `GET https://klifs.net/api/kinase_ID?kinase_name=ABL1&species=Human` (the swagger lists ABL1 as kinase ID 392; assert the lookup agrees, do not hard-code it) |
| Construct | Kinase domain only: the span covered by the selected start structures, trimmed to the residues present in both start states plus up to 5 flanking residues. No SH3/SH2, no myristoyl site |
| Numbering | ABL1 has two isoforms (1a/1b, 19-residue offset) and PDB entries number constructs differently. All residue identities come from KLIFS `interactions_match_residues` (KLIFS position ↔ X-ray position per structure) and the BoltzMaker offset reconciliation. Never type a residue number |
| Reference list | `GET /api/structures_list?kinase_ID=<id>`; keep human, X-ray, `resolution` ≤ 2.8 Å, `quality_score` ≥ 6, `missing_residues` ≤ 5 in the pocket. Store `DFG`, `aC_helix`, `Grich_distance/angle/rotation`, `ligand`, `allosteric_ligand` per entry |
| Active-like ensemble | KLIFS `DFG == "in"` and `aC_helix == "in"` |
| Inactive ensemble | KLIFS `DFG == "out"` (any αC), plus a separate "Src-like inactive" subset (`DFG == "in"` and `aC_helix == "out"`) kept as its own category, since it is a distinct inactive state in Abl |
| Start structures | Preferred seeds: **2GQG** (dasatinib, expected DFG-in/αC-in, already used in BoltzMaker) and **2HYY** (imatinib, expected DFG-out). Accept a seed only if KLIFS labels it as expected and its activation loop is ≥ 90 % modelled; otherwise take the best-ranked entry in that ensemble (resolution, then fewest missing A-loop residues). Log the choice |
| Ligands | Removed for the core runs (apo MD from both conformations). Ligand-bound runs are the optional phase in 5.7 |
| Missing residues | Fill short gaps (≤ 8 residues) with PDBFixer; if the A-loop gap is longer, pick the next seed. Flag every modelled residue in the data so the app can grey it |

### 4.3 β2AR (class A GPCR, membrane)

| Item | Rule |
|---|---|
| Identity | Human β2AR, GPCRdb entry name `adrb2_human`, UniProt P07550 |
| Generic numbering | `GET https://gpcrdb.org/services/residues/extended/adrb2_human/` returns `sequence_number`, `amino_acid`, `protein_segment`, `display_generic_number`, `alternative_generic_numbers` (use the `BW` scheme for labels). Map every motif through this, never by hand |
| Reference list | `GET https://gpcrdb.org/services/structure/protein/adrb2_human/` returns `pdb_code`, `state` (`Active` / `Inactive` / `Intermediate`), `resolution`, `preferred_chain`, `type`, `ligands[].function`, and `signalling_protein` where present. Keep X-ray and EM entries ≤ 3.6 Å. **Active ensemble:** `state == "Active"`. **Inactive ensemble:** `state == "Inactive"`. Keep `Intermediate` entries as a third reference set (they are useful on the state maps) |
| Start structures: non-fusion | The deposited β2AR structures carry T4 lysozyme or BRIL fusions, nanobodies or G proteins. Start structures must be GPCRdb **refined** models (fusions removed, thermostabilising mutations reverted, missing segments modelled). Not every PDB entry has one: `https://gpcrdb.org/structure/refined/2RH1` and `.../3SN6` both report "does not have a refined model". So **discover** at run time: for each entry in the active and inactive lists (best resolution first), request `https://gpcrdb.org/structure/refined/<PDB>`, parse the page for its model download link, and take the first entry per state that has one. Log the URL pattern found in `DECISIONS.md` |
| Fallback ladder for start structures | (1) GPCRdb refined model of an experimental entry; (2) GPCRdb state-specific structure model from `gpcrdb.org/structure/homology_models` for `adrb2_human` in the needed state (the listing currently shows an Inactive AlphaFold2 model; an Active one may or may not exist); (3) self-built: take the best experimental entry for that state, delete fusion/nanobody/G-protein chains, revert mutations using GPCRdb's construct annotation, rebuild short gaps with PDBFixer. Each rung down is logged and flagged in the app ("start model built locally") |
| Construct trimming | Start at the first TM1 residue minus 3 and end at the last helix 8 residue plus 3 (from `protein_segment`). ICL3: keep whatever the refined model contains; if it models more than 20 ICL3 residues, keep them but flag ICL3 as "modelled, flexible" and exclude it from every RMSD. No palmitoylation in v1 (log as a known simplification) |
| Ligands | Removed for the core runs, as for ABL1. This is deliberate: the active state without G protein or nanobody is expected to relax towards inactive on the µs timescale (Dror et al., PNAS 2011, found this even with agonist bound), which is exactly the motif change the example should show. Optional agonist/inverse-agonist runs are in 5.7 |

### 4.4 Reference structure preparation (both systems)

For every reference entry: download mmCIF from RCSB, extract the relevant chain (KLIFS `chain`; GPCRdb `preferred_chain`), remove waters, ligands and all other chains, map to the simulation numbering via the BoltzMaker offset reconciliation plus the 70 % gate, run full DSSP, compute every metric in Section 5 exactly as for a trajectory frame (same code path, a reference is a one-frame trajectory). Cache under `data/refs/<system>/<pdb>_<chain>.cif.gz` with a JSON sidecar of metrics.

---

## 5. Methods and metrics

### 5.1 System building

**ABL1 (water):** PDBFixer (missing heavy atoms, termini) → protonation at pH 7.4 via `Modeller.addHydrogens` (record every His tautomer) → Amber ff14SB (`amber14-all.xml`) + TIP3P (`amber14/tip3p.xml`) → truncated-octahedron or dodecahedron box, 12 Å padding → 0.15 M NaCl, neutralised.

**β2AR (membrane):**

1. **Orientation.** The membrane builder needs the receptor pre-oriented with the bilayer normal along z and the membrane centre at z = 0. Primary: superpose the start model's 7TM Cα onto an OPM-oriented entry of the same receptor (resolve the OPM download for a β2AR PDB entry at run time and log it). Fallback: geometric orientation (z axis = mean direction of the seven TM helix axes, extracellular to intracellular; membrane centre = midpoint of the aromatic belts defined by Trp/Tyr residues at the TM ends from `protein_segment`). Validate: Trp/Tyr side chains of TM ends within ±(10–17) Å of z = 0 on each leaflet, Arg/Lys at the intracellular interface below z = −13 Å. Failure is G6.
2. **Build.** OpenMM `Modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1.5*nanometer, ionicStrength=0.15*molar, positiveIon='Na+', negativeIon='Cl-')` with Amber ff14SB + Lipid17 (`amber14/lipid17.xml`) + TIP3P. POPC only in v1 (`addMembrane` does not provide cholesterol). Fallback builder: `packmol-memgen` (AmberTools) with POPC, then load with the same force field. If a POPC:cholesterol bilayer is wanted later, that is the route (Lipid21) and it is out of scope for v1.
3. **Size check.** Expect roughly 80–110 k atoms. Log the atom count and box; if more than 150 k, reduce `minimumPadding` to 1.2 nm and log it.

### 5.2 Simulation protocol

Common: PME, 1.0 nm real-space cutoff, `HBonds` constraints, hydrogen mass repartitioning to 1.5 amu with a 4 fs step (fallback 2 fs without HMR if energy drift exceeds 0.01 kT/ns/dof), Langevin middle integrator (friction 1 ps⁻¹). Frames saved every 100 ps to DCD (protein + lipid headgroups for the GPCR; a full-system DCD every 1 ns for restarts and membrane analysis).

| Stage | ABL1 (water) | β2AR (membrane) |
|---|---|---|
| Minimise | to convergence, protein heavy atoms restrained 10 kcal/mol/Å² | same, plus lipid P atoms restrained in z |
| Heat | 0 → 300 K over 500 ps, NVT | 0 → 303.15 K over 500 ps, NVT (POPC is fluid at 303 K) |
| Equilibrate | 1 ns NPT (MonteCarloBarostat, 1 bar), restraints halved every 250 ps to zero | 10 ns NPT with `MonteCarloMembraneBarostat` (1 bar, zero surface tension, XY isotropic, Z free), protein restraints released in six steps (10 → 5 → 2 → 1 → 0.5 → 0 kcal/mol/Å²), lipid z-restraint off after step 2 |
| Production | NPT, MonteCarloBarostat | NPT, MonteCarloMembraneBarostat as above |

**Protonation toggle (ABL1).** The DFG flip in Abl is strongly accelerated when the DFG Asp is protonated (Shan et al., JACS 2009). Run the active-start ABL1 replicas in two protonation arms: DFG-Asp standard (charged) and DFG-Asp protonated (`ASH`). The DFG Asp is KLIFS pocket position 81 (verify against the KLIFS residue definition in `sse_comparison` at build time), mapped through `interactions_match_residues`. This arm is the headline FlexAppeal feature for kinases: the same protein, one proton apart, different motif behaviour.

### 5.3 Run matrix and tiers

| System | Start | Arms | Replicas | Tier 2 length (Mac) | Tier 3 length (RunPod) |
|---|---|---|---|---|---|
| ABL1 | active-like (DFG-in) | Asp charged / Asp protonated | 3 each | 100 ns | 500 ns |
| ABL1 | inactive (DFG-out) | Asp charged | 3 | 100 ns | 500 ns |
| β2AR | active (refined) | apo | 3 | 100 ns | 1 µs |
| β2AR | inactive (refined) | apo | 3 | 100 ns | 1 µs |

- **Tier 1 (smoke):** 1 ns of every row, one replica. Must complete before any longer run starts.
- **Tier 2 (Mac Studio):** the lengths above, scheduled in priority order: β2AR active-start (the clearest expected relaxation), ABL1 protonated arm, then the rest.
- **Tier 3 (RunPod, optional):** never started without G2 (API key) and explicit budget approval via G1. Write `deploy/runpod/` (job script, pixi lock, resume logic) during the build but do not execute it.

### 5.4 Per-frame metrics: kinase (ABL1)

All residues resolved from KLIFS positions; KLIFS-style region names in the output.

| Metric | Definition | State rule |
|---|---|---|
| DFG spatial label | Dunbrack two-distance criterion as implemented in `sse_comparison`: D1 = DFG-Phe ring to αC-Glu+4 Cα, D2 = DFG-Phe ring to catalytic-Lys Cα, threshold 11.0 Å | `in` / `out` / `other` (as BoltzMaker) |
| DFG dihedral label (optional) | Kincore standalone on X-DFG backbone and Phe χ1 | BLAminus etc., or `N/A` |
| αC salt bridge | catalytic Lys NZ to αC-Glu carboxylate O (minimum) | `in` if < 4.0 Å, `out` otherwise, `unresolved` if atoms missing |
| αC displacement | αC Cα centroid shift and helix-axis rotation after core superposition | continuous |
| A-loop | Cα RMSD of the activation loop to the nearest active and nearest inactive reference; A-loop SSE content | continuous |
| P-loop (G-rich loop) | Report KLIFS-style G-rich distance/angle/rotation if the definitions in the KLIFS literature can be verified at build time; otherwise P-loop tip Cα to hinge distance plus P-loop Cα RMSD (log which) | continuous |
| Hinge | hinge Cα RMSD and backbone φ/ψ flags (BoltzMaker threshold 30°) | continuous |
| R-spine | pairwise side-chain centroid distances of the four regulatory-spine residues (β4 Leu, αC-Glu+4, DFG-Phe, HRD-His), resolved via KLIFS; "assembled" if all consecutive pairs < 6.0 Å | `assembled` / `broken` |

**Superposition core (kinase):** C-lobe helices αE, αF, αH defined by DSSP on the reference ensemble medoid, Cα only. Log the residue list.

### 5.5 Per-frame metrics: GPCR (β2AR)

All residues resolved through GPCRdb generic numbers; labels in Ballesteros–Weinstein form with the β2AR residue number beside them (for example `R3.50 (R131)`), the β2AR numbers coming from the GPCRdb residue table, never typed.

| Metric | Definition | Active direction |
|---|---|---|
| TM6 outward swing | Cα distance 3.50 – 6.34 (intracellular cleft opening); also report 3.50 – 6.30 | larger |
| Ionic lock | R3.50 CZ to E6.30 CD (minimum heavy-atom distance) | broken (larger) |
| TM6 rotation | helix-axis rotation of TM6 intracellular half vs reference medoids | per reference |
| NPxxY | Y7.53 χ1/χ2 and RMSD of N7.49–Y7.53 backbone to active vs inactive medoids; Y5.58 OH to Y7.53 OH distance | Y7.53 towards TM3, Y5.58–Y7.53 closer |
| PIF / connector | RMSD of I3.40, P5.50, F6.44 side chains to active vs inactive medoids | per reference |
| Toggle switch | W6.48 χ1 and χ2 | per reference |
| D2.50 / Na⁺ pocket | D2.50 to S3.39 and N7.49 distances; record whether a Na⁺ ion sits within 3.5 Å of D2.50 | Na⁺ present favours inactive |
| Helix 8 | helix fraction (DSSP) and angle to the membrane plane | continuous |
| TM kinks | kink angle at P6.50 and P7.50 (BoltzMaker kink metric) | continuous |
| GPCRdb activation index (optional) | Implement only if the published definition can be verified from the primary paper at build time; otherwise omit and log | per paper |

**Superposition core (GPCR):** Cα of TM1–TM4 residues with generic numbers x.38–x.58 present in every reference (the helical core least involved in activation). Log the list.

### 5.6 Secondary structure and the state score (both systems)

**SSE per frame:** `mdtraj.compute_dssp(simplified=False)` on every saved frame; per motif report helix fraction, strand fraction and SSE boundary positions. Compare each frame's per-residue SSE string with the reference ensembles: per residue, the fraction of active references and inactive references sharing the frame's SSE class.

**State score per motif per frame:**

- `d_act` = minimum motif Cα RMSD (after core superposition) to any active reference; `d_inact` likewise.
- `s = (d_inact − d_act) / (d_inact + d_act)`, so `s ∈ [−1, 1]`; +1 active-like, −1 inactive-like.
- Categorical label: `active` if `s > 0.2` and the motif's rule-based label (Sections 5.4/5.5) agrees; `inactive` if `s < −0.2` and it agrees; otherwise `intermediate`. Disagreement between RMSD and rule is shown as `conflict` (not silently resolved).
- Whole-protein score: residue-weighted mean of motif `s` values, with the per-motif breakdown always shown beside it.
- Smooth with a centred 2 ns window for display only; store raw values.

**2D state maps:**

- Kinase: D1 vs D2 (the Dunbrack plane) with every KLIFS reference plotted by its KLIFS label, trajectory paths coloured by time; second map αC salt-bridge distance vs A-loop `s`.
- GPCR: TM3–TM6 distance vs NPxxY RMSD-to-inactive (the plot style used by Dror et al. 2011), every GPCRdb reference plotted by its GPCRdb state; second map ionic lock vs connector `s`.

### 5.7 Optional phase: ligand-bound runs (off by default)

Behind a config flag `ligand_runs = false`. If enabled: ABL1 + imatinib from the DFG-out seed, ABL1 + dasatinib from the DFG-in seed; β2AR + carazolol (inactive) and + BI-167107 (active). Parameters with OpenFF Sage via `openmmforcefields` `SystemGenerator`, SMILES pulled from RCSB/ChEMBL (never hand-typed) and checked with RDKit as BoltzMaker does. Same analysis plus ligand RMSD and contact persistence. Not started unless Tier 2 core runs are complete.

### 5.8 Membrane quality checks (β2AR)

Per 10 ns chunk: area per lipid (expect roughly 60–68 Å² for POPC at 303 K with Lipid17; outside 55–75 Å² is G6), P–P bilayer thickness (expect roughly 36–41 Å), box z drift, water inside the hydrophobic core between TM helices (count only, as a sanity metric), receptor tilt relative to the bilayer normal. Shown in a small "membrane health" panel so a user can see the bilayer behaved.

---

## 6. App (FlexAppeal integration)

### 6.1 Design

Use FlexAppeal's existing design system and templates if they exist in the repo. If FlexAppeal has no UI yet, apply the marcdeller.com brand (header, palette, Inter + Roboto Mono) from the `marcs-vibe-coding` defaults and log that choice. Motif colours are fixed tokens shared by every view:

| Token | Kinase motif | GPCR motif |
|---|---|---|
| `--motif-1` | αC helix | TM6 |
| `--motif-2` | DFG | ionic lock (R3.50/E6.30) |
| `--motif-3` | activation loop | NPxxY |
| `--motif-4` | P-loop | PIF / connector |
| `--motif-5` | hinge | toggle switch W6.48 |
| `--motif-6` | R-spine | helix 8 |
| `--state-active` / `--state-inactive` / `--state-intermediate` / `--state-conflict` | state badges, map colours | same |

Colour-blind-safe palette; every state is also distinguishable by shape or label, never by colour alone.

### 6.2 Pages and features

1. **Examples tab:** two cards (ABL1 in water, β2AR in POPC), each with a thumbnail, atom count, simulated time, and a one-line finding generated from the data (for example "Active-start β2AR: TM3–TM6 distance fell from x Å to y Å in 2 of 3 replicas"), with numbers filled from SQLite, never written by hand.
2. **Trajectory explorer** (per system, replica selector, start-state selector, protonation-arm selector for ABL1):
   - Mol* viewer: cartoon coloured by motif tokens; the GPCR shows the membrane as transparent phosphate spheres (toggle); frame scrubber; reference overlay toggle (nearest active and nearest inactive reference superposed on the core, ghosted).
   - State badge strip: one badge per motif for the current frame.
   - DSSP ribbon: residues × time heatmap with motif bands marked, plus active/inactive reference consensus rows above it.
   - Motif timelines: Plotly, one small panel per motif metric, shaded bands for the active and inactive reference ranges (5th–95th percentile of each ensemble).
   - 2D state maps (Section 5.6) with reference clouds and the trajectory path; clicking a point jumps the viewer to that frame.
   - Motif table: per motif, start value, final 10 % mean, nearest reference PDB at end, `s`, label.
   - Membrane health panel (β2AR only).
3. **Reference browser:** table of every reference structure used, with KLIFS or GPCRdb labels, our computed labels, and agreement ticks (this is where the positive controls are visible to the user).
4. **About:** methods in plain English, software and data table (OpenMM, PDBFixer, MDTraj, MDAnalysis, gemmi, KLIFS, GPCRdb, OPM, BoltzMaker's `sse_comparison` with commit SHA), references verified at build time.
5. **Reuse hook (G-REUSE):** an "annotate my trajectory" entry point that takes a topology + trajectory + `family_type` (`kinase` / `gpcr` / `auto`) + UniProt id and runs the same module. In v1 this can be CLI-only (`flexappeal refstates annotate ...`) with the web form stubbed.

### 6.3 Performance targets

Explorer first paint < 2 s; frame scrub < 100 ms per frame (precomputed per-frame JSON, stride 1 ns for the web, full resolution in downloads); per-system web bundle ≤ 150 MB (gzipped mmCIF of the topology, compressed coordinates via Mol*'s supported trajectory format or XTC at 1 ns stride, metrics JSON).

---

## 7. Data model (SQLite, read-only at serve time)

```sql
CREATE TABLE systems     (system_id TEXT PRIMARY KEY, family TEXT, uniprot TEXT, name TEXT,
                          environment TEXT, n_atoms INTEGER, forcefield TEXT, notes TEXT);
CREATE TABLE references  (ref_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
                          pdb TEXT, chain TEXT, source TEXT, source_state TEXT,
                          our_state TEXT, resolution REAL, metrics_json TEXT);
CREATE TABLE motifs      (motif_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
                          name TEXT, token TEXT, residues_json TEXT, generic_labels_json TEXT,
                          source TEXT);
CREATE TABLE runs        (run_id TEXT PRIMARY KEY, system_id TEXT REFERENCES systems,
                          start_state TEXT, start_model TEXT, start_model_rung INTEGER,
                          arm TEXT, replica INTEGER, tier INTEGER, length_ns REAL,
                          platform TEXT, ns_per_day REAL, status TEXT);
CREATE TABLE motif_frames(run_id TEXT REFERENCES runs, frame INTEGER, time_ns REAL,
                          motif_id TEXT REFERENCES motifs, metric TEXT, value REAL,
                          s_score REAL, label TEXT,
                          PRIMARY KEY (run_id, frame, motif_id, metric));
CREATE TABLE membrane    (run_id TEXT REFERENCES runs, time_ns REAL, apl REAL,
                          thickness REAL, tilt REAL, core_waters INTEGER);
CREATE INDEX ix_mf_run ON motif_frames(run_id);
CREATE INDEX ix_mf_motif ON motif_frames(motif_id);
CREATE INDEX ix_ref_sys ON references(system_id);
CREATE INDEX ix_runs_sys ON runs(system_id);
```

Per-residue DSSP per frame lives in Parquet (`data/dssp/<run_id>.parquet`) and is served as compact per-run JSON, not in SQLite.

---

## 8. Phase plan

| Phase | Outputs | Acceptance criteria |
|---|---|---|
| **0 Setup** | `hwprobe.py`, `hardware.toml`, `tuning.toml`, pixi env, `.claude/settings.json`, `CLAUDE.md`, BoltzMaker vendored at a recorded SHA | `pixi run python -c "import openmm, mdtraj, MDAnalysis, gemmi"` passes; OpenMM platform benchmark recorded; BoltzMaker's own `sse_comparison` tests pass from the vendored copy |
| **1 References** | KLIFS and GPCRdb fetchers (through BoltzMaker's cached client), reference ensembles, prepared reference chains, per-reference metrics | ≥ 10 active-like and ≥ 10 inactive ABL1 references, ≥ 5 active and ≥ 5 inactive β2AR references; positive controls P1–P3 pass (Section 9) |
| **2 Start models** | ABL1 seeds; β2AR refined models (or logged fallback rung); OPM-oriented receptor | Both states per system exist, numbering reconciled, every modelled residue flagged; orientation check passes |
| **3 System build + smoke** | Solvated ABL1, POPC-embedded β2AR, Tier 1 runs (1 ns each row) | All Tier 1 runs finish; energy drift within limit; membrane APL within range by end of equilibration; ns/day recorded; Tier 2 projection written |
| **4 Analysis module** | `flexappeal/refstates/` (trajectory wrapper over `sse_comparison`, GPCR micro-switches, state score, DSSP, membrane QC), run on Tier 1 data | Reference structures run through the trajectory code path reproduce their own labels (P1–P3); unit tests green |
| **5 Production** | Tier 2 runs (priority order), streaming analysis | Every Tier 2 row complete or logged as degraded; analysis backfilled; P4–P6 evaluated |
| **6 App** | Examples tab, explorer, reference browser, About, CLI annotate entry point | Playwright screenshots at 1440×900 and 390×844 with no console errors; axe-core zero serious/critical; performance targets met |
| **7 Findings + deliverables** | `FINDINGS.md`, README, icons (blog mark via `marcs-vibe-icon` + favicon set), `deploy/` written not run, `deploy/runpod/` written not run | QC checklist complete; `GATE_OPEN.md` names G3 (and G1/G2 if Tier 3 is wanted) |

---

## 9. Positive controls (resolved at run time from real structures)

| ID | Control | Pass |
|---|---|---|
| P1 | Our DFG spatial label vs KLIFS `DFG` label across all ABL1 references | ≥ 90 % agreement (excluding KLIFS `na`) |
| P2 | Our αC label vs KLIFS `aC_helix` across ABL1 references | ≥ 85 % agreement |
| P3 | GPCRdb `state` reproduced by the TM3–TM6 (3.50–6.34) distance with a threshold fitted on the reference set (leave-one-out) | ≥ 90 % of Active vs Inactive entries correct; report the fitted threshold |
| P4 | Each start structure scores as its own state at t = 0 | whole-protein `s` > 0.3 for active starts, < −0.3 for inactive starts |
| P5 | Inactive-start replicas stay inactive (no spontaneous activation expected on 100 ns) | ≥ 2 of 3 replicas keep whole-protein `s` < 0 over the final 20 ns |
| P6 | Membrane integrity | APL and thickness within Section 5.8 ranges for every β2AR replica |

A P1–P3 or P6 failure after the fallback ladder is **G6**. P4–P5 failures are reported in `FINDINGS.md` as results, not hidden; note that the expected β2AR active-to-inactive relaxation (Dror 2011) occurred on the µs scale, so its absence at 100 ns is an expected limitation, stated plainly, not a failure.

---

## 10. QC

- **pytest:** numbering maps (KLIFS position ↔ residue, generic number ↔ residue) against BoltzMaker's recorded fixtures; metric functions on synthetic geometries with known answers (a hand-built salt bridge at 3.5 Å, a helix rotated by a known angle); state score symmetry (`s(active ref) = −s(inactive ref)` for a two-reference toy); membrane builder smoke test on a 5 ps run; resume-from-checkpoint test.
- **Science sanity:** P1–P6; reference distributions plotted and eyeballed in `FINDINGS.md`; any metric with all-NaN for a motif is a test failure.
- **Playwright:** every page at 1440×900 and 390×844 (preinstalled Chromium, never `playwright install`); no console errors.
- **axe-core:** zero serious or critical violations.
- **Performance:** the targets in 6.3.
- **Units everywhere:** Å, ns, degrees, Å²; threshold lines labelled on every plot.

---

## 11. Deliverables

- `flexappeal/refstates/` module + CLI (`flexappeal refstates build|run|annotate|report`).
- Two example systems with Tier 2 data, explorer pages and the reference browser.
- `FINDINGS.md`: per system, what each motif did from each start, with figures, positive-control table, limitations (timescale, POPC-only membrane, no palmitoylation, apo-only core runs).
- README to Marc's house standard (via `readme_forge.py`), not pushed.
- Icons: blog mark via `marcs-vibe-icon`; favicon set per the vibe-coding skill (if FlexAppeal already has both, leave them).
- `deploy/` (droplet rsync + nginx snippet) and `deploy/runpod/`, written but not executed.
- `DECISIONS.md`, `BUILD_LOG.md`, `PROGRESS.md`, `GATE_OPEN.md`.

---

## 12. Appendices

### 12.1 `.claude/settings.json`

```json
{
  "permissions": {
    "allow": [
      "Edit", "Write", "Read", "Glob", "Grep", "WebFetch", "WebSearch",
      "Bash(pixi:*)", "Bash(uv:*)", "Bash(python:*)", "Bash(python3:*)", "Bash(pytest:*)",
      "Bash(caffeinate:*)", "Bash(git init:*)", "Bash(git add:*)", "Bash(git commit:*)",
      "Bash(git status:*)", "Bash(git diff:*)", "Bash(git log:*)", "Bash(git clone:*)",
      "Bash(curl:*)", "Bash(wget:*)", "Bash(brew install:*)", "Bash(brew list:*)",
      "Bash(sysctl:*)", "Bash(system_profiler:*)", "Bash(df:*)", "Bash(osascript:*)",
      "Bash(tar:*)", "Bash(gzip:*)", "Bash(gunzip:*)", "Bash(mkdir:*)", "Bash(ls:*)",
      "Bash(mv:*)", "Bash(cp:*)", "Bash(gunicorn:*)", "Bash(flask:*)",
      "Bash(sqlite3:*)", "Bash(duckdb:*)", "Bash(npx playwright test:*)",
      "Bash(mkdssp:*)", "Bash(packmol-memgen:*)", "Bash(tmux:*)"
    ],
    "deny": [
      "Bash(git push:*)", "Bash(gh repo create:*)", "Bash(ssh:*)", "Bash(scp:*)",
      "Bash(rsync:*)", "Bash(sudo:*)", "Bash(rm -rf /:*)", "Bash(rm -rf ~:*)",
      "Bash(playwright install:*)"
    ]
  }
}
```

### 12.2 Project `CLAUDE.md`

```markdown
# FlexAppeal: reference-state build
- FLEXAPPEAL_REFSTATES_BUILD_SPEC.md is authoritative.
- Never stop except at Gates G1–G6; follow the gate protocol.
- Log to BUILD_LOG.md, PROGRESS.md, DECISIONS.md; every stage writes data/manifests/<stage>.jsonl and is resumable.
- Commit at the end of every phase; never push.
- Motif annotation comes from vendor/boltzmaker/sse_comparison (pinned SHA); wrap, never edit.
- No hand-typed residue numbers: KLIFS positions and GPCRdb generic numbers only.
- Hardware settings come from config/tuning.toml; OpenMM platform from the Phase 0 benchmark.
- UI uses FlexAppeal's existing design system, else the marcdeller.com brand.
- House style: British English, no em dashes, no banned words.
```

### 12.3 References (verify every DOI via Crossref at build time; drop any that fail)

- Kanev G.K. et al. KLIFS: an overhaul after the first 5 years of supporting kinase research. *Nucleic Acids Res.* 2021.
- GPCRdb: the most recent NAR database-issue paper (Taracena Herrera et al., 2025, as cited in BoltzMaker).
- Ballesteros J.A. & Weinstein H. Integrated methods for the construction of three-dimensional models of G protein-coupled receptors. *Methods Neurosci.* 1995.
- Modi V. & Dunbrack R.L. Defining a new nomenclature for the structures of active and inactive kinases. *PNAS* 2019.
- Shan Y. et al. A conserved protonation-dependent switch controls drug binding in the Abl kinase. *PNAS* 2009 (verify venue: the DFG-Asp protonation result; correct the citation to whatever Crossref returns).
- Dror R.O. et al. Activation mechanism of the β2-adrenergic receptor. *PNAS* 2011.
- Rasmussen S.G.F. et al. Structure of a nanobody-stabilized active state of the β2 adrenoceptor. *Nature* 2011 (3P0G).
- Cherezov V. et al. High-resolution crystal structure of an engineered human β2-adrenergic G protein-coupled receptor. *Science* 2007 (2RH1).
- Eastman P. et al. OpenMM 8. *J. Phys. Chem. B* 2024.
- McGibbon R.T. et al. MDTraj. *Biophys. J.* 2015.
- Dickson C.J. et al. Lipid21 / Gould I.R. et al. Lipid17 (cite whichever lipid parameter set is used, verified).
- Lomize M.A. et al. OPM database (cite the current OPM paper used for orientation).
- Deller M.C. BoltzMaker (github.com/bellcheddar/BoltzMaker), `sse_comparison` at commit `<SHA>`.
