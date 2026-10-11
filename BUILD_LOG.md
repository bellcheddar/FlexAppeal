# Build log

Timestamped, one line per meaningful event, for the reference-state build
(`FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`, section 0.1 rule 5).

## 2026-10-10

- `13:5x` Spec received and read in full (447 lines). Project root is the existing
  FlexAppeal repository, not a fresh checkout (D001).
- `13:5x` Hardware probed: Apple M2 Ultra, 16 performance cores, 8 efficiency
  cores, 128 GB. Project volume 1.3 TiB free, so Gate G5 is not open.
- `13:5x` BoltzMaker cloned to `vendor/boltzmaker/` and pinned at
  `da939b1be54a27cb8eac7c937fd28e8d1e7d008d` (D003).
- `14:0x` Dependency audit against spec section 3: eleven of seventeen already
  present. Added mdanalysis, requests, tenacity, pyarrow, duckdb (D004).
- `14:0x` Vendored `sse_comparison` would not import: its modules pull
  BoltzMaker's main module, which imports openpyxl and plotly at module scope.
  Both added rather than extracting `_align_positions`, which the spec forbids
  reimplementing (D004).
- `14:1x` BoltzMaker's pins (gemmi 0.6.5, scipy 1.13, biopython 1.84) deliberately
  not adopted: this environment runs newer versions that FlexAppeal's own tests
  are built against (D004).
- `14:1x` Four vendored tests failed on missing `boltz_cif` structures, which are
  in BoltzMaker's own `.gitignore` and cannot travel in a clone. Copied from the
  local working copy (D005).
- `14:1x` **Phase 0 acceptance met**: `import openmm, mdtraj, MDAnalysis, gemmi`
  passes; vendored `sse_comparison` suite 50/50 green; FlexAppeal's own suite 759
  passed after each environment change.
- `14:3x` Run matrix recounted on a question from Marc: Tier 2 is **fifteen**
  runs, not twelve. Section 5.3's first row is three replicas per protonation
  arm, so ABL1 active-like is six runs by itself. 1,500 ns in total.
- `14:3x` First compute projection, from this morning's measured 298 ns/day at
  25,463 atoms: 322 GPU-hours, 13.4 days continuous, 3.4x the 96 hour ceiling.
  Not treated as opening G1, because section 2.2 requires a measured budget and
  phase 3 Tier 1 provides one.
- `14:4x` Marc: drop the replicas. Tier 2 becomes five runs of 100 ns, 500 ns in
  total, projected 107 GPU-hours against the 96 hour ceiling (D006). P5 becomes
  a single observation rather than a 2-of-3 majority; recorded, not hidden.
- `15:0x` KLIFS and GPCRdb both reachable. ABL1 resolves to kinase_ID 392, which
  matches the swagger as section 4.2 requires asserting rather than hard-coding;
  158 ABL1 structures returned. GPCRdb returns 413 residues for adrb2_human and
  135 structures, 32 inactive and 84 active at 3.6 A or better.
- `15:0x` Seed check: 2GQG is KLIFS DFG-in/aC-in and 2HYY is DFG-out/aC-in, both
  2.4 A, both as section 4.2 predicts. Ligands 1N1 (dasatinib) and STI
  (imatinib).
- `15:1x` Marc: one ligand-bound and one apo run per system, in opposite states
  (D007). Matrix becomes imatinib/DFG-out plus apo DFG-in for ABL1, carazolol/
  inactive plus apo active for b2AR, with the protonation arm kept as a third
  ABL1 row. Five runs, unchanged at about 107 GPU-hours.
- `15:2x` Marc: two pairs only, and pick the ABL1 ligand most likely to hold
  (D008). Imatinib on DFG-out, because a type II binder occupies the pocket the
  motif would have to move back through, where dasatinib holds the state the
  kinase reaches most readily anyway. Protonation arm dropped. Four runs, 400 ns,
  projected 92 GPU-hours, inside the 96 hour ceiling with four hours to spare.
- `15:3x` `flexappeal/refstates/sources.py` written and fetching live. ABL1
  resolves to kinase_ID 392 and the assertion against the documented id holds.
  Ensembles: 12 active-like, 45 DFG-out inactive, 4 Src-like from 158 structures;
  b2AR 84 active and 32 inactive from 135, no intermediates. All four seed
  structures land in the ensembles they should, which is the first internal
  consistency check that the filters are not nonsense.
- `15:4x` Marc: lock Tier 2 at 50 ns (D009). Four runs, 200 ns, 46 GPU-hours,
  50 hours of headroom under the ceiling.
- `14:3x` ABL1 start structures built. Both trimmed to residues 233-498, the span
  resolved in both 2GQG and 2HYY, so the two trajectories are of the same
  molecule and a motif RMSD between them means something. 262 residues each,
  imatinib kept as a single copy in the bound system.
- `14:4x` Tier 1 smoke on ABL1 apo passed: 46,721 atoms, OpenCL single precision
  after mixed was refused, **220.4 ns/day**. That is better than the 152
  projected, and it puts all four runs at about **33 GPU-hours**, well inside the
  96 hour ceiling. G1 will not open.
- `14:4x` ABL1 apo 50 ns production launched on the GPU.
- `15:0x` **P1 failed at 11/17 and was fixed, not excused.** Six DFG-out
  references came back "other" because the vendored 11.0 A threshold cuts
  through a cluster sitting at 11.1 to 11.3 A. The two classes separate by 8.5 A
  on the difference between the two anchor distances, so the call now reads the
  nearer anchor (D010). **P1 17/17, P2 17/17, both pass.**
- `15:2x` GPCRdb generic numbering wired up and verified against the spec's own
  worked example: `R3.50` resolves to R131, which is what section 5.5 writes.
  Every motif residue the spec names resolves from the table, none typed.
- `15:3x` **P3 passes at 100/100.** TM3-TM6 cleft distance (3.50 to 6.34), with
  the threshold fitted leave-one-out: **11.39 A**, active median 15.85 against
  inactive 8.27, so the cleft opens about 7.6 A on activation. 16 of 116 entries
  were refused rather than measured, because their residue at 131 or 272 is not
  the one the catalogue names: engineered constructs renumber, and a distance
  between the wrong two residues is a number that looks fine.
- `15:3x` **Phase 1 acceptance met**: P1 17/17, P2 17/17, P3 100/100, and both
  ensembles above the spec's minimums.
- `15:4x` ABL1 apo production started: 50 ns, 12.5 M steps, OpenCL single.
- `15:5x` b2AR start models resolved on the top rung of the ladder. 2RH1 has no
  refined model, as the spec predicted; 4LDE (active, 2.79 A) and 5D5A
  (inactive, 2.48 A) do, and 5D5A carries (S)-carazolol already in place, so
  D007's planned transplant is unnecessary. URL pattern logged in D011.
- `15:5x` Both starts oriented onto OPM 2RH1 (bilayer at z = +/-15.9 A) and both
  pass section 5.1's validation. Fit RMSD 0.20 A for the inactive model and
  2.51 A for the active, which is the activation difference showing up before
  any simulation runs.
- `16:0x` **Both β2AR membrane systems build.** Three faults in the way, and the
  error message named none of them.
  - The deposited CAU has no hydrogens, so the SMIRNOFF template could not match
    the residue in the topology. Fixed the way bundles already do it, by taking
    chemistry from the CCD definition and coordinates from the deposition.
  - `addMembrane` believed gemmi's placeholder `CRYST1` of 1 1 1 and sized the
    box height from the lipid patch rather than from the 9.4 nm receptor, which
    wrapped it onto itself (D014). Not a precision problem, which was the first
    guess: it failed identically on the CPU at double precision.
  - The built system's contents exceed its own periodic box. Whole lipids are
    kept when the patch is tiled, so the bilayer overhangs by about 9 Å a side,
    and `addSolvent` then fills the box interior where those atoms already sit
    once the boundary applies. 235 overlapping pairs, 231 water against water,
    with water oxygens carrying 200,000 kJ/mol/nm. Repaired by deleting 246
    molecules of 42,022, which took the as-built energy from -2,056,143 to
    -2,213,385 kJ/mol.
- `16:1x` **The actual blocker was the lipid z restraint**, and it took seven
  diagnostics to get there because four plausible hypotheses were wrong first.
  Not precision, not minimisation depth, not the periodic overlaps, and not the
  protein restraint. Section 5.2 asks for lipid phosphorus held in z, and that
  was written as a plain `(z - z0)^2`. OpenMM wraps coordinates into the box;
  the lower leaflet sits at negative z and wraps by one box vector, so the
  restraint saw a 11.77 nm deviation on 238 atoms that had not moved. At 4184
  kJ/mol/nm^2 that is 138,000,000 kJ/mol, which is the +136,000,000 that had
  been showing up since the first membrane attempt. Written with
  `periodicdistance(0, 0, z, 0, 0, z0)` the same system reads -2,213,385 with
  both restraints at exactly zero.
- `16:1x` A reasoning error worth recording with it. The +136,000,000 was
  attributed for several hours to single precision, on the strength of a CPU
  run reading -2,060,000 for the same system. Those two runs did not differ in
  precision: one had no restraints and the other had both. Nothing was ever
  measured that showed the platforms disagreeing, and after the repair they
  agree exactly at -2,213,385.
- `16:1x` `periodicdistance` is the right form for both restraints rather than a
  plain harmonic, for the same wrapping reason: written plainly the protein
  restraint reported 5,045,427,552 kJ/mol on atoms sitting exactly at their
  reference positions, where `periodicdistance` reported 18.
- `16:0x` Repair added to the build rather than worked around in the protocol:
  solvent molecules overlapping another molecule's periodic image are deleted,
  water before ions so the neutralisation `addMembrane` chose is preserved, and
  a clash between two non-solvent molecules is reported and left alone, because
  deleting a lipid or part of the receptor is not a repair.
- `16:0x` Phase 4 analysis complete and verified on the partial ABL1 trajectory:
  43 frames scored, all six motifs, every motif RMSD between 0.16 and 4.11 Å.
  Two of the three bugs found along the way produced confident wrong answers
  rather than errors (FINDINGS section 3).
- `16:0x` Phase 6: a fifth tab reading one SQLite file, with the per-run
  explorer, the generated finding and the reference browser where P1 and P2 are
  visible to a reader as 17 of 17 with three unlabelled.
- `16:39` **Tier 1 membrane run reaches production.** The whole of section 5.2's
  membrane column executed in order for the first time: minimise with the
  protein restrained and lipid P held in z, six gentle heating stages at 1 fs
  from 50.5 to 303.1 K, the remaining 0.4 ns of NVT at 4 fs, barostat on, the
  six-step restraint ladder with the lipid z restraint released after the
  second, then unrestrained production. Minimisation took 55 s against about 5
  minutes before the clash repair, because the system no longer starts from a
  pathological energy.
- `16:3x` **P4 passes 4 of 4.** Each start structure scores as its own state
  with no simulation in between: +0.305, -0.368, +0.593, -0.934. The kinase
  rows clear the band by under 0.07 and the reason is interpretable, which
  FINDINGS records rather than leaving as a near miss.
- `16:3x` Two false claims found and fixed by running the controls rather than
  by reading the code. TM6 reported `conflict` on every inactive frame because
  the ionic lock was attached as its corroborating rule, and that lock reads
  broken in both beta2AR states. And the reference browser reported 116 of 116
  agreement for the receptor where `our_state` was a copy of the catalogue's
  own; it now comes from the TM3-TM6 cleft against P3's threshold and reads
  100 of 100, which matches P3 exactly.
- `16:4x` Section 12.3 satisfied: 15 of 15 citations resolved through Crossref,
  after four intermediate versions that each raised the pass count while making
  the output wrong. Supporting-information deposits and PDB data entries carry
  the paper's exact title with their own DOI, and a meeting abstract carries it
  in a different journal a year later.
- `16:4x` Viewer bundle at a 1 ns stride: 0.4 MB for the partial ABL1 run,
  about 1 MB for a finished row, against section 6.3's 150 MB ceiling.
- `21:0x` Both reference manifests now reproduce from source.
  `tools/build_references.py --verify` rebuilds them in memory and diffs per
  record: ABL1 20 of 20 and beta2AR 116 of 116 with nothing changed, P1 17/17,
  P2 17/17, P3 100 of 100 at 11.39 A, and `GPCR_CLEFT_THRESHOLD_A` confirmed
  against the fresh fit. Until this ran, four of the spec's acceptance numbers
  rested on artefacts no code in the repository could regenerate.
- `21:0x` The first kinase rebuild disagreed and the script was wrong, not the
  manifest: 42 rows against 20, because it deduplicated on `(pdb, chain)` where
  the original deduplicated on the deposition. KLIFS lists a row per chain, so
  ABL1 is 61 rows over 20 depositions. Per-PDB first occurrence reproduces all
  20 committed chains exactly, including the six that are not chain A. Worth
  noting that the wrong version still *passed*, at P1 37/37 and P2 36/37, so
  passing and reproducing are not the same check (D016).
- `21:1x` The per-chain set is kept as `--all-chains` and found one disagreement
  in 37: 2HZ4 chain A reads the catalytic Lys to alphaC Glu at 4.20 A against a
  4.0 A cutoff and calls alphaC `out`, where chains B and C read 2.71 and 2.74
  and call it `in`. A 0.2 A threshold miss in one copy of a three-chain crystal.
  P2 would still pass at 16/17 had that copy been listed first, which is the
  useful part: the control does not depend on catalogue row order.
- `21:1x` Asked which other committed constants are fitted rather than chosen,
  and found one more with the same problem: D015's 7.7 A R-spine threshold.
  The kinase verify now re-derives it. The conclusion holds, two of D015's
  figures were wrong, and the correction strengthens it: measured on the pair
  that actually decides the call, the published 6.0 A clears the active cluster
  by 0.04 A rather than by a tenth of an angstrom. It also surfaced 2G1T at
  6.14 A, the one reference whose label flips between the two thresholds, which
  D015 had left out of its count.
- `21:1x` Deleted `ensembles.common_core`, `traj.kinase_core` and `sse.simplify`
  with its `SIMPLIFIED` table, the three functions from the audit that were
  genuinely unreferenced rather than a missing caller.
- `21:1x` Corrected `PROGRESS.md`, which still attributed the +136,000,000
  kJ/mol to single precision and to the wrong fault. It was the non-periodic
  lipid restraint, and that correction had reached FINDINGS and BUILD_LOG but
  not the file a reader is most likely to open first.
- `21:1x` Went looking for what was still untested with Tier 2 about to produce,
  rather than asserting readiness again. The receptor ingest path had never
  executed: the database held one system, one run and 20 ABL1 references, so the
  116 beta2AR references and the five receptor motifs were implemented, fixture
  tested and never run on a real analysis. Exercised on a scratch database it
  works unchanged: 5 motifs, 100 frame rows, 116 references annotated, 100 of
  100 agreeing with GPCRdb.
- `21:1x` The bilayer QC had also never written a row. It runs at ingest rather
  than analysis time and needs `solvated.pdb` beside the trajectory, which the
  smoke run does not have because `--system-dir` put the system files elsewhere.
  Run on the smoke trajectory with the production topology: 476 lipids at
  238/238, thickness 37.62 A, area per lipid 68.9 A^2, 10.1 core waters, G6
  closed. The production directory does have its topology, so this will not
  return silently when the row lands.
- `21:1x` Misread `outside expected` beside a closed G6 as a contradiction.
  There are two windows by design, 60 to 68 A^2 from the literature and 55 to 75
  for the gate, and 68.9 falling between them is the case they exist to tell
  apart. The design was right and the reading was mine.
- `21:2x` A false alarm worth keeping. `pytest -q` from the repository root
  reported 17 failures, all in `vendor/boltzmaker/tests`, which shell out to an
  interpreter that does not exist under pixi. `pixi run test` was already scoped
  to `tests/`; `pytest.ini` now scopes the bare command too.
- `21:2x` Two process faults behind that. The run before it piped pytest into
  `tail`, so the exit status I read was `tail`'s: 0 while pytest returned 1, and
  I reported a suite as passing that had 17 failures in it. Never read an exit
  code through a pipe. And the test count in `GATE_OPEN.md` had been wrong
  twice, 833 then 837, written from memory and then from arithmetic. Measured it
  is 871.
- `21:2x` The streaming analysis loop could never stop. Its guard was
  `while pgrep -f "tools/run_md.py"`, and the script's own text contains that
  pattern, so the shell that launched it matched: the loop was held open by its
  own parent. Verified by listing every matching PID, one of which was the
  loop's grandparent. The streaming passes would have gone on working, and only
  the final pass after the loop, the one that computes the bilayer QC and G6 on
  the finished trajectory, would never have run. Replaced with the queue's own
  `queue done` marker plus a process check that excludes this script's
  ancestors by PID, both paths tested before the swap.
- `21:2x` Restarting that loop populated the receptor half of the database for
  the first time: 2 systems and 136 reference structures, against 1 and 20
  before. The receptor ingest path is no longer merely implemented.
- `21:2x` b2ar_apo production started at 21:22:36, 50 ns over 12,500,000 steps,
  within a minute of the 21:22 the restraint ladder projected.
- `21:2x` Extracted `spine_verdict` from the R-spine control so the half that
  decides pass or fail can be tested on numbers rather than on 20 structures
  and a catalogue. A verdict that inverted would otherwise still print a tidy
  table. Both manifests verify clean in one run: 20 of 20, 116 of 116.
- `21:2x` P6 evaluated on a real bilayer rather than on fixtures: the smoke
  run's chunk through ingest into a scratch database, then the controls command,
  gives `APL 68.9 A^2, 37.6 A  PASS`. The chunk reads `outside expected` and
  still passes because section 5.8 names two different windows, an advisory 60
  to 68 and a G6 gate at 55 to 75, and P6 is the gate. P5 is now the only
  control still genuinely waiting, and it needs an inactive-start row.
- `21:3x` Re-checked every explorer route against the newly two-system
  database, because until the streaming loop was fixed the pages had only ever
  been served a kinase-only one. All eight serve: the receptor reference
  browser renders its 116 rows at 40 kB, both state maps carry their axes and
  the trajectory path, both runs serve frames, DSSP and a viewer bundle, and the
  refusals still refuse (an unknown run, an unknown system and a `viewer.exe`
  all 404). Two of my own URLs were wrong in the first pass, not the routes.
- `21:4x` A resumed run had no barostat, and `abl1_apo` is queued to resume from
  4.0 ns. The force was added inside the equilibration branches, which a resume
  skips. Proved by loading the real checkpoint into a freshly built context:
  the load succeeds, reports `MonteCarloPressure` and has no barostat among its
  five forces, so 46 of abl1_apo's 50 ns would have been NVT with only a
  flat density to show it. Both call sites now share one `barostat()` helper
  and the resume adds it before the context is built (D017).
- `21:4x` Two smaller faults in the same place. The membrane ladder saves a
  checkpoint after every rung, so a kill part way down it left one that a
  resume would have treated as production: restrained, with the equilibration
  time counted against the 50 ns target. And the water branch saved before
  resetting the clock, with the same effect in a narrower window.
  `checkpoint_phase.json` now records the phase, absent meaning production.
- `21:4x` Verified after the change on the real abl1_apo checkpoint: six forces
  including MonteCarloBarostat, loaded at 4.0000 ns, barostat present. The
  running b2ar_apo is unaffected, since its process holds the old code.
- `21:4x` Made three silent-empty paths audible in ingest: a membrane row whose
  topology is not beside its trajectory, a membrane row whose topology has no
  lipid phosphorus, and a missing manifest. Each now says so, and a water row
  stays quiet so the warning keeps its meaning. All four cases tested.
- `21:5x` `abl1_imatinib` would have crashed at build time. The queue had no
  `ligand` field, its start model contains imatinib, and amber14 has no template
  for a drug: the row is fourth, so this was about twenty hours away. Passing
  the CCD SDF is not enough and fails as "similar to CHL1", naming cholesterol,
  because the generator matches a chemical graph and a PDB ligand's bonds are
  inferred from distance. `run_md.py` now imports `ligand_molecule` from
  `build_membrane.py` rather than keeping a second copy (D019). System built in
  place: 48,185 atoms, formal charge zero.
- `21:5x` A guard I added was wrong in a useful way: `getUnmatchedResidues`
  reports against static templates only and ignores a registered generator, so
  the refusal rejected a working setup. It now runs only when nothing is
  generating templates, which is the real hazard.
- `21:5x` The test written to stop that recurring found a second undeclared
  residue, and a more important one. `abl1_apo_start.pdb` carries PTR 393,
  activation-loop phosphotyrosine, and PDBFixer rewrites it to TYR: both ABL1
  rows simulate a dephosphorylated kinase domain started from a phosphorylated
  structure, and nothing said so. Kept, since amber14 ships no phosphotyrosine
  and the choice is the usual one, but now logged by both builders and written
  up as a limitation (D018, FINDINGS 4j). Tyr393 phosphorylation is one of the
  main things that holds ABL1 active, so this is not a packaging detail.
- `21:5x` Noted while checking it: the built systems are renumbered 1 to 262
  where the start model is 233 to 498, and the analysis restores the deposited
  numbering through `traj.deposited_numbering` and `traj.renumber`. Verified
  rather than assumed, since every motif in this package is addressed by its
  deposited number.
- `22:0x` A bug of my own, caught by writing the test for the previous one: I
  changed `--ligand` to take a CCD code and left the queue passing a path to the
  SDF. Nothing would have caught it until a row had to build a system, because
  a pre-built one is reused and the argument is never read. Fixed, and the test
  now asserts both sides mean the same thing.
- `22:0x` All four rows verified buildable and self-consistent: b2ar_carazolol
  187,848 particles with lipids present, abl1_imatinib 48,185 without, both
  with their system particle count matching their topology.
- `22:1x` The bilayer QC loaded the whole solvated trajectory. At a finished
  50 ns membrane row that is 500 frames of 194,032 atoms, about 1.2 GB held at
  once for data used one 10 ns window at a time. It needs the solvent so it
  cannot drop atoms the way `traj.load` does, but the chunks were already
  independent, so it now streams with `iterload` and holds one chunk. Verified
  to give identical numbers on the smoke trajectory, and the cross-chunk
  bookkeeping checked at five chunks: contiguous windows, every frame counted
  once. A test builds its own two-plane fixture rather than committing one.
- `22:1x` Corrected the About page, which still said the published R-spine
  threshold clears the active range by "a tenth of an angstrom". It is 0.04 A,
  and the page now also names 2G1T as the reference whose label turns on the
  choice.
- `22:1x` Added a "What these runs are not" section to the About page. The
  dephosphorylation was written up in FINDINGS and DECISIONS but not on the
  page a reader of the results actually sees, along with the single replica,
  the 50 ns against a microsecond relaxation, and the scores being relative to
  a reference set.
- `22:2x` Armed `tools/watch_stall.sh`. The existing watchers cover completion
  and failure, both of which the queue writes to its log; neither covers a
  hang, where every process stays alive and every log goes quiet. Over a
  19-hour row that is the failure most worth catching. Stall is defined on the
  state data, not the trajectory file, whose size can sit unchanged while a
  write is buffered, and the patience is 20 minutes against a 2.4-minute frame
  interval. It follows whichever row the queue says it is on rather than naming
  one. Both branches tested.
- `22:2x` G4 opened. Marc authorised the push and 82 commits went to
  `origin/main`, the whole reference-states build. Checked first, because the
  remote is public: nothing private tracked, no structure string in any tracked
  file, and every directory holding coordinates, trajectories, caches, built
  systems, manifests or the explorer database is ignored. Five files under
  `data/` are tracked and none holds a coordinate.
- `22:2x` That check found a correction worth more than the push. The reference
  manifests are gitignored and were never committed, and my own write-up had
  called them committed. It makes the gap `tools/build_references.py` closes
  worse rather than smaller: a fresh clone has no manifests at all, so before
  that script the reference browser, both state maps and P3's threshold could
  not be reproduced on any other machine. Corrected in FINDINGS 4g and in the
  script's own wording, which now says "on disk" rather than "committed".
- `22:2x` Measured the analysis and the explorer at the size a finished row will
  be, since every performance figure so far came from a partial run. Tiled the
  receptor trajectory to 500 frames: analysis 19.8 s, viewer bundle 50 frames at
  1.4 MB against a 150 MB ceiling, analysis JSON 1.9 MB, and every endpoint
  between 1 and 34 ms with the largest payload 161 KB. Nothing to do.
- `22:2x` Nearly recorded a number I had not measured. `frames.json` came back
  at 12 KB and 8 ms, which I took for the 500-frame figure; it reads the SQLite
  table and the database still held the real run's ten frames. Ingesting the
  tiled analysis first gives 97.5 KB and 17.1 ms. Pointing a test at a large
  input is not the same as the code under test reading it.
- `22:5x` G-REUSE, one of section 1's six goals, was annotating every kinase
  against ABL1. `annotate` accepted `--protein` and then called
  `kinase_setup(start, reference_pdb)` without it, so CDK2 resolved KLIFS
  kinase 392 with 43 structures, which is ABL1, and reported a clean
  `intermediate` with no conflicts. The receptor half had the same hole against
  beta2AR. Never caught because the hook had never once been run on a third
  protein.
- `23:0x` Fixing it exposed the second half: with the lookup corrected, CDK2
  failed on `--reference-pdb`, which defaulted to the ABL1 entry 2GQG and
  raised a bare `KeyError: '2GQG'` out of a dictionary of CDK2's structures. It
  now picks one of the kinase's own and says which, and refuses with a message
  naming valid codes when given one that is not that kinase's.
- `23:1x` Then the ensembles were uncapped, so CDK2 tried to fetch all 337 of
  its active structures and was still going after fourteen minutes. Capped at 40
  per state for the reuse path, best first by quality or resolution; the
  production path stays uncapped so every committed number is unchanged (D020).
  The same run now takes 90 seconds.
- `23:1x` And of those 40, the core alignment kept 1. A score resting on one
  reference reads exactly like one resting on forty-five, so a state below five
  usable references now writes a note into the analysis payload and the log
  (D021). CDK2's run says so in its own output.
- `23:1x` Verified the production path is untouched: P4 still 4 of 4 at +0.305,
  -0.368, +0.593 and -0.934, ABL1 still 12 active and 45 inactive.
- `23:2x` Destroyed and restored the acceptance table in `GATE_OPEN.md`. It
  lived inside the G4 section, and rewriting that section to record the push
  deleted the whole thing: P1 to P6, the test count, the axe-core and Playwright
  results and every budget. Caught because a later edit asserted on a line that
  was no longer there, which is the only reason it did not go unnoticed.
  Restored above the gates as its own section, so a gate rewrite cannot take the
  results with it.
- `23:2x` Chasing the CDK2 alignment drop led back to the worked example and to
  a number I had been quoting wrongly all build. ABL1's "12 active and 45
  inactive references" are catalogue rows: KLIFS lists one per chain and the
  ensemble is keyed on the deposition, so the scores rest on **5** active and
  **14** inactive structures. Nothing is lost on the way, since the builder
  records a reason for every row it does not keep and recorded none; the rows
  collapse rather than drop, which is the same choice D016 made for the
  manifest. What was wrong was the reporting, which counted the flattering
  number. Both counts are now logged, D020 is corrected, and the About page
  says it too, because an active score being the distance to the nearest of
  five crystal structures is something a reader needs.
- `00:3x` Fourteen streaming passes had printed `membrane 0 rows` on a POPC run.
  That part was deliberate: the streaming loop passes `--skip-membrane`. What it
  meant was that the final pass had never once run on this trajectory, so I
  rehearsed it against a scratch database and found three faults waiting for the
  finished row. The tilt column was a dash and always would have been, because
  `analyse_chunk` takes `tm_ca_indices` optionally and nothing ever passed it.
  With it passed the first number was 25 degrees, for a receptor built square in
  the membrane: the measurement fitted a principal axis to a Calpha cloud that
  is 40 by 45 by 37 angstrom, ratio 1.15, so it was fitting noise. Section 4.3
  already had the right definition, the mean of the seven TM helix axes, which
  gives 6.7 degrees as built and 13.5 over the first 7.6 ns. The old function
  was correct for one helix, which is what both its tests fed it, so it was
  split in two rather than deleted. And the panel showed five numbers with no
  verdict against two windows, which is not a panel a user can read the
  bilayer's health off.
- `00:4x` `python -m flexappeal.refstates.cli` had never been run. It dies on
  `NameError: name 'controls' is not defined` for any argument, because the
  `__main__` guard sat above the subcommand functions. 896 tests never saw it,
  since importing the module defines everything first, and nor did the
  `FlexAppeal.py refstates` wrapper, which imports `main` rather than running
  the script. The only broken path was the one a user would take. `build` was
  missing from the surface too.
- `00:5x` FINDINGS had no figures, which section 11 asks for by name, and no
  limitations section, which it also names. `tools/make_figures.py` now draws
  them from the explorer database, going through the page's own map endpoint so
  a figure cannot disagree with the page. The kinase map paid for itself
  immediately: with the published 11.0 A cut drawn on both axes, which is where
  the criterion applies it, D010's six disagreements are visibly sitting a
  tenth of an angstrom the wrong side of the horizontal line. The receptor map
  showed something I had not measured: `b2ar_apo` runs at a TM3-TM6 distance of
  11.17 to 15.32 A against an active cluster at 15.3 to 16.5, and one frame is
  already below the 11.39 A threshold. A single-cut classifier would call that
  frame inactive on a receptor that is active-like on every other motif, which
  is the clearest argument yet for scoring against ensembles rather than lines.
- `01:0x` Asked what the explorer would show if a row died overnight. It would
  show `running`, indefinitely: `run_md.py` writes `run.json` only on success
  and the ingest reads an absent one as a row still in flight. Section 0's
  ladder ends in `status = "failed:<reason>"` and nothing wrote it for a run.
  The queue now writes it, because the queue is what survives an OOM kill.
- `01:0x` Rehearsed the two rows that have never produced a frame, by writing a
  twelve-frame trajectory from each built system into the scratch directory.
  Both analyse first time, `abl1_imatinib` at -0.281 and `b2ar_carazolol` at
  -0.662, and P5 printed a result for the first time in the build. Nothing
  about the physics, everything about whether 24 GPU-hours were about to be
  spent discovering a crash.
- `01:1x` Checked section 6.2's viewer line against the page. Four features
  asked for, one present. The cartoon colouring is absent and the page says so.
  The membrane toggle and the reference overlay were absent with nothing saying
  so. Built the membrane toggle: the reason it was missing was that the viewer
  bundle drops the bilayer on size grounds, and the bilayer phosphorus is 476
  atoms, 38 KB, against 0.6 MB of protein. The size argument was about the
  water. Verified in a browser; the spheres themselves could not be confirmed
  because this Chromium has no WebGL, and that is recorded rather than glossed.
- `01:2x` Section 7's `systems` table has `uniprot`, `n_atoms` and `forcefield`
  and nothing had ever written any of them. The Examples card already had the
  atom count in its markup behind a truth test, so it rendered nothing and read
  as a card that was never meant to have one. Accessions taken from the
  catalogues the build already queries rather than from memory, atom counts
  measured and agreeing with section 5.1's independent figures. Adding the
  force field to `run.json` nearly cost a run: the list is chosen in `build`
  and written in `main`, different scopes, so it would have raised NameError
  after the production loop finished and before anything was recorded. Module
  constant now, with an `ast` test that was checked against the bug.
- `01:3x` Built section 6.2's thumbnail: a Calpha trace from the built system's
  own coordinates, with the 476 lipid phosphorus for the receptor. The first
  one looked fine and was wrong, with the receptor apparently half out of the
  membrane, because the coordinates are wrapped and the bilayer straddles the
  periodic boundary. Centred with `centre_bilayer`, which returns its shift for
  exactly that. A picture can be confidently wrong too, and is harder to check.
- `01:3x` The bilayer toggle used `md-check` without the span the pattern needs.
  The house CSS hides the real input and styles the span, so the control was
  invisible, had no checked state and no focus ring, and read as a line of
  prose that happened to respond to a click. Corrected and then clicked in a
  browser, which is the user path, rather than driven through the DOM.
- `01:4x` Checked the About page against the nine items section 6.2 names for
  the software table. Seven present, the SHA present, the references verifying
  through Crossref. PDBFixer was missing, which matters because it is the tool
  behind D018's dephosphorylation and the page carries that limitation without
  naming what causes it. MDAnalysis is absent because this build does not use
  it: it appears only in the downloadable run bundle's environment, and listing
  a library that is not used would be worse than the gap.
- `01:4x` Checked section 6.2's motif table columns. Five of six: the nearest
  reference PDB at the end was missing, which is the column D023 had just cited
  as the reason not to build the viewer's reference overlay. The justification
  for declining one thing rested on another thing that was also absent, and I
  wrote it without checking. Column added, reading the payload the view already
  opens, both references shown rather than the nearer one. Also recorded that
  the motif timelines are not the chart section 6.2 describes.
- `02:2x` Built section 6.2's reference consensus rows on the DSSP ribbon. The
  data was in the analysis payload the whole time; what was missing was taking
  the modal code and resolving the motif position to a residue number. Drawn as
  two columns left of t = 0 so they share the residue axis exactly. Took three
  goes to be readable: at one frame step the columns were 0.1 ns slivers, and
  the tick labels printed over each other. Width is now a fraction of the run,
  which is constant in pixels at any run length.
- `02:4x` Built section 6.2's state badge strip, server-rendered at frame 0 and
  then following the scrubber. Needed `start_label` in the summary, which was
  not there; rendering the final label instead would have disagreed with the
  viewer from first paint. Then found `STATE_SHAPES` was computed, passed to
  every run render and written into the DOM as `data-shape`, and read by
  nothing: there was no `.md-state` rule in the only stylesheet. Section 6.1's
  "shape or label" was being met by the label alone. The marker is drawn now,
  from the state class, and the redundant map is gone.
- `02:5x` The ribbon's "motif bands marked" clause was unbuilt, and building it
  turned up the colour system behind it. `--motif-1` to `--motif-6` are stored
  on every motif row and were defined nowhere in the only stylesheet, and the
  map assigning them is section 6.1's own table, which does not match the
  motifs the build ended up scoring: `c.l`, `TM5` and `TM7` had no colour,
  while `r_spine`, the ionic lock, `PIF` and `toggle` had one and are measured
  as pocket states, micro-switch metrics and corroborating rules rather than as
  motifs with ensembles of their own. Palette defined in Okabe-Ito, and the
  table honoured wherever it names a motif the build renders, with only the
  free tokens reused. I had also moved helix 8 off its spec token for tidiness
  and put it back: reading the clause rather than working from the shape of the
  problem is what caught that. The ribbon bands, the state timeline and the
  figures now all read the token. Plotly's autorange ignores shapes, so the
  first version clipped every band it drew.
- `03:0x` Section 6.1's table has a last row naming `--state-active` and three
  siblings for "state badges, map colours". Not defined either, and the state
  markers I had written an hour earlier hard-coded the four values into the
  rules: the same fault as the motif tokens, in the same stylesheet, while I
  was writing the motif tokens up as a fault. Defined, referenced, and read by
  the figure generator for the state map. Also re-measured section 6.3's three
  budgets after the explorer grew: first paint 432 ms of 2 s, scrub handler
  0.1 ms of 100 ms, bundle 3.5 MB of 150 MB at 50 ns.
- `02:5x` Audited section 2, which is the hardware section and had never been
  checked. It names five files and none existed: `tools/hwprobe.py`,
  `config/hardware.toml`, `config/tuning.toml`, `config/benchmarks.toml` and
  `data/manifests/md.jsonl`. The settings themselves were right, because
  `run_md.py` derives the platform at run time and does it better than the
  spec's table asks, but nothing recorded what the machine is or what it
  managed, and section 2.2 says every compute budget is derived from these
  measurements rather than assumed. Probe written, three files generated,
  benchmarks refreshed in the final pass as rows finish. The queue now appends
  to md.jsonl on start, finish and failure; the campaign already running
  cannot, since Python does not reload a live module, so the file opens with
  two lines marked reconstructed rather than pretending or staying empty.
- `02:5x` The probe test caught me writing the same `relative_to` crash that
  4o records from make_figures: `--out` outside the project kills the log line.
  Same bug, same shape of function, three hours apart.
- `03:0x` Checked section 3's package table by grepping for every name rather
  than reading D004 and believing it. Three of the five packages D004 added are
  unused: `mdanalysis` for membrane analysis and helix axes, both of which are
  hand-rolled in `membrane.py`, and `pyarrow` and `duckdb` for a per-residue
  DSSP store "too large for JSON" that is 46 KB of JSON. The estimate behind
  two dependencies was never checked against a file. D004 amended, 4z written,
  packages left installed because removing them from a lock that serves the
  deployed app saves nothing.
- `03:1x` Three page tests failed on a loaded machine and all sixteen passed in
  isolation. The failure was the screenshot, not the assertion: the test falls
  back from full-page to viewport capture and the fallback had no fallback, so
  a clean page was reported broken. Capture is best effort now, which is what
  its own comment had already argued for.
- `03:1x` Section 12's two files. `CLAUDE.md` written, with the PROGRESS.md
  cadence and the "never push" line changed to match the user's overrides and
  both changes marked in the file. `.claude/settings.json` left undone and
  recorded in D024: its deny list blocks `git push`, which the user authorised,
  and amending a permissions policy for future sessions is their call.
- `03:1x` Section 10 names five pytest items and four existed. Wrote the fifth,
  the membrane builder smoke on a 5 ps run, against the system the builder
  produced rather than a fixture. It skips with a stated reason while a
  production run holds the GPU, so the body was proved separately on the CPU
  platform with ten steps: 187,848 particles, which matches build.json to the
  atom, -2,353,259 kJ/mol, every coordinate finite. A test that has only ever
  skipped is not a test.
- `03:2x` Ran section 8's phase 0 acceptance criteria as written for the first
  time. The import line passes. BoltzMaker's own `sse_comparison` tests pass
  50 of 50 from the vendored copy, which I had never checked: what I had
  measured was the whole vendored suite, whose 17 failures are elsewhere and
  are about its own interpreter path. `pytest.ini`'s comment implied the
  vendored tests were unusable; sharpened, with the command to repeat it.
- `03:3x` Went after the CDK2 limitation 4l recorded and found it was not what
  4l guessed. The cores are near-uniform: 28 of CDK2's 29 references give a
  42-residue KLIFS scaffold and one gives 41. What broke it is that
  `align_to_common_frame` took its target as the first reference in sorted
  order, which was the odd one, so every other reference failed the shape check
  and the ensemble collapsed to the target alone. The "core residual mean 0.0 A"
  in the log was one structure fitted onto itself and should have been the tell.
  Target is now the first with the modal core length: CDK2 goes 1 of 29 to 28 of
  29, core agreement 0.49 A mean. Production is bit-for-bit unchanged, which is
  the check that mattered: ABL1 19/19 onto 2F4J at 0.568, beta2AR 100 onto 3P0G
  at 0.829, P4 still +0.305, -0.368, +0.593, -0.934.
- `03:4x` `abl1_apo` resumes unattended tonight and the only cover for 4i's
  barostat fix was a source-order test that says so in its own docstring. Eight
  argon particles on the Reference platform turn out to be enough to
  demonstrate the trap: a checkpoint written with a barostat loads without
  complaint into a system without one, and the context still reports
  MonteCarloPressure. The source test stays; the two are different claims.
- `03:5x` First real per-motif production reading, at 14.6 ns. The
  whole-protein score has not moved, -0.022, while H8 falls 0.176 and NPxxY
  rises 0.101 underneath it: the first case in real data where the headline
  number alone would have said nothing happened. Checked before writing it up,
  which was worth doing: H8 split into thirds is 0.328, 0.153, 0.143, a step
  then a plateau with scatter nearly as large as the change, and its helix
  fraction is 0.909 throughout. "H8 is drifting" is not supported; NPxxY's
  approach to 0.42 A from the nearest active reference is the cleaner signal.
- `04:0x` Put the halves and the scatter into `store.summarise` so section 5
  does not need hand computation per row, with the page and the report reading
  the same definition. Two faults while adding it, both caught by the test: the
  count read 294 frames for a 147-frame run, because motif_frames holds one row
  per metric, which is 4m's rows-against-things error made again; and the first
  test asserted on the shared two-frame fixture, which is below the four frames
  summarise needs, so it passed by skipping everything.
- `04:0x` The third fault in that change broke the run page. The halves needed
  a midpoint and I named it `cut`, which already held the final-tenth time
  boundary computed before the motif loop. First motif correct, every one after
  it querying for frames past 73 ns of a 14.6 ns run and getting NULL: nine of
  eleven lost their final-tenth mean and the page 500'd. The suite caught it in
  ten tests as well. The regression test asserts on every motif, because an
  error that starts at the second element is invisible to one that checks the
  first.
- `04:1x` Fixing the shadowing took the suite from ten failures to four, and
  the four were the new column: `row.sd_s is not none` does not catch a missing
  key, because Jinja gives Undefined rather than None and `Undefined is not
  none` is true. `summarise` omits the scatter for a run too short to halve,
  which is section 0's rule working as intended, so the template has to test
  for absence. `row.get("sd_s")` does.
- `04:0x` Three faults in one small change is a signal about the change, so the
  next thing added was a check of a different kind rather than more code: the
  analysis JSON and the SQLite summary are independent paths to the same seven
  numbers and must agree. Against live data all eleven motifs across both rows
  agree on every one, zero disagreements. Now a test, checked by putting the
  shadowing back: it fails with the bug and passes without it, and it would
  have caught the fault on the first run instead of when the page 500'd.
- `04:0x` Made the round-trip check run on real data, not only a fixture:
  `report --verify` compares every stored number against the analysis file the
  run was built from, and the streaming loop's final pass calls it so the
  unattended landing verifies itself. Reported rather than fatal, so a
  disagreement does not cost the figures and the controls. There is a test that
  it can fail as well as pass, because a verifier that cannot fail is worth
  nothing. Live right now: every number agrees, 11 motifs, 2 rows.
- `04:0x` Rehearsed the streaming loop's final pass end to end, with the loop
  stopped so nothing competed. Every piece had been checked alone; the sequence
  had not, and it is the part that runs with nobody watching. It completes, exit
  0. Two firsts: P6 passes on the live database, 2 chunks at 67.7 A^2 and 38.6
  A, so five of section 9's six controls are now decided; and the bilayer figure
  drew, showing area per lipid falling 68.0 to 67.2 into the expected window,
  which is what 5.2 predicted from the Tier 1 value of 68.4.
- `04:0x` The rehearsal found the verification step announcing a failure on
  every success. `||` on a pipeline takes grep's status, and `log` writes to
  stderr so the "agrees" line never reached the pipe. The same construct would
  have swallowed a real disagreement, because grep matching "disagree" exits 0.
  Rewritten to test the command's own status; both branches run rather than
  reasoned about.
- `04:1x` Tried to read abl1_apo's checkpoint to find exactly where tonight's
  resume will restart, and learned something better: checkpoints are bound to
  the platform that wrote them. `platform()` falls back to CPU when OpenCL is
  missing, which is right for starting a row and a trap for resuming one, so a
  GPU that disappears overnight becomes an OpenMM traceback with the row's
  nanoseconds still on disk and nothing saying they are good. Guarded, and
  deliberately not recovered from: restarting would discard the checkpoint's
  progress unasked. The original question stays open, because reading the
  checkpoint needs the GPU production is using.

### 11 October, 04:19 to 04:55

Stopped a background shell I had read as a stale poller. Its command ended in
`exec caffeinate -dimsu pixi run python tools/queue_md.py`, so the shell was the
production queue and `b2ar_apo` went from 17.2 ns of 50 to no process. Restarted
it at 04:20: resumed from the checkpoint at 17.000 ns, barostat re-added,
2 trajectory frames and 2 state rows past the checkpoint dropped, 33.0 ns
remaining, ETA about 17:35. Cost 0.2 ns of sampling and two minutes of idle GPU,
and bought the first unattended run of the resume path on a membrane row.

Two faults found while checking the restart was safe, both fixed:

- `run.json` recorded `production_ns`, the target the row was asked for, and
  `queue_md.complete` compared it against the same target, so a row that
  stopped short certified itself. Now records `produced_ns` from the context and
  judges on that, with `status = "short:40.00_of_50ns"` for a row more than
  50 ps short, and the target as the fallback for rows written before the field
  existed. Verified against the bug: the old expression returns True on a
  40-of-50 row.
- The "already at target, nothing to do" exit returned 0 without writing
  `run.json`, so a row resumed to exactly its target would have read as in
  flight for good. Both exits now go through one writer, which omits
  `ns_per_day` and `wall_hours` rather than zeroing them.

Also found: a row in flight runs the code it started with. `b2ar_apo` had no
`checkpoint_phase.json` because the marker was added at 21:47 and the row
started at 16:07. Harmless, because an absent marker means production, but had
the default been `equil` the restart would have discarded 17 ns and silently
repeated a 10 ns equilibration.

FINDINGS 4ag. Seven tests added, 939 passed, 1 skipped.

### 11 October, 04:28 to 05:00

`tools/queue_md.py --preflight`: reports, for each queued row, whether its
system is built, whether its checkpoint is from equilibration or production,
where the production clock has reached, how far the trajectory ran on past the
checkpoint, which barostat the resume adds and which platform the checkpoint
needs. Read-only, no OpenMM context, safe beside a production row, non-zero exit
on a row that cannot run.

Also fixed: the row log was opened `"w"`, so a resume erased the row's own
history. b2ar_apo's five hours of equilibration went that way at 04:20 and
survive only in production.log. Now appended, with a separator. The live queue
holds the old code, so the four existing row logs are copied to
`data/md/rowlogs/` before it truncates them.

Found by running it: both unstarted rows are already built, ligands included
(45 CAU atom lines, 69 STI), so neither can fail at build time, which was the
risk the `ligand` field's comment records. And the first version of the
frame-drop figure was arithmetic from the current checkpoint interval, which put
abl1_apo at 4.2 ns when its checkpoint is at 4.0; replaced with the measurable
thing, how long the trajectory went on after the checkpoint was written.

FINDINGS 4ah. Five tests added, 944 passed, 1 skipped.

### 11 October, 05:00 to 05:40

`tools/stream_analysis.sh`: the liveness check looked once, and the queue runs
each row as a fresh subprocess, so a look taken between one runner exiting and
the next starting would have ended the streaming after the first row and run the
final pass, with the bilayer QC and G6 in it, against a partial trajectory. The
04:19 restart opened the same window for a minute and the loop survived only
because it was asleep. Now checks for the queue as well as the runner, and
requires an absence to persist over twelve looks five seconds apart.

Made testable to prove it: the script returns early when sourced, as the Python
tools guard `__main__`, and the grace is two environment variables. Four tests
drive `producing` against a stubbed `pgrep`. With the grace set to one look,
which is the old behaviour, the gap case returns false; with two, true.

Loop restarted on the fix at 04:38 and its first pass picked up the
post-truncation frames coherently: b2ar_apo at 177 frames against production at
17.7 ns.

FINDINGS 4ai. Four tests added, 948 passed, 1 skipped.

### 11 October, 05:40 to 06:20

Giving the new `short:` status a colour found that its class, `md-tag`, was
defined in no stylesheet. Scanning every template for `md-*` tokens and
differencing against brand.css found five more: `md-stats`, `md-scrub`,
`md-readmore`, `md-sub`, `md-viewer`. `md-stats` is the figures block at the top
of every index card and every run page, rendering as a browser-default
description list instead of a row of labelled figures. All six defined, and the
difference is now a test with five named exceptions, one of them a genuine
JavaScript hook. Confirmed by screenshot at desktop and phone width.

Also checked and ruled out: disk headroom. 1.3 TiB free against 2.4 GB of
trajectory still to be written across the four rows.

FINDINGS 4aj. Three tests added, 951 passed, 1 skipped.

### 11 October, 05:17 to 06:05

Audited section 11's deliverable list. The CLI has all four subcommands plus
`controls`; section 7's limitations name all four the spec requires, plus three
more; the favicon set is present; `deploy/` is written. `deploy/runpod/` was not
empty as I first read it, a truncated `ls` on my part: `job.sh` and `sync.sh`
were written on 10 October. The pixi lock already covers linux-64 with a CUDA
build of OpenMM 8.5.2.

Added the README and `.env.example` that were missing, and reading the scripts
to write them found a real fault. `job.sh` ran the queue, then `status=$?`, then
pushed to durable storage; under its own `set -euo pipefail` a failing command
exits immediately, so the assignment was unreachable and the push never ran on
the one path that needed it. A crashed row or a reclaimed pod lost the
checkpoint it had just earned. Now `|| status=$?`, with the push outside any
success branch.

Also added `--require-platform` to `tools/run_md.py` and `tools/queue_md.py`,
and `job.sh` passes `CUDA`. The fallback ladder is right on a machine you own
and a trap on a pod billed by the hour, where it drops a 1 µs row to a fiftieth
of the speed and goes on billing.

FINDINGS 4ak. Five tests added, 956 passed, 1 skipped.

### 11 October, 05:50 to 06:25

Audited section 10's QC list. Five of six items have tests: numbering maps,
metric functions on synthetic geometries, state-score symmetry, the membrane
builder on 5 ps, and the resume. The sixth, "any metric with all-NaN for a
motif is a test failure", had nothing, and the check that looks as though it
would cover it cannot: an all-NaN metric passes the round-trip verifier
perfectly, because the analysis file says NaN too.

`verify_live_metrics` added, reporting three failures: a metric NaN in every
frame, a value with no state score derived from it, and a motif that
contributed no rows at all. The last needed its own query, because a motif with
no rows leaves nothing for a `GROUP BY` to find. Some NaN frames are
deliberately not reported; all of them is the failure.

Live database: 2,536 rows, no NULL in any, every motif metric with a reading.
The final pass echoes both verdicts rather than only the round trip's.

FINDINGS 4al. Four tests added, 960 passed, 1 skipped.

### 11 October, 06:22 to 06:55

Section 10's "threshold lines labelled on every plot": the state score plot
shaded an active and an inactive band and never said where the edge was. Now
drawn as a dotted line labelled `s = +0.20` and `s = -0.20`.

Fixing the label found the larger fault. `state.DECISIVE_S` was named once and
written out as a literal three more times: the summary label in `store.py`, the
pair verdict in `web.py`, and the band shading in the run template. All three
read the constant now, the template through the view. The test probes the
behaviour either side of the threshold rather than asserting the source
contains a name.

FINDINGS 4am. Three tests added, 963 passed, 1 skipped. Screenshot confirms the
labelled lines.

### 11 October, 06:54 to 07:30

First look at the live beta2AR page, served from the production database rather
than the fixture site. It renders: 200, no console errors, every panel
populated, membrane health 68.0 then 67.3 A^2 with working verdict badges, tilt
13.6 then 16.1 degrees.

The structure viewer was empty. Partly a timing artefact, the membrane system
taking longer than two and a half seconds to build, but with a longer wait the
receptor appears a fifth of its proper size at the bottom of the panel. The
trajectory's bounding box is 186.7 A in y against a 131.2 A box: the receptor
sits on the periodic boundary and ten of twenty-four frames are written into
the other image.

Established first that it affects nothing measured: the protein is whole in
every frame, largest y span 58.3 A, and every metric here is within-frame. The
numbers are unaffected.

`periodic_shifts` added to tools/analyse_run.py, applied to the viewer bundle
and passed to the lipid bundle so the bilayer moves with the receptor rather
than being left a box length away. Collapses y from 186.7 to 60.8 A on the live
trajectory with every per-frame span bit-identical. Returns zeros for a
non-orthorhombic box rather than guessing.

FINDINGS 4an. Eight tests added, 970 passed, 1 skipped.

### 11 October, 07:29 to 08:00

The streaming pass at 07:23 rebuilt the bundle on the new code and logged
"11 of 25 frames were in a different periodic image". Measuring the result
rather than assuming it showed the receptor fixed, y span 60.8 A, and the
bilayer made worse: y span 269 A, centre alternating 65/196.

Applying the protein's shift to the lipids was the wrong operation. The bilayer
was already in one piece, so half its frames were moved a box length. Replaced
with the minimum-image convention about the receptor's centre, which also fixes
lipids wrapped in z that the first attempt never touched: y centre now tracks
the receptor and z span falls from 120 A to 62.

FINDINGS 4an corrected rather than appended to. 970 passed, 1 skipped.

Confirmed at 07:53 on the rebuilt bundles: receptor y span 60.8 A (was 186.7),
bilayer z span 62.3 A (was 120.1), bilayer centre following the receptor's
frame by frame, every phosphorus within half a box of it. On the page the
receptor draws as a seven-helix bundle down the membrane normal, and with the
section 6.2 toggle on it sits inside a contiguous field of phosphorus. No
console errors either way.

Also checked: GATE_OPEN.md names all four gates that apply, which is section
8's phase 7 requirement. G4 open and pushed, G3 waiting, G1 will not open,
G2 not needed.

### 11 October, 08:05 to 08:45

README: the positive-controls table listed five of the six controls, missing P6
(passes) and P5 (pending). Both added, with the pending one stated as pending,
and a test that the table names all six.

Live reference browser: `resolution` was NULL in all 136 rows under a column
the page renders and section 7's schema declares. Both catalogues return it and
both filters already parse it, for the 2.8 A and 3.6 A cuts; it was simply not
carried into the manifest. Now carried, omitted rather than nulled where
unreadable. Rebuilt: 136 of 136 rows have one, P1/P2/P3 unchanged at 17/17,
17/17, 100/100 at 11.39 A.

One suite run failed `test_an_upload_is_kept_when_it_is` and passed on a re-run.
Not flaky. `keep_run` built ids with `token_urlsafe`, whose alphabet contains
`-` and `_`, so one id in thirty-two began with a hyphen while `_RUN_ID`
requires an alphanumeric first character. Such a run was written to disk,
reported as kept, and then listed nowhere and reachable at no URL. Measured
3.03% over 20,000 draws. The first character is drawn from the alphanumerics
now; the regression test makes 3,000 draws, because one draw passes 31 times in
32.

FINDINGS 4ao. Six tests added, 977 passed, 1 skipped.

### 11 October, 08:22 to 08:50

Confirmed the resolutions landed: 136 of 136 in the database, ABL1 1.53 to
2.80 A and beta2AR 2.40 to 3.60 A, which are exactly the two filter cuts, and
the column renders on the live page.

DECISIONS.md had fallen behind: three of tonight's choices were in FINDINGS and
not there. D026 (a row is judged on what it produced), D027 (the platform
fallback is turned off on rented hardware), D028 (the periodic correction is
for the picture, and the bilayer is imaged rather than shifted).

Three documentation guards added, two of them for faults already hit. FINDINGS
section ordering, which I broke three times tonight by anchoring an insert on
the wrong heading. Decision numbers in sequence with no gaps or repeats. And
every D0NN cited in FINDINGS existing, which immediately found that D001 to
D011 were written `## D001 - Title` and everything after `## D012: title`, so a
check written against one form reported three real decisions as missing. All
eleven normalised to the colon form.

980 passed, 1 skipped.

### 11 October, 08:50 to 09:15

`tools/make_figures.py` runs only in the final pass, so nothing had executed it
since the CSS edits, and it reads brand.css for its palette. Ran it: seven
figures, exit 0, palette parsing intact.

The first real per-motif production figure showed the fault. It shaded its
active band from zero, captioned "above zero: active-like", while
`state.DECISIVE_S` is 0.2 and the label rule and the explorer both cut there.
A motif at +0.1 was drawn active in the document and called intermediate on the
page. Now uses the constant, with both edges drawn as labelled dotted lines and
the caption naming the intermediate region. The white band is where TM5 spends
most of the second half of this run.

Adding the import broke the tools-import check, which read the trailing
`# noqa: E402` as part of the imported name. A fault in the check; fixed to
drop anything after a `#`.

FINDINGS 4ap. 980 passed, 1 skipped.

### 11 October, 09:27 to 09:55

The bilayer figure is right: labelled window bands on both panels, units on
every axis. The page beside it read "68.0 A^2  outside expected" against a note
saying the window is 60 to 68. The measurement is 68.05 and the window is
inclusive to 68.0, so the badge is correct and one decimal place rounded away
the digit that decides. The two cells that carry a verdict are at two decimals
now.

The note also wrote the three windows out as prose while the badges read the
constants. Same shape as 4am and 4ap. The view passes them and the note renders
them, so the sentence a reader checks is the one the badge was decided by.

FINDINGS 4aq. Two tests added, 982 passed, 1 skipped.

### 11 October, 09:58 to 10:20

Read the per-motif numbers again at 30.8 ns, 309 frames, against 5.0b's reading
at 14.6. Whole protein +0.319 to +0.266, still active. TM6 falls from +0.534 to
+0.377 and H8 from +0.474 to about +0.21 in the first few nanoseconds and then
holds; both are about one standard deviation of their own scatter. NPxxY rises
slightly, TM7 flat. Nothing has crossed a state boundary.

And a correction to myself. I described TM5 as drifting down and crossing zero
from about 12 ns on, twice, having read it off the figure. Its halves are -0.034
and +0.016, so it is flat: it crosses zero because it sits on zero. The halves
were already in the summary and I looked at the picture instead.

FINDINGS 5.0c. 982 passed, 1 skipped.

### 11 October, 10:28 to 10:50

Drew each motif's half-means on the per-motif figure, as a thick flat segment
over each half in the motif's colour. The figure without them misled its own
author last tick: TM5's trace wanders through +/-0.3 about zero and reads as a
drift, and its two segments sit level. Two segments rather than a fitted slope,
because a slope invites being read as a rate.

Verified against the stored summary rather than by eye: the figure's halves and
`summarise`'s agree to three decimals on the live run, and a test asserts the
two paths agree on a fixture. The figure's query already used SELECT DISTINCT,
so it counts frames and not rows; that now has a test, because without it the
means would start to drift the day the two metrics of a motif carried different
scores.

985 passed, 1 skipped.

### 11 October, 10:59 to 11:20

The half-mean segments are on the explorer chart now as well as the figure. The
page is the deliverable; the figure is only where I noticed the problem. Same
split, `len // 2` of the plotted series, hidden from the legend, explained in
the note under the chart. Confirmed on the live page: TM5's pair level, TM6's
stepping down, no console errors.

The note under that chart wrote "+0.2" as prose beside bands drawn from
`DECISIVE_S`. It reads the constant now. Same shape as 4am and 4aq, this time
in the sentence explaining the bands.

987 passed, 1 skipped.

### 11 October, 11:30 to 11:55

Made the threshold-in-prose fault a check instead of finding it a fourth time
by hand. It reads every reference-state template, strips scripts, Jinja
comments, expressions and statements, and tags with their attributes, and
searches what a reader sees against ten named constants.

Found four more in three templates: the index giving the pass band as +/-0.2,
About giving the R-spine contact as 7.7 A, the reference browser giving the
KLIFS cut as 2.8 A and GPCRdb's as 3.6 A. All four correct today, all four
would have stayed unchanged if the constant moved. All four render the constant
now.

Only forms with a decimal point are searched: IONIC_LOCK_A is 4.0 and the run
page says "alphaC-Glu+4 Calpha", which is a residue offset, and a check that
reported that would be turned off within a week. A second test proves the scan
can fail, against a sample holding the same number in an attribute, a comment,
an expression, a script and a sentence. It finds the sentence.

FINDINGS 4ar. Two tests added, 989 passed, 1 skipped.

### 11 October, 12:01 to 12:25

Read the About page end to end for the first time and found a wrong sentence:
"the receptor is better served, with 72 active and 28 inactive structures"
where the ensembles hold 84 and 32. The ABL1 half, 5 and 14, is right.
Reproducing both needed the rule `sources` builds ensembles with rather than
the catalogue's column: ABL1 active is DFG-in and alphaC-in. The page counts
them from the database now and drops the comparison when there is none.

Named the published DFG cut `kinase.DFG_PUBLISHED_A`, since the state map's
dashed line and the comparison rule were two copies of 11.0. Adding it to the
scanned list immediately found two more sentences holding it, on the About and
run pages. Both render it now.

That edit was half-applied: the patch matched on a line that had just stopped
being unique, so it stopped after editing both templates and before the view,
and the run page returned 500 with `'dfg_published_a' is undefined`. Caught by
rendering all five pages; `tests/test_pages.py` catches it too.

FINDINGS 4as. Three tests added, 992 passed, 1 skipped, 12 page tests green.

### 11 October, 12:38 to 12:55

Checked what G3 would send. `deploy/deploy.sh` is right about the hard part: it
excludes data/md, data/cache and data/refs and keeps data/refstates.db,
data/analysis and data/systems, by listing paths rather than filtering on
.gitignore, since the database is itself gitignored and the app renders an
absent one as an empty section.

Its size claim was stale: "under 5 MB" against 7.3 MB measured, with one row a
third done. It prints the measured figure before syncing now. Two tests, the
second scoped to non-comment lines because the script's explanation of why not
to use the gitignore filter contains the filter.

Also checked and clean this tick: the index cards count references from the
database (20 and 116, the "4" I had seen was the fixture site), and
config/benchmarks.toml is written by hwprobe in the final pass so its staleness
between passes is expected.

FINDINGS 4at. Two tests added, 994 passed, 1 skipped.

### 11 October, 13:09 to 13:35

The committed figures in `data/figures/` were drawn at 04:04 by a version of
`tools/make_figures.py` that no longer exists: bands from zero rather than
`DECISIVE_S`, no half-means. The final pass regenerates them and does not run
for two more days, so the repository was showing superseded output that
FINDINGS describes. Regenerated: four of seven changed, three byte-identical
because abl1_apo has produced no frames while it waits in the queue.

A test now compares the last commit date of the script against the last commit
date of the directory. Written first against file mtimes, which a checkout
destroys, so it meant whatever order git unpacked things in. It failed on being
written, correctly, and passes once both go in together.

GATE_OPEN's acceptance table updated: it said six figures where there are
seven, and it now records tonight's guarantees as their own rows, including
section 10's all-NaN rule, the threshold labelling and the prose check, the
short-row recording, --require-platform, the streaming grace and the periodic
correction.

FINDINGS 4au. One test added. 994 passed with the new test red until this
commit lands, which is the state it is asserting about.

### 11 October, 13:36 to 14:00

The user asked for all four of the remaining options, including the reference
overlay I had recommended against. First of them: section 6.2's cartoon
coloured by motif, now built.

No API route exists: the vendored Mol* bundle is minified and exports no
MolScript, so a component per motif is unreachable from the page. The colour
theme registry does expose `uncertainty`, which reads the B-factor, so the
number travels in the structure. `tools/analyse_run.py` writes each atom's
`--motif-N` number there and the page colours by it from brand.css. Display
only; `loaded` is untouched. Overlaps go to the smaller motif, so xDFG keeps
380-383 inside the activation loop.

Wrong twice, both caught by looking: the colours were CSS strings where Mol*
wants integers, which drew a near-black cartoon with no error; and the scale
runs backwards, which gave every motif its opposite number's colour and looked
deliberate. The second was confirmed by predicting the whole mapping from atom
counts and matching it against the render rather than by flipping it.

D023 marked as superseded in part, D029 records the method. 1001 passed,
1 skipped.

FINDINGS 4av. Six tests added.

### 11 October, 14:00 to 14:30

Second of the four: section 6.2's reference overlay, which D023 had declined.
Built with the check that answers D023's argument rather than merely being
careful about it. RMSD is invariant under a rigid transform applied to both
sets, so superposing the reference onto the run must reproduce the distance the
analysis got superposing the run onto the reference, and the per-frame record
names which reference was nearest for each motif. Disagreement beyond 0.005 A,
or nothing to check against, means no overlay and a logged reason.

The check changed the design. Choosing by lowest mean motif RMSD picked 5MO4
for abl1_apo's inactive state, which is nearest for no single motif, so the
guard refused it. Selection is now by the residue weight a reference is nearest
for. abl1_apo takes 2HZ4 active and 2HZ0 inactive, both verified.

`representationParams` is ignored by this bundle, so the ghosts arrived in
default chain colours at full opacity; the component refs are differenced
across the load and rethemed after it, at alpha 0.3 in the state colours.

D030 records it, D023 marked superseded. 1007 passed, 1 skipped.

FINDINGS 4aw. Six tests added.

### 11 October, 14:30 to 15:05

Third of the four: D025's deferred core intersection, built. `intersect_core`
takes the KLIFS scaffold positions the run and every reference all resolve and
resolves them back per structure, so every core is the same length and the i-th
coordinate is the same place in the fold. CDK2 goes from 28 of 29 aligned to
29 of 29 on 41 shared positions; the intersection across all 113 CDK2
structures is 38, D025's estimate.

Checked rather than assumed: among the 109 CDK2 structures with 42 positions
there is exactly one distinct position set, and the same for ABL1's 19, so the
old same-length cores were comparing like with like. This is more coverage, not
a corrected measurement. ABL1 re-analyses bit-identical.

Uniform cores then left modal length deciding nothing, so the target fell
through to sorted order and picked 1F5Q on CDK2, residual 0.94 A. Target is the
medoid now: CDK2 onto 3BHT at 0.40, ABL1 onto 2HYY at 0.425 against 2F4J's
0.568.

I was wrong that the target could not move the scores: each reference is fitted
onto it separately, so it is not one rigid motion. Measured on abl1_apo over 43
frames, largest distance change 0.0146 A, largest score change 0.0016, no label
moves, P4 four of four. abl1_apo's P4 margin is thin at +0.304 against +0.300,
as it was before at +0.305.

An existing test pinning the sorted-order tiebreak failed and was right to;
rewritten to assert what D025 protected, that a uniform set loses nobody.

D031 records it, D025 marked superseded. 1011 passed, 1 skipped.

FINDINGS 4ax. Four tests added.

### 11 October, 14:40 to 15:15

Audit tick, and the overlay check paid for itself three more times.

Verify is green: every stored number agrees with its analysis file and every
motif metric has a reading, so the streaming pass has re-ingested both rows on
the new code.

beta2AR's overlays were both refused, TM5 measuring 1.6735 A against a stored
1.6797. Mine, not the data's: `traj.analyse` superposes the frame once onto one
target core and measures every reference in that common frame, where I was
fitting each reference onto the run directly. The ghost now follows the same
two-hop path, own frame to common to run, so the distance is identical by
construction. A separate test checks the composition against doing it the long
way.

That exposed the next: the inactive reference chosen was 9CHU, nearest for TM5
at 1.6 A and 27 A away averaged over its motifs. A motif read in the wrong
place never wins a minimum so it never reaches the analysis; an overlay draws
the whole structure. Selection is now best mean among references nearest for at
least one motif. beta2AR takes 8GE8 at 1.46 A and 9RKI at 2.62 A.

Open and recorded: abl1_apo's inactive ghost is still refused, 6NPV's hinge
measuring 0.5593 against 0.5053. The active ghost on the same run verifies, and
both beta2AR ghosts verify without being the alignment target, so the
composition is right in general and this is specific. Not found. The tolerance
stays at 0.005 A: a check tuned until it stops refusing is not a check.

1012 passed, 1 skipped.

### 11 October, 15:05 to 15:25

Chased the open overlay discrepancy and ruled out the leading candidate. The
own-frame to common-frame hop was refitting the reference onto the target core,
where `align_to_common_frame` fits every reference onto the medoid, a different
structure. It is taken from the reference's own stored core now, which recovers
the alignment's transform exactly whatever target it used. abl1_apo's hinge
moved 0.5593 to 0.5600, so that was not the cause; the change stays because it
removes a dependency that should not have existed.

A diagnostic written to isolate the rest was unfaithful, loading the raw
solvated trajectory rather than the prepared one, so it never reproduced the
analysis's core. Stopped there rather than building a bespoke harness: the
system refuses what it cannot confirm, which is the behaviour that matters.

Two test fixtures corrected along the way. The second gave every reference its
own frame as the common one, which is incoherent when the overlay takes the
first code as the target; it now fits every core onto the first, as the real
alignment does.

1012 passed, 1 skipped.

### 11 October, 15:30. Spec re-read against the code, for what the landing data needs

Queue state and projection from measured speeds (membrane 59.7 ns/day, water
202 ns/day), equilibration included for the two unstarted rows:

    b2ar_apo        44.7 of 50 ns   2.1 h   Sun 11 Oct 17:38
    abl1_apo         4.0 of 50 ns   5.5 h   Sun 11 Oct 23:06
    b2ar_carazolol   0.0 of 50 ns  24.1 h   Mon 12 Oct 23:13
    abl1_imatinib    0.0 of 50 ns   6.1 h   Tue 13 Oct 05:20

Sections 5.4 and 5.5 re-read metric by metric against the stored payloads. Both
are complete: the kinase has the DFG plane, alphaC salt bridge, alphaC
displacement, A-loop, P-loop (fallback definition, logged), hinge backbone
flags and R-spine; the receptor has TM3-TM6 and its alternate pair, the ionic
lock, TM6 rotation and centroid shift, NPxxY and toggle chi1/chi2, the PIF
connector, helix 8's angle to the membrane plane, the kinks at 6.50 and 7.50
and the sodium pocket distances.

**Three gaps found, one of which fires when the remaining rows land.**

1. **Ligand RMSD and contact persistence (5.7) are not implemented.** Both
   remaining rows are ligand-bound, b2ar_carazolol with CAU and abl1_imatinib
   with STI, and 5.7 asks for "the same analysis plus ligand RMSD and contact
   persistence". `flexappeal/measure.py` already has `residue_contacts`,
   `fingerprint` and the pose machinery, and `traj.kept_atoms` keeps the ligand
   deliberately, so this is wiring rather than a new measurement.

2. **The sodium pocket cannot report its ion.** `gpcr.sodium_pocket` takes
   `ions` and reports `sodium_present` only when given positions, which is
   right. But `traj.load` drops solvent, ions and bilayer on load, so no ion
   positions ever reach it and 5.5's "record whether a Na+ ion sits within
   3.5 A of D2.50" is unreachable as the pipeline stands.

3. **The motif timelines are not 6.2's small multiples.** The spec asks for one
   small panel per motif metric with shaded 5th-95th percentile bands from each
   reference ensemble; what exists is one combined state-score chart. Recorded
   as a deviation in 4u and still open.

Optional and explicitly deferrable: the Kincore DFG dihedral label (5.4 says
"or N/A") and the GPCRdb activation index (5.5 says implement only if the
published definition can be verified).

Checked and already working for the landing data: the pair comparison needs a
second run per system and both `start_state` values are already in
`tools/ingest.py`; clicking a state-map point jumps the viewer; the final pass,
P5 and P6 are built and rehearsed.
