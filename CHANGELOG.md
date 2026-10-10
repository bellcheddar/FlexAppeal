# Changelog

Dated history. Goals and non-goals live in `PROJECT_PLAN.md`; the roadmap is the
checkbox list in `README.md`.

## 2026-10-10

### Changed

- **Every sub-score row on the grade meter has a bar, and its number in one
  column.** The energy row had the kcal/mol inside the bar cell and an em dash in
  the score column, so it read as a gap beside five filled rows, and a
  convergence score of 0.0 drew an empty track that looked like a missing panel
  rather than a measured nought. The measurement now sits under the label, the
  bar cell holds only a bar, and the score column holds exactly one number per
  row.

- **The example intro now describes the run it is showing.** It still said 1AKI,
  hen egg-white lysozyme, 10 ns, 19,433 atoms, an M1 Max and a 200 ps
  equilibration. It is 181L, T4 lysozyme L99A with benzene, 1.80 A, 100 ns,
  25,463 atoms with 7,600 waters, an M2 Ultra and a 1 ns equilibration, and the
  ligand is parameterised with OpenFF 2.2.1 and AM1-BCC charges.
- **The trajectory clip's caption no longer claims to be this page's run.** The
  clip was rendered from the 10 ns apo run and is decorative, so it now says it
  is illustrative rather than naming the Example tab's own data. Re-rendering it
  from the committed run would let it say otherwise.

- **The Example tab is now a 100 ns run of T4 lysozyme L99A with benzene.**
  Previously a 10 ns apo run of hen egg-white lysozyme, which renders the
  whole-system panels and not one binding-site panel, so eleven panels and the
  2D depiction were shipping untested by the page whose job is to prove the
  product works. The new example is PDB 181L, 25,463 atoms with water, 100 ns of
  production in 8.05 h at 298.3 ns/day on OpenCL, and it exercises every
  binding-site panel the Analysis tab can draw.

  Two panels stay empty, and the page says why where each would have been.
  Benzene has no nitrogen, oxygen or sulphur, so there are no polar contacts to
  occupy; and it has no rotatable bond, so there is no torsion to flip. Both are
  facts about the molecule, not gaps in the measurement.

  The 10 ns apo example is removed rather than kept beside it. Two committed
  examples means two multi-megabyte results files in a public repository for
  one Example tab, and the apo case is already covered by a test fixture
  that asserts a run with no ligand shows no binding section.
- **The example's length guard reads the example's own `config.json`.** It
  asserted 10 ns literally, so swapping the committed example failed a test whose
  message named the wrong problem. It now compares the results file against the
  production duration the committed bundle was built to do, which is the thing
  worth guarding, and keeps the 95% floor against a truncated run.

### Added

- **The MM/GBSA term panel, which was the last one DeepFlex had and this did
  not.** A horizontal bar of the difference block -- van der Waals,
  electrostatics, the GB desolvation penalty, the surface term, and the two sums
  -- with the standard deviation over the sampled frames as error bars rather
  than the standard error, because the spread is the interesting part and a
  standard error would shrink with frame count and hide it.

  Audited against DeepFlex rather than guessed at. Of its fifteen panels, twelve
  were already here under different names, two need a package and cannot exist
  for one complex (the grade distribution and the league table), and this was the
  only real gap. FlexAppeal also draws three DeepFlex does not: interaction
  occupancy per residue, the pose clusters and the per-atom ligand fluctuation.

  The bonded terms are deliberately not drawn. They cancel to zero in a single
  trajectory, so five bars at zero would say nothing; that they cancel is a check
  on the calculation, which `energy._check_cancellation` makes, not a finding.

- **`flexappeal/energy.py`: MM/GBSA binding free energies.** The last of the six
  sub-score inputs FlexAppeal could not produce. A ligand bundle now computes it
  at the end of its analysis, and `flexappeal package --mmgbsa` does it from the
  command line.

  It shells out to AmberTools' `MMPBSA.py` rather than computing GB in OpenMM,
  which is more moving parts and the right choice: the existing cohort's energies
  came from `MMPBSA_igb_2` on a water-stripped GAFF2 topology, and a percentile
  across a cohort measured two ways is not a percentile. Single trajectory,
  `igb=2`, mbondi2 radii, salt at the run's own ionic strength.

  Two failures found while building it, both silent:

  A topology written from an OpenMM system carries `RADII` and `SCREEN` all zero,
  because the system it came from had no implicit-solvent force to hold them.
  MMPBSA.py then divides by them and reports every `EGB` as `nan` without
  failing: the first real run produced a complete, well-formatted results file
  whose every total was `nan`. `changeRadii mbondi2` fixes it and the parser now
  refuses a `nan` and names the cause.

  The simulated system constrains every bond to hydrogen, so those bonds carry no
  force-field term and ParmEd refuses to write a topology from it. The energy
  calculation builds its own unconstrained system, which is the right physics for
  a single-point energy anyway.

  The result is checked rather than trusted. In a single trajectory the bonded
  terms must cancel exactly, because the receptor and ligand are read out of the
  complex's own coordinates; a residual means the three topologies disagree about
  the molecule, which produces a plausible number rather than an obvious failure.
  FlexAppeal refuses the energy in that case.

  Benzene in T4 lysozyme L99A: **-9.8 kcal/mol**, SEM 0.66 over 200 frames, in
  about 80 seconds. Van der Waals -11.3 dominates with electrostatics at -1.1,
  which is what a non-polar ligand in a buried non-polar cavity should look like;
  experiment is about -5.2, and single-trajectory MM/GBSA with no entropy term is
  expected to be over-negative by roughly that much. DeepFlex's ingester now
  draws 50% of the energy page for the one-complex file against 0% before, with
  twelve energy terms.

  And scored, which took two goes. `S_E` is normally a percentile and a lone
  complex has no cohort, so it now falls back to absolute thresholds in DeepFlex's
  `grade.toml` (D078): full marks at -60 kcal/mol, zero at 0. Zero is physics; -60
  is the strong end of the 27 complexes that pipeline has measured. The fallback
  runs only where there is no cohort, is marked on the row, and the limitation
  travels with it -- within one pocket these energies do not separate good
  compounds from poor ones, so nothing anchors -60 beyond the shape of its own
  distribution. Proved inert on DeepFlex's real package: no score, grade,
  sub-score, gate or rank moved on any of the 27.

- **A ligand run now writes the version 2 package itself.** Spec section 0 was
  updated today: a run with a ligand emits version 2, including a single system,
  because a one-complex package is a first-class graded run and version 1 has
  nowhere to put a ligand measurement. Previously only `flexappeal package` could
  produce one, which meant a bundle user needed a FlexAppeal checkout to get a
  grade at all.

  A ligand bundle now carries `fxalib/`: `fxa.py`, `fxa2.py`, `measure.py`,
  `options.py` and `package.py`, copied byte for byte out of the installed
  package rather than templated, with a test asserting the bytes match. They are
  stdlib-only at import time and a ligand bundle's environment already declares
  mdtraj, scipy and rdkit, so nothing new is pulled. An apo bundle carries none
  of it and stays 80 KB; a ligand bundle goes to 188 KB.

  The version 1 file is still written and is not a duplicate. It carries the
  whole-system analysis -- RMSD, radius of gyration, SASA, secondary structure,
  the contact map, essential dynamics -- plus the topology and a decimated
  trajectory for the viewer, none of which a version 2 package has anywhere to
  put. Measured directly: a one-complex v2 of the benzene run is 2.5 MB of
  binding-site material across 23 members and has no whole-system panel in it.
  The two are different reports of one run.

  The package is measured from the run directory, not from the version 1 file:
  that payload is decimated for viewing, and a grade must not rest on thinned
  data. Verified against DeepFlex's own ingester, which accepts the one-complex
  file as `graded (fxa/2)` and scores it.

- **A grade meter on the Analysis tab, and `flexappeal/grade.py` behind it.** The
  page used to carry a paragraph saying there was no grade on it, that this was
  deliberate, and that a study package was the way to get one. That was
  unnecessary. DeepFlex already grades a complex whose sub-scores are incomplete
  by renormalising the missing weight across the rest rather than filling it with
  a neutral value, so one complex with no cohort is a case its own arithmetic
  covers. Withholding the number and pointing at another site was worse for a
  reader than running the same arithmetic here.

  The module is a transcription of DeepFlex's `config/grade.toml` and
  `pipeline/grade.py`: the six weights, the five ramps, the bands with their
  sub-bands, the sparse-fingerprint fallback that may only lower a score, and
  both hard gates. `tests/test_grade.py` pins every one of those numbers against
  the values published in the spec, because a transcription is worth having only
  while it agrees with the original.

  What stays DeepFlex's is what needs the cohort. The energy sub-score is a
  percentile among the other ligands against the same receptor, so it is pending
  and its 15% is spread over the five that were measured; the meter shows the
  struck-through 15% beside the weight each sub-score actually bore. Both rank
  columns need the cohort too.

  Nothing from this module is written into a `.fxa`. Spec section 7 reserves the
  grade to DeepFlex on ingest so a file carries one copy of the judgement; this
  one is rendered and discarded.

  The committed example grades **A, 83.6**: pose 73.5, occupancy 100,
  interactions 100, convergence 40.4, conformation 100.

- **RDKit on the deployed site.** It was left out of `requirements.txt` on the
  reasoning that a chemistry toolkit was too much to carry for two charts on a
  small box. Both halves of that were wrong. It is three panels and the 2D
  depiction, which is most of what a reader comes to a ligand run for, and the
  box has 89 GB free against a wheel of about 35 MB. The note in that file
  describing a 3.8 GB droplet shared with four other apps was stale by some
  margin.

  The guarded import in `measure.py` stays, because a generated bundle solves its
  own environment and may not have it.

- **`rebuild.sh` captures the analysis half too.** It called
  `scripts/terminal_capture.py` without `--analysis-log`, so a rebuilt example
  rendered the run's terminal output and nothing of what turned it into a results
  file. Nine captures now rather than seven.
- **A test that the committed example is a bound run.** It asserts all eleven
  binding panel ids, the depiction and the mobility legend are on `/example`, so
  a swap back to an apo example cannot quietly stop testing them.

### Removed

- **The two measurement notes shown under that paragraph.** The bond-order
  inference caveat and the no-rotatable-bond note. They are still recorded on the
  measurement and still travel into a package, where DeepFlex handles them; they
  were caveats about how a number was made, printed above the numbers, which is
  not where a reader needs them.

- **The 10 ns apo example.** Two committed examples means two multi-megabyte
  results files in a public repository serving one Example tab, and the apo case
  is already covered by a fixture-backed test asserting that a run with no ligand
  shows no binding section. `examples/README.md` and the capture script's usage
  now describe the run that exists.
- **Both binding heat maps rendered in narrow cards** while every version 1 heat
  map is wide. A square matrix in a half-width card is drawn 340 pixels tall and
  then letterboxed, which wastes the width it was given and makes a block pattern
  harder to read than the same data in the wide card the version 1 maps get.


### Fixed

- **The pose RMSD charged a symmetric ligand for its own symmetry.** Benzene's
  six carbons are chemically identical, so a 60 degree ring rotation moved every
  atom by the ring radius while changing the pose not at all. The RMSD now
  minimises over the automorphism group of the ligand's hydrogen-stripped graph,
  which is the standard correction and can only lower the number because the
  identity is always in the group.

  On the committed example: median pose RMSD **2.16 to 0.83 A**, `conv_js_mean`
  **0.41 to 0.048**, pose clusters **6 to 1** at full dominance. `S_pose` 73.5 to
  100 and `S_conv` 0 to 100. The convergence score of nought was almost entirely
  benzene reorienting into its own symmetry images and then stopping.

  Implemented in DeepFlex in the same change (D079), with the same two functions
  and the same eight reference counts pinned by a test in each suite. 24 of that
  package's 27 ligands carry symmetry, so it was re-measured there before
  landing: no grade, rank or gate moved and the worst score shift was 0.78
  points.
- **A rigid ligand reported zero strain and collected full marks for it.** The
  method restrains each rotatable torsion, minimises, minimises again free, and
  takes the difference. With no rotatable bond there is nothing to restrain, so
  the restrained force field *is* the free one: both minimisations start from the
  same coordinates and reach the same minimum. Benzene measured -4.7e-07
  kcal/mol.

  Reported as 0.0 that is zero excess strain, which is full marks on the
  conformation ramp -- 5% of a grade for a quantity nobody computed, and exactly
  the stand-in the format forbids. Both conformation inputs are now omitted, the
  sub-score goes pending and its weight renormalises. The strain and torsion
  panels are absent rather than flat, which is what led to the question.

- **The grade meter said the binding energy was not measured, when it was.** Two
  faults behind one wrong sentence. The energy was computed into the version 2
  package only, and the Analysis tab reads the version 1 payload and nothing else,
  so the number never reached the page it was wanted on. And the row's one piece
  of copy covered both "no measurement" and "measurement that cannot be ranked",
  which are different things.

  The run now computes MM/GBSA before it packs either file, so the version 1
  metrics carry it too and both files report one number rather than two samples of
  the same trajectory. The row reads `-9.8 kcal/mol ± 6.6, measured, but ranked
  only within a cohort`, and keeps `not measured here` for a sub-score that really
  has no measurement.

- **The convergence divergence was not on DeepFlex's ruler.** `conv_js_mean` was
  computed over an adaptive bin count, holding samples-per-bin constant to remove
  the finite-sample bias that inflates a short window. The bias is real -- a
  stationary 400-frame window reads 0.086 against a threshold whose full marks
  are 0.05 -- but the correction cost more than it bought.

  DeepFlex bins with `np.histogram_bin_edges(sp, bins=20)`, fixed, and the ramp
  this value is scored on is calibrated against exactly that. A divergence over
  any other number of bins is not comparable to it. The benzene example made that
  concrete: the Analysis tab scored this sub-score 40 while DeepFlex scored the
  same trajectory 0, purely because one read 500 frames and the other 2,000.

  Twenty fixed bins now, matching DeepFlex, with the bias reported rather than
  corrected: where a block holds fewer than ten frames per bin the note says so
  and says which way it pushes. The example's grade moves from A, 83.6 to
  **A-, 78.9**, which is the number DeepFlex gives for the same run -- the two now
  agree to a tenth of a point, the residual being 500 frames against 2,000.

- **The deploy could not repair the permissions it had just broken.** `rsync -a`
  preserves the source's modes, and a repo checked out on a Mac can be 0700, so
  the sync landed `/opt/flexappeal` as `drwx------` owned by the syncing user.
  The service user could then not traverse its own directory, and the deploy died
  on `env: './.venv/bin/pip': Permission denied` -- which names neither the
  directory nor the mode.

  The chown that would have fixed it ran *after* the pip install, so the deploy
  exited before reaching it, and every subsequent deploy failed identically until
  someone chowned by hand. The chown now runs first, the top directory is made
  traversable explicitly, and the sync lands files 644 and directories 755 rather
  than whatever the Mac had.

- **The Analysis page carried every frame of the run, and a 100 ns run made that
  obvious.** The first render of the new example was 24 MB: 18.7 MB of it the
  DSSP heat map, which carries a cell and a state name per residue per frame, and
  another 2.2 MB spread over the six convergence panels. It rendered correctly
  and the whole suite passed, which is why it needed measuring rather than
  reading.

  The version 1 panels now decimate for the page the way the binding panels
  already did, with the helpers moved up beside the other shared ones so both
  halves of the file reach them. Line panels keep 600 points, the DSSP heat map
  400 columns, and the residue axis is never touched because that is the axis a
  reader reads. The page is 1.7 MB, the 10 ns example came down too, and nothing
  on screen changed: 600 points is more than a chart a few hundred pixels wide
  can resolve.

  The file is untouched. Every frame is still in the `.fxa` and `/reanalyse`
  still works from the full arrays, which is the distinction: this is about the
  wire, not about the data.

### Measured

- **The OpenCL memory growth plateaus and then reverses.** The README has carried
  an extrapolation since the 10 ns run: 8 GB up over 10 ns, so about 80 GB for
  100 ns of the same system. The 100 ns run says otherwise. Resident memory grew
  linearly at 0.78 GB per percent to 46 GB at 64% of the run, then fell to
  34.4 GB by the end, adding 8.4 MB of swap on a 96 GB machine.

  This matters because the run's own warning fits a line to the early slope and
  acted on it: it reported memory growing by 1.3 GB per percent, projected
  105 GB still to come, and advised stopping and restarting from the checkpoint.
  Had that been followed it would have cost eight hours of a run that finished
  comfortably. The warning needs refitting against this trace before its advice
  is worth taking.

## 2026-10-09

### Added

- **`flexappeal/fxa2.py`: the `.fxa` version 2 format.** Version 1 is one system
  and one replicate and DeepFlex renders it as a flexibility run with no grade;
  version 2 carries a package of receptor-ligand complexes with the six sub-score
  inputs and is rendered as a fully graded report. It is a separate module rather
  than a branch in `fxa.py` because it is a different kind of file, not a bigger
  version of the same one, and version 1 has to keep working untouched.
  Specified in `DEEPFLEX_FXA2_SPEC.md`, which is the contract.

  This commit is the format alone: the salted identifiers, the member layout, the
  archive limits, the chemistry bans, and a reader that reads back everything the
  writer wrote. Nothing measures anything yet.
- **A writer that refuses what DeepFlex would discard.** Section 8 of the spec
  bans structure strings and chemistry members outright, and DeepFlex drops them
  on ingest and tells the uploader it did. Refusing them here instead means the
  mistake is caught by this repository's own test suite rather than by somebody's
  upload: `.sdf`, `.mol2`, `.prmtop`, `.xtc` and the rest are refused by suffix,
  and any manifest field whose name promises a structure string is refused at any
  depth, on the substring list in `fxa2.CHEMISTRY_ROOTS`. Compound names travel in
  clear, which is the line: a name says *which* compound, not what it is.
- **The conformance checklist from section 10 of the spec, as a test.** All
  fourteen checks, including the one that catches a het code left behind on a
  `TER` record. A `TER` carries the residue name in the same columns as an `ATOM`
  record, so renaming only the atom lines leaves the real code in the file; that
  exact bug shipped once on the DeepFlex side, so the check looks for a `TER`
  naming a residue no `ATOM` record names.
- **`FXA2_PROGRESS.md`**, the running log: current step, what the validator last
  reported, and every decision taken with its reason.
- **`flexappeal package`: a study, assembled from finished runs.**

      flexappeal package -o study.fxa run-dir-1 run-dir-2 ...

  A version 2 file needs many complexes and FlexAppeal runs one system per job, so
  packaging is a separate, deliberate act rather than something a job does. The job
  model is untouched, no version 1 file is altered, and a study can be assembled
  from runs done weeks apart. It takes either a finished run directory or a version
  1 `.fxa`, detected per argument, and expands a directory of either.

  The run directory is the better source: it carries the complete trajectory rather
  than a decimated payload, and it is the only one from which the ligand strain can
  be measured. The command says which it was handed and what that costs.
- **`flexappeal/measure.py`: the pocket, and the crop.** The binding site is
  reported twice, because the two are used for different things and conflating them
  is how a headline metric goes quietly wrong. Raw is the 10 A crop shell, which is
  what the viewer loads; filtered is the 6 A selection with the chain termini
  dropped, which is what the ligand RMSD is superposed on. A superposition frame
  containing a terminus wanders, and then every ligand RMSD measured against it is
  wrong in a way that looks exactly like the ligand moving.

  Checked against the chemistry rather than only the shape: the test asserts that
  ASP189, SER190, TRP215 and GLY216 are in trypsin's filtered pocket, which is the
  S1 site benzamidine actually binds. If the selection ever drifts, those are the
  residues that stop appearing.
- **The ligand is renamed to `LIG` on every record that names it**, including the
  `TER`. A `TER` carries the residue name in columns 18 to 20, the same columns as
  an `ATOM` record, so a rename that walks only the atom lines leaves the real het
  code sitting in the file under a hashed filename. Matching on the residue-name
  columns regardless of record type is what makes that impossible rather than
  merely unlikely.
- **`kind` is decided from what was measured.** Section 9 of the spec refuses a
  `graded` file with required inputs missing, so the writer checks every analysable
  complex against the required list after measuring and downgrades the file to
  `flexibility` if any is short, naming what was missing. A run that loses a
  measurement downgrades itself rather than claiming a grade it cannot support.
- **`--energies`, for the one thing FlexAppeal cannot compute.** `mmpbsa_gb_dg` and
  the interaction-energy terms need an MM/GBSA or LIE rescoring of every frame,
  which means rebuilding the parameterised system in implicit solvent, and nothing
  in the pipeline produces it. The keys are omitted rather than written as zero,
  which leaves the energy sub-score pending and DeepFlex says so on the page. Where
  a study has the numbers from elsewhere this carries them in, keyed on the compound
  name rather than the complex id, because the id is a salted hash that is fresh per
  file and nobody can write it down in advance. An absolute MM/PBSA block is
  refused: the three absolute blocks are each around -100,000 kcal/mol and only
  their difference is a binding quantity.
- **The scored window, and the two sub-scores that are half the grade.**
  `equil_start_frame` is chosen by Chodera's rule, taking the truncation that
  leaves the most effectively uncorrelated samples. That is the quantity pymbar's
  `detect_equilibration` maximises, so this is the same method without the
  dependency rather than a cheaper substitute for it. Evaluated at about fifty
  candidate truncations because the curve is smooth and the search is otherwise
  cubic in the frame count, which on a 2,000-frame run costs minutes to produce one
  frame index.

  Searched unconstrained and clamped to half the run afterwards, not searched
  inside the cap: clamping a bounded search cannot tell a run that equilibrated
  just under the cap from one that never equilibrated at all, and those two deserve
  different words on the page. On synthetic series it finds a drift ending at frame
  60 at 57, one ending at 20 at 18, noise at 0, and reports the ceiling for a ramp
  that never settles.

  Every scored scalar is then cut from that frame, which is tested by holding a
  ligand in the pocket for the first half of a run and 25 A away for the second:
  the occupancy has to come out at 0, because the window is the half where the
  ligand is gone. Averaged over the whole run it would read 0.5, which is the
  error that test exists to catch.
- **Three superpositions, each against its own reference copy.** The pose RMSD fits
  the pocket Calpha and measures the ligand, the in-place RMSD fits the ligand and
  measures the ligand, and the receptor RMSD fits and measures the receptor
  Calpha. The first two together are what separate the ligand moving from the
  ligand changing shape.

  Each takes a second copy as its reference rather than passing the trajectory as
  its own. MDTraj's `superpose` rewrites coordinates in place, so a trajectory
  fitted against itself moves frame 0 first and then fits every later frame against
  the moved version: a rigid 2 A translation came out as 1.88 A that way. Wrong,
  plausible, and one unit test away from shipping. The test asserts the exact
  number on a jitter-free trajectory, where a unit or reference error has nowhere
  to hide.
- **Centroid drift is measured from the ligand's own frame-0 position.** Not from
  the pocket centroid, which is the average of the residues lining the site and
  sits several angstrom from the ligand centre even in a perfectly held pose:
  measuring against it adds a constant offset, and the 5 A hard gate then fires on
  ligands that have not moved. The offset is still reported, as a diagnostic.
- **The periodic-boundary check is measured, not defaulted.** A ligand wrapped
  across a boundary sits a box length from its receptor in the stored coordinates,
  and then every distance spanning the two is wrong by that much.  `pbc_ok: false`
  is not a failure, but it makes every distance in the file suspect and the report
  says so, so it is a measurement: the largest ligand-centroid step between
  consecutive frames, against a 10 A threshold. The contacts are computed before
  any superposition for the same family of reason: `superpose` rotates coordinates
  and leaves the unitcell vectors alone, so a periodic distance taken afterwards
  applies stale box vectors to rotated coordinates.
- **`series/<cid>.npz` and the trajectory the viewer can step through.** The
  per-frame series covers the **whole** trajectory rather than the scored window,
  because DeepFlex draws the discarded head in a different colour and that is how a
  reader sees what was thrown away. A test asserts the scalar and the series agree:
  the median of `pose_rmsd_a` over the window has to equal `ligand_rmsd_median`, or
  the chart and the number beside it describe different trajectories.
- **One multi-model PDB, not one file per frame.** Around a hundred models, cropped
  to the pocket and written from the pocket-fitted copy so the receptor holds still
  and what moves is the ligand. The viewer steps the model index, which is a
  parameter update; a file per frame makes it tear the structure down and rebuild
  it, so the canvas blanks and the camera jumps on every step. Without the member
  the viewer shows a single frame and the trajectory control is absent.
- **A time axis that is checked rather than assumed.** An `.xtc` carries its own
  frame times and they are used. Some writers leave them at zero, and a time axis
  of all zeros collapses every chart to a single point without anything raising, so
  the fallback is the configured frame spacing and the last resort is the frame
  index. Both fallbacks add a note saying which, because an axis labelled
  nanoseconds that is really frame numbers is an error a reader cannot see.
- **An interaction fingerprint, typed by geometry rather than by chemistry.**
  Polar contact, ionic and hydrophobic, from heavy-atom distances. DeepFlex uses
  ProLIF for this, which needs bond orders, formal charges and ring perception;
  FlexAppeal has none of those for a ligand at this point, because the SDF the
  bundle fetched stays on the user's machine and never enters the package.
  **Nothing is claimed as a pi-stack**, and the interaction record says so in a
  note rather than leaving a reader to assume the vocabulary matches a `.dfx`.

  No angle term on the hydrogen bond, deliberately: a decimated payload may carry
  no hydrogens at all, and a criterion that silently changes definition depending
  on which source the trajectory came from is worse than one that is coarse
  everywhere.

  Checked against the chemistry: benzamidine's amidinium to trypsin's ASP189 comes
  out as both a polar contact and ionic, at 100% occupancy, present at frame 0.
  That is the interaction the complex is known for, so a fingerprint that misses it
  is measuring the wrong atoms whatever else its numbers look like.
- **Interactions that never happened are not columns.** The candidate set is every
  pair the topology makes geometrically possible, which for a pocket is mostly
  pairs that never come within their cutoff. Keeping them cost two numbers:
  `hbond_occupancy_mean` read 0.07 instead of 0.44, because it was the mean over
  dozens of zeros rather than over the hydrogen bonds that exist, and the occupancy
  table carried eighty-two rows claiming contacts at 0%.
- **`summary/` and `barcode/`.** The occupancy table records a mean lifetime beside
  each occupancy, which is what separates one long contact from a hundred flickers
  at the same fraction. The barcode is capped at twenty-eight rows, busiest first,
  and four hundred columns, because twenty-eight rows is already a tall panel and
  the tail is contacts present in a handful of frames.
- **The energy gap, named rather than papered over.** With no `--energies` file the
  command now says in one warning that FlexAppeal cannot compute MM/GBSA or an
  interaction energy, why (every frame needs rescoring in implicit solvent from a
  rebuilt parameterised system, and nothing in the pipeline produces one), and
  exactly what is lost: the energy panel, the MM/GBSA column of the league table,
  and the 0.15 energy weight, which DeepFlex renormalises away and records as
  pending. A partly supplied cohort is flagged too, because a percentile needs a
  cohort and a complex without an energy takes no percentile rather than a last
  place it did not earn.
- **A per-term energy decomposition, carried through.** `terms` in the energies
  file fills the energy breakdown panel from somebody else's `MMPBSA.py` output.
  Carried rather than interpreted: these are not FlexAppeal's numbers and it has no
  basis for arguing with them beyond refusing an absolute block.
- **`S_conv`: the four-block divergence, with a bin count that is not fixed.** A
  histogram divergence between two finite samples is biased upward by roughly the
  bin count over the sample size, so at a fixed twenty bins a **stationary**
  400-frame window measures 0.086 against a threshold whose full marks are 0.05: a
  perfectly converged short run scored as a diverging one, with nothing about the
  output looking wrong. Holding samples-per-bin constant instead keeps the bias
  constant, and the divisor is chosen so a long run lands on twenty bins, which is
  the resolution the spec's thresholds were calibrated against. Measured on
  stationary noise the divergence now stays under 0.02 from 60 frames to 1,800.

  It is the divergence and not the Jensen-Shannon distance, which is its square
  root: a 0.06 divergence would read 0.24 as a distance, the difference between
  full marks and most of the way to none.
- **Pose clusters, and the alternative-pose flag.** Average-linkage clustering on
  pairwise ligand RMSD, cut at 2 A, with one pocket PDB per cluster centroid so the
  panel can link each row to a real frame. The flag compares the dominant cluster
  of each half of the scored window and requires the second half's winner to
  actually dominate it, because comparing the first and last frame fires on
  ordinary thermal motion and the flag caps the grade at B.
- **`S_conf`: MMFF94 strain in excess of the frame-0 pose, never the absolute
  number.** The coordinates come from an OpenFF or GAFF force field and are being
  scored by MMFF94, and a 298 K snapshot of a drug-sized ligand carries tens of
  kcal/mol in its internal modes. Both offsets are constant per ligand and cancel
  when frame 0 is subtracted; an absolute MMFF number is not comparable between
  ligands and would be scored as if it were.

  The strain itself is torsional rather than a single point: minimise with every
  rotatable torsion restrained to its observed value, minimise again with nothing
  restrained, take the difference. That is the cost of holding the observed
  torsions, and it cancels the bond and angle terms inside each frame as well as
  between frames, which matters because the MD constrained its bonds under a
  different force field.
- **The ligand's bond orders are inferred from the trajectory, not read from the
  SDF.** The chemical definition the bundle fetched is still beside a finished run
  and is deliberately not used. A version 1 `.fxa` carries no SDF, so a strain that
  needed one would be a column that silently means different things depending on
  where the trajectory came from; and with no code that opens a chemistry file,
  there is no path by which one could reach the package, which section 8 of the
  spec bans outright. It needs hydrogens, and a payload without them leaves the
  keys out and the sub-score pending. Verified on the committed trypsin run:
  benzamidine comes back with two rotatable torsions and no flips, which is what a
  near-rigid ligand should measure.
- **`rdkit` declared in `pixi.toml`.** Not a new dependency: it already arrives
  with `openff-toolkit`, so the solve does not change. Named for the same reason as
  `rich` and `psutil` beside it, and imported lazily, so a build without it leaves
  the conformation sub-score pending rather than failing.
- **`panels/`: the receptor's own flexibility.** Calpha RMSD, radius of gyration,
  per-residue RMSF, a frame-against-frame RMSD matrix, a residue cross-correlation
  matrix and the first two principal components, all from the receptor-fitted copy
  so they describe internal motion rather than the whole molecule tumbling.

  Both matrices are reduced to at most 500 on a side before being written, because
  a 2,000-frame 2D RMSD is four million floats and the panel is a few hundred
  pixels wide. The cross-correlation decimates *residues* rather than cropping the
  matrix, so a receptor with more than 500 of them still gets a panel spanning the
  whole chain rather than its first quarter. The components come from an SVD of the
  centred coordinates rather than an eigendecomposition of their covariance, which
  for a 200-residue receptor would build a 600 by 600 matrix to extract two
  vectors.

  `rmsf_resid` and `rmsf_a` share a length and an order, and a test asserts the
  numbering agrees with `pocket/<cid>.json`: the report marks pocket residues on
  this plot by looking their numbers up, so a disagreement marks the wrong
  residues or none, and nothing about either file on its own would look wrong.
  Trypsin's radius of gyration comes out at 16.3 A, which is the published value
  and reads 1.6 if the unit conversion is missed.
- **`annotations/`, behind `--annotate`.** A receptor's UniProt accession, entry
  name and protein name, from the RCSB entry's first polymer entity. Off by
  default because it is the only part of this format that needs the internet and a
  packaging run should work on a train; a test asserts nothing reaches the network
  without the flag. A failure is caught per receptor and costs that page its
  accession rather than aborting the package. No consensus numbering:
  `resid_to_generic` is a GPCR convention built from an alignment to a reference
  receptor, and FlexAppeal has neither, so `n_generic_mapped: 0` says the mapping
  found nothing rather than leaving the field out and implying it was never tried.
- **`depictions/`: a 2D ligand diagram, shaded by per-atom mobility.** Blue to
  orange, two hues rather than a rainbow: the quantity is one-dimensional and a
  reader has to tell which end is which without a key, and a red-green ramp is the
  one choice that fails for a common form of colour blindness. The RMSF travels as
  an RDKit atom property so it survives `RemoveHs` attached to the atom it belongs
  to; mapping by position afterwards would map onto a renumbered molecule.
  `painted.json` lists only the diagrams that have the shading, because the
  mobility legend is a claim about colours.
- **Both readers now refuse the other's file by its number.** A version 2 package
  has a `metrics/<cid>.json` per complex and no `metrics.json` at all, so the
  version 1 reader rejected one as truncated: a message about a missing member that
  said nothing about the real problem. In the other direction a version 1 file
  carries `traj.xtc` quite legitimately, and the version 2 reader rejected it for a
  banned member rather than for its version. Both now check the declared version
  first, after the hostile-archive checks and before anything that assumes a
  layout, and the version 1 message says where a package does belong. No existing
  version 1 test changed.
- **The command in the README**, with the reason a study has to be a package
  rather than a bigger single-system file.
- **`n_exit_events` is counted over the scored window**, as section 4 asks
  ("within the window"), not over the whole trajectory. An exit in the
  equilibration head is a frame the file itself says was discarded, so counting it
  would put a departure on the page that nothing scored. The test asserts the
  scalar against the series rather than against a fixed number, because which
  frame the detector picks is its business: what has to hold is that the count
  starts at `equil_start_frame`.

### Added

- **The structure viewer in DeepFlex's neon presentation.** Near-black ground,
  the ligand's carbons pinned to one vivid colour, emissive sticks and a
  translucent shell that Mol\*'s bloom pass turns into a glow, with bloom in
  `emissive` mode so only the ligand blooms and the protein behind it stays
  readable. The constants are DeepFlex's, which got them from BINMAN, and each
  is tuning somebody already paid for: carbons are pinned because a
  chain-coloured ligand can come out the same hue as its own shell and vanish
  into the glow, and the sticks are thickened because thinner ones read as a
  cloud of unconnected dots with the bonds present and invisible. Every step is
  guarded, so a viewer that cannot apply it shows the default presentation
  rather than nothing.
- **The 2D ligand diagram, over the viewer, with a mobility legend.** Each atom
  is shaded by how far it moved, which is the question three dimensions cannot
  answer: a mobile tail is hidden behind whatever it is moving in front of. Six
  swatches rather than a gradient, because the shading is banded and a reader can
  match an atom to a band and recover a number from it. The legend reads its
  colours from the module that paints the diagram, so it cannot come to describe
  colours that are not on it.

  The bands are **absolute**, not per-ligand. A min-to-max normalisation was the
  first version and it is actively misleading: it paints the stillest ligand in a
  study exactly like the twitchiest one, so the colour says where an atom sits
  within its own molecule and nothing about how far it moved. Fixed half-angstrom
  edges to 3 Å, matching DeepFlex, so a diagram here and one in a graded report
  mean the same thing. Benzamidine in trypsin comes out entirely in the lowest
  band, which is what a tightly held ligand should look like.

  Drawn light on the viewer's own ground rather than as a white card floating on
  it. One bond and one atom label came out pure black on near-black, because
  `RemoveHs` keeps the hydrogens on a heteroatom and the palette default did not
  reach them.
- **A frame slider, linked to every time plot.** The viewer's blurb already
  promised one; there wasn't one. It replaces both halves of what was there:
  Mol\*'s own model-index animation, and a 250 ms poll that watched the frame
  change so the plot cursor could follow. Two things driving one model index
  fight, and the slider's position then disagrees with what is drawn; and a poll
  could only read the frame, never set it, so there was no way to go to a frame a
  plot made you curious about.

  Stepping a model is a parameter update, so dragging drives it directly and the
  camera is never touched: zoom and orientation stay where they were left.
  Loading a frame as a structure instead is what makes a scrub flash and jump.
  Dragging pauses the loop, because a loop that keeps advancing under the hand on
  the control is fighting it.
- **The binding panel ids are namespaced**, because the id becomes a DOM id and
  these panels share a page with the whole-system ones. `clusters` collided:
  conformational clustering above, pose clustering below, two elements with
  `id="plot-clusters"`. Both still drew, because the bootstrap walks nodes rather
  than ids, which is exactly why it went unnoticed; the document was invalid and
  anything reaching for one by id found whichever came first. A test now asserts
  no plot id appears twice on a page, so the next collision fails rather than
  renders.
- **The cursor is appended to each panel's own shapes, not substituted for
  them.** Several binding panels carry their own: the shaded equilibration head,
  the 5 Å drift gate, the contact threshold. Replacing the list wiped all of them
  the first time the slider moved.

### Fixed

- **A ligand in the file is now simulated by default.** Three things stood between
  the defaults and a runnable complex, all found by driving the real entry points
  rather than reading them.

  `./FlexAppeal.py bundle 3ptb.pdb` set `chains` and never `keep_heteroatoms`, so
  it silently produced an apo trypsin run and dropped the benzamidine.

  **No bundle carried a ligand's chemistry at all.** `bundle.collect_ligands` had
  no caller outside the test suite, so `ligands/<NAME>.sdf` was never in a bundle.
  A PDB `HETATM` record has no bond orders and the run script transfers them from
  that member, so a kept ligand warned "no chemical definition bundled" and then
  failed at system creation with a missing residue template, several minutes into
  somebody's run. Both entry points now fetch them at prep time, where the network
  already is, so a fetch failure is a clear refusal there rather than a crash hours
  later on a machine the user is waiting on.

  And the web app's auto-selection filtered the heteroatom list by a prefix match
  on help text, which could not see that a haem has an iron in it: haemoglobin
  defaulted into a configuration `collect_ligands` then correctly refuses.
- **`structure.default_heteroatoms` owns which heteroatoms count**, for both entry
  points, by category rather than by matching names. Real ligands and cofactors are
  kept. Crystallisation additives are not, because glycerol, PEG and sulfate are
  there because of how the crystal was grown and parameterising one costs an
  AM1-BCC charging run for something nobody is studying. Ions are not, because the
  protein force field carries its own ion parameters and listing one as a ligand
  would imply it needed a template. Anything bearing a metal a small-molecule force
  field cannot type is not, because such a system cannot be parameterised at all,
  and defaulting to a configuration that cannot build is worse than defaulting to
  one that is merely incomplete.

  Every omission is reported rather than silent, so a reader is told their cofactor
  was left out and why instead of noticing later. All of it stays tickable, and an
  explicit `--config` still wins on the command line.
- **A ligand-less trajectory selection is named as the cause.** A run written with
  `traj_selection` of `protein`, `protein_heavy`, `backbone` or `ca` has no ligand
  coordinates, and packaging it failed with "no ligand in the trajectory": true,
  and silent about why. The packager now says which selection the run used and what
  to use instead.

### Added

- **A `Study panels` section in the README**, naming each panel a version 2
  package populates and the measurements behind it, plus the one gap (the binding
  energies) and how `--energies` fills it.

### Changed

- **A default-length run no longer needs gigabytes to package.** FlexAppeal's
  default job is 100 ns at 10 ps, which is 10,000 frames, and three
  superpositions each hold a copy of the coordinates: measured at 1,022 MB for
  2,000 frames of a 3,000-atom system, so five gigabytes for a default run and
  tens for a membrane one, where `protein_ligand` keeps the lipids as well.

  Three changes, in order of how much each bought. Sources are now reduced to
  protein plus the kept heteroatoms with the selection and the stride applied by
  `md.load` itself, so the large array never exists; nothing measured here reads a
  water, an ion or a lipid, and section 6.1 bans all three from the crop anyway.
  Frames are capped at 2,000, with the frame spacing in the manifest multiplied by
  whatever was skipped, because every scalar is a statistic over the scored window
  and 2,000 frames agree with 10,000 to more digits than anybody reads. And a
  superposition reference is now one frame rather than a whole second copy:
  `superpose` reads only `reference.xyz[frame]`, and the fitted RMSDs are
  bit-identical either way, which was checked rather than assumed.

  Peak is 805 MB for 2,000 frames and bounded there for a run of any length.
- **Lipids are stripped from the pocket crop and flagged by the conformance
  check.** Section 6.1 bans them alongside the waters and the counter-ions, and
  check 7 names only the latter two; a bilayer in a binding-site view is noise
  that hides the site. Both spellings are listed, because a PDB residue name is
  three characters wide so POPC is written and read back as POP, and selecting on
  the configured name alone matches nothing silently.

## 2026-07-27

### Added

- **The ASCII banner on the Example tab**, so the wordmark now leads all three
  rather than two. It sits in `_example_intro.html` rather than the results
  branch that tab shares with Analysis, because a real upload should open on the
  run's own header instead.

### Changed

- **The clip sits on white, not cream.** Its background was `#faf9f5`, a shade
  off the panel it sits in. Keyed to the panel's own white by RGB distance
  (`colorkey`) rather than by a levels curve: the yellow ribbon's red channel is
  244 against the background's 250, close enough that any brightness-based
  remap would have shifted the artwork too.
- **The clip panel is as tall as the panel beside it**, with the video absorbing
  the difference and the caption pinned underneath it. `height: auto` on the
  video is load-bearing and its absence is silent -- the `height` attribute is a
  presentation hint that otherwise pins the element to 552px, beating both the
  intrinsic ratio and an explicit `aspect-ratio`. Invisible on desktop, where
  flex-grow resizes it anyway; it only showed up stacked, as a tall column of
  white. Measured in the browser rather than reasoned about, after two wrong
  guesses. Tested.

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
