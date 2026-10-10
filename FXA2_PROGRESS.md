# `.fxa` version 2: progress

Running log for the version 2 package writer, specified in `DEEPFLEX_FXA2_SPEC.md`
and driven by `FXA2_BUILD_PROMPT.md`. Rewritten at every step boundary, so an
interrupted session can be resumed from this file alone.

**Status: all ten steps of section 12 are built. Nothing is pushed.**

## What the validator reports

`pixi run python tools/validate_package.py <pkg>.fxa --verbose`, run from the
DeepFlex checkout, on a four-complex package built by `flexappeal package` from
the committed trypsin fixture with `--energies` supplied:

| Group | |
|---|---|
| grade, pose, occupancy, interactions, energy, convergence, conformation, receptor panels, ranking | **100%** |
| payload | **80%** |

> every sub-score populated: this file would produce a full report.

**The one short group, and why.** `payload` is missing `pocket_traj_path` alone.
DeepFlex fills that column from `structures/<cid>_pocket.xtc`, and section 8 of the
spec bans `.xtc` members outright. There is no way to fill it without writing a
member the spec forbids, so it is left NULL. The fix belongs on the DeepFlex side:
the column should accept the multi-model `_scrub.pdb` that section 6.2 asks for,
which this writer does produce and which the viewer already loads by name.
Divergence 5 below.

Without `--energies` the energy group is 0% and grade 62.5%. That is the one
measurement FlexAppeal has no way to produce, documented under **What FlexAppeal
cannot compute** below.

### The cohort, checked rather than assumed

A package of two receptors by two ligands ingests and grades, and the rows show
the comparative columns working: grades A and A-, `pending_subscores` empty on
every complex, `gb_cohort` equal to the receptor rather than to the package (which
is the ruling the spec calls out), `gb_cohort_n` 2, and ranks 1 and 2 within each
receptor.

## The shape of the work

Three new modules, so the format, the measurements and the orchestration can be
tested apart from each other:

| Module | Job |
|---|---|
| `flexappeal/fxa2.py` | the format: identifiers, archive limits, the chemistry bans, `Package`, `read()`, `conformance()` |
| `flexappeal/measure.py` | one complex's measurements from an MDTraj trajectory, in angstrom |
| `flexappeal/package.py` | `flexappeal package`: collect finished runs, derive ids, assemble the file |

Plus `./FlexAppeal.py package` and `tests/test_fxa2.py`.

Version 1 is untouched. `flexappeal/fxa.py` keeps `FXA_VERSION = 1` and its own
reader; version 2 is a separate module because it is a different kind of file, not
a bigger version of the same one.

## Build order

Section 12 of the spec, which is also the commit boundaries.

- [x] 1. Writer skeleton: manifest, member layout, the zip, the limits, round trip.
- [x] 2. `metrics/`, `pocket/`, one pocket PDB, one complex. **Validated.**
- [x] 3. The scored window, `S_pose`, `S_occ`.
- [x] 4. `series/` and `traj/<cid>_scrub.pdb`.
- [x] 5. `S_int`, `summary/`, `barcode/`.
- [x] 6. `S_E`, the energies.
- [x] 7. `S_conv`, the block split and the pose clusters.
- [x] 8. `S_conf`, the strain excess and the torsions.
- [x] 9. `panels/`, the receptor matrices.
- [x] 10. `annotations/` and `depictions/`.

## Does a default job produce all of this?

Checked against the option registry and driven end to end, rather than assumed.

| Default | Value | Consequence |
|---|---|---|
| `keep_heteroatoms` | **the file's real ligands and cofactors** | a complex, with no option touched |
| `has_ligands` | **derived from the above** | the small-molecule force field is wired in |
| `traj_selection` | `protein_ligand` | the ligand **and its hydrogens** are in the trajectory, so `S_conf` works |
| `enforce_periodic_box` | `True` | the ligand can be wrapped away from its receptor, so every source is re-imaged before anything is measured |
| `production_duration` / `traj_interval` | 100 ns / 10 ps | 10,000 frames, strided to `MAX_MEASURE_FRAMES` |
| `payload_tier` | `standard` | 500 frames, protein plus the kept heteroatoms, hydrogens kept: enough for every panel |
| `traj_format` | `xtc` | read by MDTraj, as are `dcd` and `hdf5` |

`./FlexAppeal.py bundle tests/fixtures/3ptb.pdb` with nothing else set now produces
`keep_heteroatoms: ["BEN"]`, `has_ligands: true`, `LIGAND_RESIDUES = ["BEN"]`,
`traj_selection: protein_ligand` and a `ligands/BEN.sdf` member. That is a complete
receptor-ligand run, and a package built from one fills every panel except the two
gaps named below.

### What that needed fixing first

Three things stood between the defaults and a packageable run, all found by driving
the real entry points rather than reading them:

1. **The command line never kept any heteroatom.** `cmd_bundle` set `chains` and
   not `keep_heteroatoms`, so `./FlexAppeal.py bundle 3ptb.pdb` silently produced
   an apo trypsin run and dropped the benzamidine.
2. **No bundle carried a ligand's chemistry.** `bundle.collect_ligands` had no
   caller outside the test suite, so `ligands/<NAME>.sdf` was never in a bundle. A
   PDB `HETATM` record has no bond orders and the run script transfers them from
   that member, so a kept ligand warned "no chemical definition bundled" and then
   failed at system creation with a missing residue template, several minutes into
   somebody's run. Both entry points now fetch them at prep time, where the
   network already is, and a fetch failure is a clear refusal there rather than a
   crash hours later on the user's machine.
3. **The web app's auto-selection could not see a metal.** It filtered the
   heteroatom list by a prefix match on help text (`"crystallisation"`, `"ion"`),
   which let a haem through; `collect_ligands` then refuses it, correctly, so
   haemoglobin defaulted into a configuration that cannot build.
   `structure.default_heteroatoms` now owns the decision for both entry points.

### Which heteroatoms count

Real ligands and cofactors, and nothing else. Decided by category rather than by
matching names:

| Left out | Why |
|---|---|
| crystallisation additives | glycerol, PEG and sulfate are there because of how the crystal was grown, and parameterising one costs an AM1-BCC charging run for something nobody is studying |
| ions | the protein force field carries its own ion parameters, so an ion needs no small-molecule template and listing it as one would imply it did |
| anything with an unparameterisable metal | GAFF, OpenFF and espaloma have no transition-metal parameters, so a haem cannot produce a runnable system at all. Defaulting to a configuration that cannot build is worse than defaulting to one that is merely incomplete |
| waters and membrane markers | removed before simulation anyway |

Each omission is reported rather than silent, so a reader is told their cofactor
was left out and why instead of noticing afterwards. Everything stays tickable:
this is a default, not a policy, and an explicit `--config` wins on the command
line.

**An apo run still cannot be a complex.** `flexappeal package` skips one with a
reason and refuses a package where every source is apo. That is inherent rather
than a gap.

**A ligand-less trajectory selection is now named as the cause.** A run written
with `traj_selection` of `protein`, `protein_heavy`, `backbone` or `ca` has no
ligand coordinates at all, and the failure used to arrive as "no ligand in the
trajectory": true, and silent about why. The packager now says which selection the
run used and what to use instead.

## What a default job costs to package

Three superpositions each hold a copy of the coordinates, so 10,000 frames of a
3,000-atom system came to roughly five gigabytes, and a membrane system (where
`protein_ligand` keeps the lipids as well) to tens. Measured, not estimated: 2,000
frames of the trypsin fixture peaked at 1,022 MB before the fix.

- **Sources are reduced to protein plus the kept heteroatoms**, with the stride
  and the atom selection applied by `md.load` itself, so the large array never
  exists. Nothing measured here reads a water, an ion or a lipid, and section 6.1
  bans all three from the crop anyway.
- **`MAX_MEASURE_FRAMES = 2000`**, with the frame spacing in the manifest
  multiplied by whatever was skipped. Every scalar is a statistic over the scored
  window and agrees between 2,000 frames and 10,000 to more digits than anybody
  reads; 2,000 is also the length of DeepFlex's own runs.
- **A one-frame superposition reference** rather than a whole second copy.
  `superpose` reads only `reference.xyz[frame]`, and the fitted RMSDs are
  bit-identical either way.

Peak is now 805 MB for 2,000 frames and bounded there for a run of any length.

## Front end

**No changes needed.** Packaging is a command, not a page: the web app prepares
bundles and reads version 1 files, and neither path is involved. The one place the
app can meet a version 2 file is the Analysis upload, and that refuses it by
version with a message naming `flexappeal package` instead of reporting a missing
member, with a `tests/test_webapp.py` case driving it through the real route.

Worth offering later, not needed now: a line on the Analysis tab saying what a
version 2 package is and where it goes.

## What is left, and what it would cost

Nothing blocks a full report. Three things are known gaps rather than bugs:

1. **The binding energies**, below: the one measurement FlexAppeal cannot make.
2. **`pocket_traj_path`**: unfillable without writing an `.xtc` member that
   section 8 forbids.
3. **The tracked-distance panel** (divergence 16): no FlexAppeal option
   corresponds to it, so inventing one would put a number on that panel nobody
   asked for.

Worth raising with the DeepFlex side, in order of how much each costs a reader:
divergence 5 (`pocket_traj_path`), divergence 4 (the two layer manifests), then
1 to 3, which cost nothing now because this writer works around them.

## The benzene example, and what benzene cannot show

`examples/lysozyme_benzene_100ns` is T4 lysozyme L99A with benzene in its buried
cavity (PDB 181L), the canonical cavity-binding model system. Verified end to end
against the live run directory at 19% of the production run: `flexappeal package`
reads it, DeepFlex accepts it as a graded run, and the numbers are the ones this
system should give.

| Measurement | Value | Reading |
|---|---|---|
| `ligand_rmsd_median` | 1.97 A | benzene rattling in the cavity |
| `inplace_rmsd_median` | 0.06 A | rigid: a flat symmetric ring has nothing to flex |
| `in_pocket_fraction` | 1.00 | never leaves a buried site, which is the point of the system |
| `com_drift_mean_a` | 0.61 A | |
| `n_rotatable_bonds` | 0 | |
| `radgyr_mean` | 16.7 A | T4 lysozyme, 162 residues |
| pocket | 15 of 54 residues | the cavity, out of the 10 A crop shell |

**Three columns benzene structurally cannot fill**, and each is the "omit, never
zero" rule doing its job rather than a gap in the writer:

1. `hbond_occupancy_mean`. Benzene has no N, O or S, so it can form no hydrogen
   bond, so the mean occupancy of its hydrogen bonds is undefined rather than
   zero. The interactions group reads 80%.
2. `flips_per_ns_per_torsion`. No rotatable bonds, so the torsion-flip term is
   vacuous and `S_conf` rests on the strain alone. The conformation group reads
   67%, and the torsions panel does not draw at all: 11 of the 12 binding panels
   rather than 12.
3. `mmpbsa_gb_dg` and the rest of the energy block, which is the documented gap
   below and nothing to do with benzene.

So the example demonstrates the pose, occupancy, cavity, convergence and receptor
panels on a system the literature has an opinion about, and cannot demonstrate
hydrogen bonding or torsional flexibility, because this ligand has neither. A
second complex with a polar, flexible ligand would fill those two and give the
package a cohort for `S_E`; that is a decision for Marc rather than one to take
while he is asleep.

## Decisions taken

### D1. A separate `package` command, not a change to how jobs run

`flexappeal package --out study.fxa <source> ...`, taking either a finished run
directory or a version 1 `.fxa`. Marc's recommendation, adopted unchanged: a
version 2 file needs a cohort, FlexAppeal runs one system per job, and making the
package a deliberate act lets a study be assembled from runs done weeks apart.
Nothing about an existing job changes and no version 1 file is altered.

Both source kinds are accepted because they cost almost nothing apart: a run
directory carries the complete trajectory and the ligand's chemical definition, and
a version 1 `.fxa` carries whatever its payload tier packed. The run directory is
the better source and the command says so when handed the lesser one.

### D2. Receptor token, in the order the spec gives

`config.pdb_id`, then the prepared structure's filename stem, then the job name.
Uppercased and reduced to `[A-Z0-9_-]`. A collision between two different receptors
landing on the same token is an error rather than a silent merge.

### D3. The ligand id and the salt

`HMAC-SHA256(salt, compound_name)[:16]`, 32 fresh random bytes per file, never
recorded in it. The compound name is the PDB chemical component code the job kept
(`config.keep_heteroatoms`), which is the only compound identity a FlexAppeal job
has. It travels in clear in `manifest.ligand_names`, as section 2 allows.

### D4. The pocket, raw and filtered

Section 6.5 asks for both and does not say how to derive either, and the filtered
set is load-bearing: it is the superposition frame for the ligand RMSD, so a
selection that includes a flexible loop or a terminus makes the fit wander and
every ligand RMSD measured against it is wrong in a way that looks exactly like the
ligand moving.

- `pocket_residues_raw` is section 6.1's crop shell: every whole residue with any
  atom within **10 A** of the ligand at frame 0. This is what the pocket PDB holds.
- `pocket_residues_filtered` is residues with a heavy atom within **6 A**, less the
  first and last **5** residues of every chain. Under **12** survivors the shell is
  widened to **8 A** and `used_fallback` records it.

Those four numbers are DeepFlex's own, taken from its pocket layer rather than
chosen here. A pocket defined differently from the cohort it is ranked against is
not comparable to it, and the same goes for the 5 A contact cutoff, the quarter of
native contacts for occupancy and the 2 A pose-cluster linkage cutoff used later.

`added_vs_raw` is therefore always empty for a FlexAppeal package: the filter is a
subset of the shell, so it can only drop. The key is still written, because an
empty list is a measurement and an absent key is not.

### D5. `kind` is decided from what was measured, not from what was intended

Section 9 refuses a `graded` file with required inputs missing, and says
`flexibility` always works and is never wrong, only less. So the writer checks every
analysable complex against section 4's required list after measuring, and downgrades
the whole file if any is short, naming what was missing. A run that loses a
measurement downgrades itself rather than claiming a grade it cannot support.

### D6. Residue numbers are not unique, and are kept that way

3PTB uses chymotrypsin numbering with insertion codes, which MDTraj drops, so
residue 184 is both a GLY and a TYR in the trypsin fixture. Pocket lists and
`rmsf_resid` report `resSeq` because that is the number a reader sees in a viewer
and the number DeepFlex looks up to mark pocket residues on the RMSF plot.
Duplicates are kept rather than merged: dropping one would silently shorten the
pocket. Atom selection is always by residue *index*, never by number, so a shared
number cannot pull in the wrong residue.

### D7. `detect_equilibration` without pymbar

Chodera's rule: take the truncation that leaves the most effectively uncorrelated
samples, `(N - t) / g(t)`, where `g` is the statistical inefficiency. That is what
`pymbar.timeseries.detect_equilibration` maximises, so this is the same quantity
rather than a cheaper substitute, and pymbar is not in `pixi.toml`.

Two things it does differently, both deliberate:

- **About fifty candidate truncations, not every frame.** The curve is smooth and
  the search is otherwise cubic in the frame count, which on a 2,000-frame run
  costs minutes to produce one frame index.
- **Searched unconstrained, then clamped to half the run.** Clamping a search that
  was already bounded cannot tell a run that equilibrated just under the cap from
  one that never equilibrated at all, and those two deserve different words on the
  page. Over the cap the method reads `detect_equilibration->ceiling`, which
  DeepFlex renders as `capped` (it looks for `ceiling` as a substring).

Measured on synthetic series: a drift ending at frame 60 is detected at 57, one
ending at frame 20 at 18, pure noise at 0, and a never-settling ramp reports the
ceiling.

### D8a. `n_exit_events` is counted over the window

Section 4 says "times the ligand left and did not return **within the window**",
so it is counted from `equil_start_frame` onward like every other scored
quantity. DeepFlex's binding layer counts the transitions over the whole
trajectory, which includes the equilibration head: an exit there is a frame the
file itself says was discarded, so counting it puts a departure on the page that
nothing scored. Same shape of difference as D8 below, and the same remark
applies: it affects the displayed count and not the sub-score.

### D8. `n_native_contacts` counts residues, as the spec says

Section 4 defines it as "receptor residues within the contact cutoff of the ligand
at frame 0", and that is what is written: a closest-heavy-atom distance per
residue, under the minimum image convention, so a residue counts once however many
of its atoms are near. DeepFlex's own binding layer counts atom *pairs* under the
same column name, which is why its numbers run to the low hundreds where these run
to the low tens. The sub-score is unaffected, because `in_pocket_fraction` is the
fraction of frames retaining a quarter of whatever was counted at frame 0 and S_occ
is a ramp on that fraction rather than a percentile against DeepFlex's cohort. The
displayed count will differ from a `.dfx`-sourced run, and that is worth one line
on the DeepFlex side deciding which definition the column means.

### D9. The time axis, when the trajectory has none

An `.xtc` carries its own frame times and they are the right answer. Some writers
leave them at zero, and a time axis of all zeros collapses every chart to a single
point without anything raising, so the fallback is the configured
`traj_interval` and the last resort is the frame index. Both fallbacks add a note
saying which was used, because an x axis labelled nanoseconds that is really frame
numbers is the kind of error a reader cannot see.

### D10. The interaction fingerprint is typed by geometry, not by chemistry

DeepFlex's own interaction layer uses ProLIF, which needs bond orders, formal
charges and ring perception. FlexAppeal has none of those for a ligand at this
point: the SDF the bundle fetched stays on the user's machine and never enters the
package, and re-deriving chemistry from coordinates alone is how a fingerprint
acquires types it cannot justify. ProLIF would also be a new scientific
dependency, which the build prompt says to stop and ask about.

So three types, all of which heavy-atom geometry supports:

| Type | Criterion |
|---|---|
| `HBond` | ligand N/O/S within **3.5 A** of receptor N/O/S: the standard heavy-atom polar-contact proxy |
| `Ionic` | ligand N within **4.0 A** of an ASP/GLU carboxylate oxygen, or ligand O within 4.0 A of LYS NZ, ARG NE/NH1/NH2 or a HIS ring nitrogen |
| `Hydrophobic` | ligand C within **4.5 A** of a carbon in ALA, VAL, LEU, ILE, PHE, MET, PRO, TRP, CYS or TYR |

**Nothing is claimed as a pi-stack**, and the interaction record says so in a note
rather than leaving a reader to assume the vocabulary matches a `.dfx`. The
occupancies themselves are honest; the type vocabulary is coarser.

No angle term on the hydrogen bond, deliberately. A decimated payload may carry no
hydrogens at all, a ligand's hydrogens have no template in an MDTraj topology, and
a criterion that silently changes definition depending on which source the
trajectory came from is worse than one that is coarse everywhere.

Validated against the chemistry: benzamidine's amidinium to trypsin's ASP189 comes
out as both `HBond` and `Ionic` at 100% occupancy, present at frame 0. That is the
interaction the complex is known for, so a fingerprint that misses it is measuring
the wrong atoms whatever its occupancies look like.

`n_key_interactions` uses DeepFlex's threshold (occupancy at or above **0.3** over
the scored window) rather than the spec's wording ("of those [at frame 0], the ones
you treat as key"), so the two sides agree on which contacts carry the grade.

### D11. Interactions never seen are not columns

The candidate set is every pair the topology makes geometrically possible, and for
a pocket that is mostly pairs which never come within their cutoff. Those are
dropped before anything is computed. Keeping them cost two numbers: with them
`hbond_occupancy_mean` read 0.07 instead of 0.44, because it was the mean over
dozens of zeros rather than over the hydrogen bonds that exist, and the occupancy
table filled with eighty-two rows claiming contacts at 0%. A ProLIF dataframe has
no such columns either, so dropping them is also what keeps these comparable.

The occupancy table drops one more class: a contact seen only in the discarded
head and absent at frame 0, which would be a row at 0% and reads as a contact that
never formed. It is still in the barcode, which covers the whole trajectory, and
that is where a reader should see it.

### D12. The convergence bin count follows the window length

This one is a real trap rather than a tuning choice. A histogram divergence
between two finite samples is biased upward by roughly the bin count over the
sample size. At DeepFlex's fixed twenty bins, a **stationary** 400-frame window
measures 0.086, against a threshold whose full marks are 0.05: a perfectly
converged short run would be scored as a diverging one, and nothing about the
output would look wrong.

So the bin count holds samples-per-bin constant instead (`round(n / 90)`, floored
at 3 and capped at 20). The divisor is chosen so a long run lands on twenty bins,
which is DeepFlex's resolution and therefore what the spec's thresholds were
calibrated against: a stationary 1,800-frame window reads 0.017 either way.
Measured on stationary noise the divergence now stays under 0.02 from 60 frames to
1,800, where a fixed twenty bins ranges from 0.09 to 0.017.

Below about 60 frames the floor of three bins is reached and the bias is
unavoidable. The record then carries a note saying the number reads worse than a
longer run of the same trajectory would, which penalises rather than flatters a
short run.

It is the **divergence**, not the Jensen-Shannon distance. The distance is its
square root, so a 0.06 divergence would read 0.24 as a distance: the difference
between full marks and most of the way to none.

### D13. `changed_pose_cluster` compares halves, and the detector gets there first

The flag caps the grade at B on the DeepFlex side, so a false positive is
expensive. It fires when the dominant cluster of the second half of the scored
window differs from the first half's **and** the new one holds at least half of
that second half.

Worth knowing while testing it: a pose change exactly at the midpoint of a run is
indistinguishable from equilibration finishing there, and `detect_equilibration`
reasonably discards everything before it, leaving one cluster and no flag. That is
correct behaviour rather than a miss, and it is why the test moves the pose three
quarters of the way in.

### D14. The ligand's chemistry is inferred from the trajectory

The MMFF94 strain needs bond orders, and the obvious source is the SDF the bundle
fetched from the chemical component dictionary, which is still beside a finished
run. It is not used. Bond orders are inferred from the trajectory's own
coordinates and hydrogens instead (RDKit's `DetermineBondOrders`, trying a small
range of total charges and taking the first assignment that sanitises). Three
reasons, in order of weight:

1. **It works for both source kinds.** A version 1 `.fxa` carries no SDF, so a
   strain that needed one would be a column that silently means different things
   depending on where the trajectory came from.
2. **No chemistry file is opened anywhere in this writer.** Section 8 is absolute
   about chemistry not reaching the package, and the cheapest way to honour a
   rule is to have no code that could break it.
3. Matching an SDF's atom order onto an MD topology's is its own failure mode, and
   a mismatch produces a plausible wrong number rather than an error.

It needs hydrogens, because the valences are what make the inference possible and
MMFF cannot type a molecule without them either. A payload that dropped them (a
deposited structure, or a selection that stripped them) leaves the keys out and
the sub-score pending. Verified on the committed trypsin run: benzamidine comes
back as `N=C(N)c1ccccc1` with two rotatable torsions.

`rdkit` is now declared in `pixi.toml`. It is not a new dependency: it already
arrives with `openff-toolkit`, so the solve does not change, and the repo's own
comment beside `rich` and `psutil` gives the reason for naming it anyway. The
import is lazy and optional.

### D15. The strain is torsional, computed twice per frame

Not a single-point MMFF energy. Per sampled frame: minimise with every rotatable
torsion restrained to its observed value, minimise again with nothing restrained,
and take the difference. That is the cost of *holding the observed torsions*, and
doing it this way cancels the bond and angle terms inside each frame as well as
between frames, which matters because the MD constrained its bonds under a
different force field. `strain_kcal` is then the median of that over the scored
window less its frame-0 value, which is DeepFlex's own definition.

## Divergences between the spec and the reader that actually ingests these files

Found by reading DeepFlex's `pipeline/ingest.py`, `pipeline/store.py` and
`pipeline/grade.py` against the spec. Each one is a case where a file written exactly
to the spec leaves a panel empty, so each is listed here for the DeepFlex side rather
than worked around silently. In every case this writer does **both**: the spec's own
spelling, and whatever the reader reads today.

1. **`harvest/<cid>.json` is required and the spec does not mention it.**
   `ingest_dfx` refuses a package with no `harvest/` ("this bundle carries no
   analysis"), and reads the binding energies, `box_volume_a3`, `pocket_n_raw` and
   `is_reference` from it rather than from `metrics/`. The spec's section 1 layout
   has no `harvest/` entry at all. DeepFlex's own conformance fixture
   (`tests/test_fxa.py::make_fxa2`) writes the same record to both, which is what
   this writer does.

2. **`is_reference` is required.** `grade.score_package` indexes
   `h["is_reference"]` directly, so a harvest record without it raises. Not in the
   spec. Written as `false`: a FlexAppeal study has no reference compound.

3. **`status: "ok"` is required inside `metrics/<cid>.interactions.json`.**
   `score_package` skips the interaction sub-score unless the record carries it. Not
   in the spec, which lists only the seven measurements and `complex_id`.

4. **Two layer manifests the spec never names, and the columns that need them.**
   `store.metric_columns` reads the receptor flexibility rows (`rmsd_ca_mean`,
   `radgyr_mean`, `rmsf_ca_mean`, `n_rmsf_residues`) and the MM/PBSA term
   decomposition from `manifests/panels.jsonl`, and the four strain columns
   (`strain_kcal`, `strain_absolute_kcal`, `strain_t0_kcal`, `n_rotatable_bonds`,
   `n_torsion_flips`, `flips_per_ns_per_torsion`) from
   `manifests/conformer.jsonl`. Section 4 puts all of those in
   `metrics/<cid>.json`, and section 4's own note on `S_conf` says explicitly that
   "DeepFlex reads these four from `metrics/<cid>.json`, which is where this spec
   puts them: you do not need a separate conformer manifest".

   It half does. `ingest.conformer_records` falls back to the metrics record for
   the *grade*, so `s_conf` is scored, but the six `complex` columns stay NULL and
   the conformation panel does not draw: measured at 0% with the metrics keys
   alone and 100% once `manifests/conformer.jsonl` is added. The panels manifest
   has no fallback at all. Both are written, as well as the spec's own spelling in
   `metrics/`, so a DeepFlex that later reads them there finds them. Note the
   spelling: that manifest calls it `n_rotatable`, not `n_rotatable_bonds`.

5. **`pocket_traj_path` wants a file the spec forbids.** The column is filled from
   `structures/<cid>_pocket.xtc`, and section 8 bans `.xtc` members outright. This
   writer does not write one, so that one column stays NULL and the validator's
   `payload` group tops out at 80%. The fix belongs on the DeepFlex side: the
   column should accept the multi-model `_scrub.pdb` that section 6.2 asks for,
   which is what the viewer itself already loads by name.

### D17. Annotations are opt-in, and carry no consensus numbering

`annotations/<receptor>.json` is the only part of this format that needs the
internet, so it is behind `--annotate` and off by default: a packaging run should
work on a train, and a test asserts nothing reaches the network without the flag.
The accession comes from the RCSB entry's first polymer entity, and the entry name
from UniProt; verified against the spec's own example, 3PTB resolves to P00760 /
TRY1_BOVIN. A failure is caught per receptor and costs that page its accession
rather than aborting the package.

The first polymer entity only. An entry with two different chains has two
accessions, and a package keys its annotations by receptor rather than by chain,
so picking one and saying so in a note beats inventing a field for the rest.

`resid_to_generic` is omitted and `n_generic_mapped` is `0`. Consensus numbering is
a GPCR convention built from an alignment to a reference receptor, and FlexAppeal
has neither the reference nor the alignment. The zero says the mapping was
attempted and found nothing to map, which is true; leaving the field out would
imply it was never attempted.

### D18. The depiction is shaded only where there is something to shade

Atoms are coloured by their RMSF over the scored window on a blue-to-orange ramp:
two hues rather than a rainbow, because the quantity is one-dimensional and a
reader has to tell which end is which without a key, and because a red-green ramp
is the one choice that fails for a common form of colour blindness. The RMSF
travels as an RDKit atom property so it survives `RemoveHs` attached to the atom
it belongs to; mapping by position afterwards would map onto a renumbered
molecule.

`depictions/painted.json` lists only the diagrams that carry the shading, because
the mobility legend is a claim about colours and would otherwise be drawn beside a
plain diagram. A ligand whose chemistry could not be inferred gets no diagram at
all rather than an unshaded one.

### D16. No tracked distance

`distance_a`, `distance_select_a`, `tracked_distance_mean` and
`tracked_distance_select_mean` are all omitted. In DeepFlex's pipeline the tracked
distance is a pair of atoms somebody chose for that study; FlexAppeal's option
registry has no equivalent, and inventing one (the ligand centroid to the pocket
centroid, say) would put a number on the tracked-distance panel that nobody asked
for and that is already in `series/` as `pocket_offset_a`. The cost is that one
panel, which the spec's own table names.

## What FlexAppeal cannot compute

**The binding energies, and so `S_E`.** `mmpbsa_gb_dg`, `mmpbsa_gb_sd`,
`lie_total_mean`, `lie_total_sd` and `lie_drift_slope` all need an MM/GBSA or
interaction-energy rescoring of each frame, which means rebuilding the
parameterised system in implicit solvent. FlexAppeal's pipeline does not produce
any of it and adding it would change how a job runs, which the build prompt rules
out without asking.

**Cost, in the spec's own terms:** the energy page (`mmpbsa_gb_dg`,
`gb_percentile`, `gb_cohort_n`, `lie_drift_slope`), the MM/GBSA column of the
league table, and `S_E`, whose 0.15 weight DeepFlex renormalises away and records
as pending. The grade, the ranks and the other five sub-scores are unaffected.

So the keys are omitted, never written as `0`, and the command says what that costs
rather than leaving a reader to find an empty panel.

Where a study has the numbers from elsewhere, `flexappeal package --energies
<file.json>` carries them in and the energy page draws in full. The file is keyed on
the compound name, or `<compound>/<receptor>` where one compound was run against
several, because the complex id is a salted hash that is fresh for every package
file and nobody can write it down in advance:

```json
{
  "benzamidine/3PTB": {
    "mmpbsa_gb_dg": -38.4, "mmpbsa_gb_sd": 2.6,
    "lie_total_mean": -28.9, "lie_total_sd": 1.4, "lie_drift_slope": 0.12,
    "terms": {"gb": {"VDWAALS": [-22.4, 1.9, 0.1], "EEL": [-61.2, 4.1, 0.2]}}
  }
}
```

Only the seven keys section 4 lists under `S_E` are taken, plus the optional
per-term block, so a sidecar cannot be used to smuggle a grade into the package. A
`*_dg` over 1,000 kcal/mol in magnitude is refused: the three absolute MM/PBSA
blocks are each around -100,000 kcal/mol, and a percentile over absolute blocks
ranks the size of the system rather than the affinity of the ligand.

With energies supplied, every group except the receptor panels and
`pocket_traj_path` reaches 100%.

## Two files deliberately left untracked

`DEEPFLEX_FXA2_SPEC.md` and `FXA2_BUILD_PROMPT.md` are both in the working tree and
neither is committed. `flexappeal/fxa2.py` cites the spec by filename, so a reader
of this repository has a dangling reference, and whether to resolve it is a
decision for Marc rather than one to take quietly. The spec is a format document;
the build prompt is a working note. `git add DEEPFLEX_FXA2_SPEC.md` is the whole
fix if the citation should resolve.
