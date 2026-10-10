# `.fxa` version 2: progress

Running log for the version 2 package writer, specified in `DEEPFLEX_FXA2_SPEC.md`
and driven by `FXA2_BUILD_PROMPT.md`. Rewritten at every step boundary, so an
interrupted session can be resumed from this file alone.

**Current step:** 3 of 10 (the scored window, `S_pose`, `S_occ`).
**Validator last reported:** accepted. A real two-receptor package built by
`flexappeal package` from the trypsin fixture ingests as a run with 2 complexes,
2 receptors and 42 pocket residues. Every measurement group is at 0%, which is
correct for step 2: nothing is measured yet. A hand-built package carrying the
step 4 to 8 scalars reaches 100% on pose, occupancy, convergence and ranking, so
the shapes those steps have to fill are settled.

---

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
- [ ] 3. The scored window, `S_pose`, `S_occ`.
- [ ] 4. `series/` and `traj/<cid>_scrub.pdb`.
- [ ] 5. `S_int`, `summary/`, `barcode/`.
- [ ] 6. `S_E`, the energies.
- [ ] 7. `S_conv`, the block split and the pose clusters.
- [ ] 8. `S_conf`, the strain excess and the torsions.
- [ ] 9. `panels/`, the receptor matrices.
- [ ] 10. `annotations/` and `depictions/`.

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

4. **The receptor flexibility summary rows come from `manifests/panels.jsonl`.**
   Section 4 lists `rmsd_ca_mean`, `radgyr_mean`, `rmsf_ca_mean` and
   `n_rmsf_residues` as optional keys in `metrics/<cid>.json`, but
   `store.metric_columns` reads them from a JSON-lines layer manifest the spec never
   names. Written in both places.

5. **`pocket_traj_path` wants a file the spec forbids.** The column is filled from
   `structures/<cid>_pocket.xtc`, and section 8 bans `.xtc` members outright. This
   writer does not write one, so that one column stays NULL and the validator's
   `payload` group tops out at 80%. The fix belongs on the DeepFlex side: the
   column should accept the multi-model `_scrub.pdb` that section 6.2 asks for,
   which is what the viewer itself already loads by name.

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

So the keys are omitted, never written as `0`. Where a study has these numbers from
elsewhere, `flexappeal package --energies <file.json>` carries them in, and then
the energy page draws in full.

## Two files deliberately left untracked

`DEEPFLEX_FXA2_SPEC.md` and `FXA2_BUILD_PROMPT.md` are both in the working tree and
neither is committed. `flexappeal/fxa2.py` cites the spec by filename, so a public
reader currently has a dangling reference, and that is a decision for Marc rather
than one to take quietly:

- The **spec** is a format document and reads clean: benzamidine, 3PTB, TRY1_BOVIN,
  GIPR and GLP1R are all public, and no compound from the unpublished series appears
  in it. Committing it is defensible and would make the citation resolve.
- The **build prompt** says in as many words that DeepFlex is private because it
  discusses an unpublished compound series. That is Marc's to publish or not.

`git add DEEPFLEX_FXA2_SPEC.md` is the whole fix if he wants the citation to resolve.
Nothing is pushed either way.
