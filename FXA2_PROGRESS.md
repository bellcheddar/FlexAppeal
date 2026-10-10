# `.fxa` version 2: progress

Running log for the version 2 package writer, specified in `DEEPFLEX_FXA2_SPEC.md`
and driven by `FXA2_BUILD_PROMPT.md`. Rewritten at every step boundary, so an
interrupted session can be resumed from this file alone.

**Current step:** 2 of 10 (`metrics/`, `pocket/`, one pocket PDB, one complex).
**Validator last reported:** accepted, on a synthetic package from the step 1
format module alone. pose, occupancy, convergence and ranking 100%; grade 62.5%
(`s_int`, `s_e`, `s_conf` pending); interactions, energy, conformation, receptor
panels and payload 0%, which is what steps 5 to 9 fill.

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
- [ ] 2. `metrics/`, `pocket/`, one pocket PDB, one complex. **Validate here.**
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
