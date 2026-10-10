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
