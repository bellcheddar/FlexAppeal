# Progress: reference-state examples

Phase checklist for `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`. Rewritten after every
phase and at least every 30 minutes during long compute (section 0.1 rule 5).

**Last written:** 2026-10-10, end of phase 0.
**Gates open:** none.

| Phase | State | Notes |
|---|---|---|
| **0 Setup** | **done** | Hardware probed, BoltzMaker vendored and pinned, dependencies resolved, both suites green |
| **1 References** | **done** | Ensembles built; P1 17/17, P2 17/17, P3 100/100 at an 11.39 A threshold |
| **2 Start models** | in progress | ABL1 seeds built and trimmed to a shared 233-498 span; β2AR refined models next |
| **3 System build + smoke** | in progress | ABL1 apo solvated (46,721 atoms) and running at 220 ns/day; β2AR membrane next |
| **4 Analysis module** | not started | `flexappeal/refstates/` over the vendored `sse_comparison` |
| **5 Production** | not started | Tier 2, four runs at 50 ns (200 ns), two ligand-bound |
| **6 App** | not started | Two entries in the existing Runs library, trajectory explorer |
| **7 Findings** | not started | `FINDINGS.md`, deliverables |

---

## Phase 0, in detail

| Acceptance criterion | Result |
|---|---|
| `pixi run python -c "import openmm, mdtraj, MDAnalysis, gemmi"` | passes: OpenMM 8.5.2, MDAnalysis 2.10.0, gemmi 0.7.5 |
| Vendored `sse_comparison` tests pass | **50 of 50** |
| OpenMM platform benchmark recorded | **outstanding**, carried into phase 1 |
| BoltzMaker vendored at a recorded SHA | `da939b1be54a27cb8eac7c937fd28e8d1e7d008d` |

The platform benchmark is deliberately deferred. Section 2.2 says to benchmark on
"a 20 ps run of the solvated ABL1 system", and that system does not exist until
phase 3. The 100 ns benzene run finished this morning on OpenCL at 298 ns/day for
25,463 atoms, which is the figure the phase 3 projection starts from; the real
benchmark replaces it as soon as there is a system to measure.

Regression guard: FlexAppeal's own suite was run after every environment change
and stayed at **759 passed**. This environment also serves the deployed
application, so a dependency solve that broke it would be a live incident rather
than a build problem.

---

## What phase 1 has to settle

Two things the spec leaves to run time and that shape everything after them:

1. **Compute budget (Gate G1): settled, with real headroom.** Four runs of 50 ns
   (D007, D008, D009): ABL1 with imatinib from DFG-out and apo from DFG-in, β2AR
   with carazolol from inactive and apo from active. 200 ns in total.

   | Row | Atoms | ns/day | Hours | Runs | Hours |
   |---|---|---|---|---|---|
   | ABL1, water | ~50 k | ~152 | 7.9 | 2 | 16 |
   | β2AR, POPC | ~95 k | ~80 | 15.0 | 2 | 30 |
   | | | | | **4** | **46** |

   **1.9 days against a 96 hour ceiling, with 50 hours spare.** Every run could
   be redone once and the build would still come in at 61 hours, which matters
   because one replica per row means a failed run costs its whole row. Phase 3
   measures both systems and replaces these projections.

2. **β2AR start models.** The spec's primary route is a GPCRdb *refined* model,
   and it says plainly that not every entry has one (2RH1 and 3SN6 both report
   none). The fallback ladder runs through state-specific homology models to
   building one locally. Which rung is reached is discovered, not assumed, and is
   flagged in the app for the reader.


---

## Where the build actually stands (16:1x)

**Done and verified.**

| Phase | State |
|---|---|
| 0 Hardware probe | done; OpenCL single is the only GPU option, 206 against 29.5 ns/day on the CPU |
| 1 References | done; P1 17/17, P2 17/17, P3 100/100 |
| 2 Start models | done; both ABL1 seeds and both β2AR refined models, oriented and validated |
| 3 System building | done; ABL1 46,721 atoms, β2AR 194,032 and 188,xxx with the bilayer |
| 4 Analysis | done; kinase and receptor paths both verified against real structures |
| 6 App | done; fifth tab, per-run explorer, reference browser, 28 tests |
| 5 Production | Tier 1 passed at 59.8 ns/day; Tier 2 running, b2ar_apo through its restraint ladder |
| 7 Deliverables | FINDINGS.md and DECISIONS.md current; the production section is open |

**The four things that cost the most time, all in phase 3.** Each presented as
`Particle coordinate is NaN` and none of them was a timestep problem:

1. The lipid z restraint was a plain `(z - z0)^2`. OpenMM wraps coordinates into
   the box, the lower leaflet sits at negative z and wraps by one box vector, so
   the restraint saw an 11.77 nm deviation on 238 phosphorus atoms that had not
   moved: about **+138,000,000 kJ/mol** on a system at rest. Written as
   `periodicdistance(0, 0, z, 0, 0, z0)^2` the same system reads
   **-2,213,385** with both restraints at exactly zero. This was the real
   blocker and it took four wrong hypotheses to reach.
2. gemmi's placeholder `CRYST1` of 1 1 1 made `addMembrane` size the box from
   the lipid patch rather than from the 9.4 nm receptor (D014).
3. The barostat was on during what section 5.2 defines as NVT heating.
4. The built system's contents exceed its own periodic box, so `addSolvent`
   placed 246 waters where overhanging lipid and patch-water atoms already sat.
   Repaired by deleting them, preferring water over ions so neutrality holds.

**One of those four wrong hypotheses is worth keeping.** The +136,000,000 was
attributed to single precision for several hours, on the strength of a CPU run
reading -2,060,000 for what was described as the same system. The two runs did
not differ in precision: one had no restraints and the other had both. Nothing
measured anywhere in this build ever showed the two platforms disagreeing. A
difference between two runs is only evidence about the variable that actually
differed, and an earlier version of this file carried the wrong version of it.

**The membrane ns/day is now measured, not projected.** The Tier 1 run gives
**59.8 ns/day** with the barostat on, against the 86 the heating stages had
suggested: the `MonteCarloMembraneBarostat` costs about 30 per cent. That was
the one number the compute budget turned on, and it closes G1 rather than
opening it. Including the equilibration that section 5.2 puts in front of every
row, the four rows come to **59.4 GPU-hours** against a 96-hour ceiling, with
36.6 spare. Both membrane systems came in near 194 k atoms rather than the
spec's 80 to 110 k (D012), which is why the projection was worth replacing with
a measurement.
