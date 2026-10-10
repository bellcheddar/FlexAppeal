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
