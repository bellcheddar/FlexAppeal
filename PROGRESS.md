# Progress: reference-state examples

Phase checklist for `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`. Rewritten after every
phase and at least every 30 minutes during long compute (section 0.1 rule 5).

**Last written:** 2026-10-10, end of phase 0.
**Gates open:** none.

| Phase | State | Notes |
|---|---|---|
| **0 Setup** | **done** | Hardware probed, BoltzMaker vendored and pinned, dependencies resolved, both suites green |
| **1 References** | next | KLIFS and GPCRdb fetchers, reference ensembles, per-reference metrics |
| **2 Start models** | not started | ABL1 seeds (2GQG, 2HYY); β2AR refined models from GPCRdb |
| **3 System build + smoke** | not started | Solvated ABL1, POPC-embedded β2AR, Tier 1 (1 ns per row) |
| **4 Analysis module** | not started | `flexappeal/refstates/` over the vendored `sse_comparison` |
| **5 Production** | not started | Tier 2, twelve runs at 100 ns |
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

1. **Compute budget (Gate G1).** Tier 2 is twelve production runs of 100 ns:
   ABL1 in water across two protonation arms and two start states, three replicas
   each, plus β2AR in a POPC bilayer from two start states, three replicas each.
   The membrane system is expected at 80 to 110 k atoms against lysozyme's 25 k,
   so its ns/day will be a fraction of this morning's 298. Until phase 3
   measures it, the total is unknown and the default ceiling is 96 hours. If the
   measured projection exceeds that after scope reduction, G1 opens.
2. **β2AR start models.** The spec's primary route is a GPCRdb *refined* model,
   and it says plainly that not every entry has one (2RH1 and 3SN6 both report
   none). The fallback ladder runs through state-specific homology models to
   building one locally. Which rung is reached is discovered, not assumed, and is
   flagged in the app for the reader.
