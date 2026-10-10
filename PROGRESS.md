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
| **5 Production** | not started | Tier 2, fifteen runs at 100 ns (1,500 ns) |
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

1. **Compute budget (Gate G1), and it is very likely to open.** Tier 2 is
   **fifteen** production runs of 100 ns, not twelve: section 5.3's first row is
   three replicas *per protonation arm*, so ABL1 active-like is six runs on its
   own. Six ABL1 active-like, three ABL1 inactive, three β2AR active, three β2AR
   inactive: **1,500 ns in total**.

   Projected from this morning's measured 298 ns/day at 25,463 atoms on this
   machine and platform, scaling inversely with atom count:

   | Row | Atoms | ns/day | Hours per run | Runs | Hours |
   |---|---|---|---|---|---|
   | ABL1, water (three rows) | ~50 k | ~152 | 15.8 | 9 | 142 |
   | β2AR, POPC (two rows) | ~95 k | ~80 | 30.0 | 6 | 180 |
   | | | | | **15** | **322** |

   That is **13.4 days of continuous GPU against a 96 hour ceiling, 3.4 times
   over**. Scope levers the spec allows before the gate: two replicas rather than
   three gives 215 hours, halving Tier 2 to 50 ns gives 161, both together give
   107, still just over.

   This is an extrapolation and section 2.2 is explicit that budgets come from
   measurement, so the gate is not opened on it. Tier 1 in phase 3 measures both
   systems for real and the number is settled then. Two things could move it
   either way: a membrane system under `MonteCarloMembraneBarostat` is usually
   slower than atom count alone predicts, and `concurrent_md = 1` means none of
   this overlaps.
2. **β2AR start models.** The spec's primary route is a GPCRdb *refined* model,
   and it says plainly that not every entry has one (2RH1 and 3SN6 both report
   none). The fallback ladder runs through state-specific homology models to
   building one locally. Which rung is reached is discovered, not assumed, and is
   flagged in the app for the reader.
