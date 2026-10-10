# Decisions

Choices made during the reference-state build (`FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`),
each with the alternatives, the reason, and how to reverse it. Section 0.1 rule 1:
where the spec is ambiguous, pick the most defensible option and record it here
rather than stopping.

---

## D001 - The build happens inside the existing FlexAppeal repository

**2026-10-10, phase 0.** The spec names
`/Users/dellboy/Documents/Vibe_Coding/FlexAppeal/` as the project root and says to
use the existing FlexAppeal code if it lives elsewhere. It does live here, and it
is not a skeleton: a working three-tab Flask app with 759 passing tests, a
deployed droplet, an `.fxa` reader and writer, a measurement module and a panel
builder.

So this is an addition to a live application rather than a new project, which
changes three things the spec allows for:

- **The design system exists** (section 6.1), so the motif tokens join
  `brand.css` rather than bringing the marcdeller.com defaults in fresh.
- **There is already a Runs tab**, added earlier today, which is exactly the
  "examples tab" of section 6.2 item 1. The two new systems become entries in
  that library rather than a second listing beside it.
- **`.claude/settings.json` and `CLAUDE.md`** in appendices 12.1 and 12.2 are
  written for a fresh checkout. Applied as given they would narrow what is
  already permitted here and describe only this build. Recorded as **D002**.

*Reverse:* nothing to reverse; the alternative was a separate repository, which
would have meant a second copy of the `.fxa` reader, the panel builder and the
deployment.

---

## D002 - The spec's CLAUDE.md is added as a section, not as the whole file

**2026-10-10, phase 0.** Appendix 12.2 gives a `CLAUDE.md` for a project whose
only job is this build. Written over the top it would drop whatever the repository
already tells a reader and would state, wrongly, that this spec is authoritative
for everything here. The spec's rules are added as their own section, scoped to
the reference-state work, and the house rules that already applied continue to.

The same applies to `.claude/settings.json`: the deny list is kept (it forbids
`git push`, which is Gate G4 and matters), and the allow list is not narrowed to
the appendix's set, which omits tools this repository already uses.

*Reverse:* delete the section and restore the appendix verbatim.

---

## D003 - BoltzMaker is vendored from GitHub at `da939b1`, not symlinked locally

**2026-10-10, phase 0.** `sse_comparison` also sits at
`~/Documents/Vibe_Coding/BoltzMaker`, which would have saved the clone. The spec
says to clone at the current `main` commit and record the SHA, and a working copy
on this machine is not a fixed thing: it can be edited, rebased or left dirty, and
then this build's motif annotation would quietly depend on an unpublished state.

Vendored at **`da939b1be54a27cb8eac7c937fd28e8d1e7d008d`** into
`vendor/boltzmaker/`, read-only, gitignored. Wrapped, never edited (section 4.1);
anything that needs a change inside it becomes a FlexAppeal subclass or adapter
and is logged here as a candidate upstream change.

*Reverse:* `rm -rf vendor/boltzmaker` and re-clone at a different SHA.

---

## D004 - Five dependencies added, four of the spec's left out

**2026-10-10, phase 0.** Section 3 lists the packages. Most were already here
(openmm, pdbfixer, mdtraj, gemmi, biopython, pandas, scipy, numpy, rdkit, parmed).
Added: `mdanalysis` (membrane analysis and helix axes), `requests` and `tenacity`
(the KLIFS and GPCRdb clients, where the spec's fallback of bare `urllib` would
mean hand-rolling retries and caching), `pyarrow` and `duckdb` (the per-residue
DSSP store, which is too wide for SQLite and too large for JSON).

Two more were added for a reason that only appeared on trying it, **`openpyxl`
and `plotly`**: the vendored `sse_comparison` modules import BoltzMaker's main
module for `_align_positions`, and that module imports a spreadsheet writer and a
plotting library at module scope. Neither is used by FlexAppeal. The alternative
was to extract `_align_positions`, which is the numbering-offset reconciliation
the spec explicitly says not to reimplement.

BoltzMaker's own pins are deliberately not adopted. It holds gemmi at 0.6.5,
scipy at 1.13 and biopython at 1.84 to match boltz's pip metadata, and this
environment runs 0.7.5, 1.17 and 1.87, which FlexAppeal's own code and its 759
tests are built against. Taking the pins would have downgraded three libraries
under a working application to satisfy a vendored helper. Only the two
import-time packages were taken, and the full suite was re-run after each solve.

Left out, each on the spec's own fallback ladder:
- **`mkdssp`** has the documented fallback `mdtraj.compute_dssp(simplified=False)`,
  which is already installed and is what the per-frame path uses anyway. Adding a
  second DSSP implementation would invite the two to disagree on the reference
  structures, which is where section 9's controls are measured.
- **`packmol-memgen`** is the *fallback* membrane builder; OpenMM's
  `Modeller.addMembrane` is primary (section 5.1). It arrives with ambertools,
  which is already a dependency, so it can be reached if the primary route fails.
- **`parasail`** is a fallback for biopython, which is present.

The lock was backed up before the solve and the full suite re-run after it, because
this environment also serves the deployed application.

*Reverse:* remove the five lines from `pixi.toml` and restore the backed-up lock.


---

## D005 - The vendored clone is given back its gitignored test data

**2026-10-10, phase 0.** Four of BoltzMaker's 50 `sse_comparison` tests failed in
the fresh clone with `FileNotFoundError` on
`examples/adrb2_gs_panel/boltz_cif/*.cif` and the EGFR equivalent. Those files are
in BoltzMaker's own `.gitignore` (`boltz_cif/`), so no clone of that repository can
carry them and the failure is not about this build.

The three CIFs were copied from the working copy at
`~/Documents/Vibe_Coding/BoltzMaker` into the vendored tree: 1.2 MB of structures,
which the spec names as the Phase 1 smoke-test data anyway (section 4.1, "its
`adrb2_gs_panel` example and test fixtures are the first smoke-test data"). All 50
tests then pass.

This is data, not code, so it does not breach the rule against editing the
vendored copy. It does mean the vendored tree is not reproducible from the SHA
alone, which is recorded here rather than left to be discovered.

*Reverse:* delete `vendor/boltzmaker/examples/*/boltz_cif/`; the four tests then
fail again with a clear file-not-found.
