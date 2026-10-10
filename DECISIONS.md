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

---

## D006 - One replica per row, not three

**2026-10-10, phase 0.** Ruled by Marc after the first compute projection. Tier 2
goes from fifteen runs to **five**, and from 1,500 ns to 500: ABL1 active-like
with the Asp charged, ABL1 active-like with it protonated, ABL1 inactive, β2AR
active, β2AR inactive, one trajectory each.

Projected cost falls from 322 GPU-hours to **107**, 4.5 days rather than 13.4.
That is still 1.12 times the spec's 96 hour default, which is close enough that
the phase 3 measurement decides it: at 100 ns/day for the membrane rows it comes
in at 95 hours and under the ceiling, at 60 ns/day it is 127 and over. Section
0.2 opens G1 only when the budget exceeds the ceiling *after* scope reduction,
and this is that reduction, so the gate stays shut unless the measured number
says otherwise.

**What it costs, stated plainly because it is a real loss and not a tidy-up.**

- **Positive control P5 cannot be evaluated as written.** It reads "≥ 2 of 3
  replicas keep whole-protein `s` < 0 over the final 20 ns". With one trajectory
  there is no majority to take. It becomes a single observation, reported as
  such: either that trajectory stayed inactive or it did not.
- **No estimate of run-to-run variability.** A conformational change in MD is a
  stochastic event, and three replicas exist so that "the active-start β2AR
  relaxed" can be distinguished from "one trajectory did something". With n = 1
  every finding in `FINDINGS.md` and every sentence on the page is one
  trajectory's behaviour and must say so. No claim of reproducibility can be
  made from this data.
- **A single run that goes wrong costs the whole row.** With three replicas, one
  trajectory that unfolds, drifts out of the box or hits a bad start is survivable
  and the other two carry the row. Here it is the row. The checkpoint-and-resume
  machinery (section 0.1 rule 4) matters more, not less.

**What it does not cost.** Every positive control except P5 is measured on the
reference structures or on a single frame, so P1, P2, P3, P4 and P6 are
unaffected. The annotation module, the state score, the motif metrics and the app
are all per-frame work and do not care how many trajectories exist.

*Reverse:* set replicas back to 3 in the run matrix and re-run; the per-run
pipeline is identical and the analysis keys on `run_id`, so added replicas appear
as extra rows rather than needing anything rebuilt.

---

## D007 - One ligand-bound run and one apo run per system, in opposite states

**2026-10-10, phase 1.** Ruled by Marc: "1 ligand should be included for each
system, the other apo partner should be in the opposite state so we can see the
differences." This replaces section 4.2 and 4.3's "ligands removed for the core
runs" and promotes half of the optional phase 5.7 into the core.

| # | System | Start | Ligand | State | Expectation |
|---|---|---|---|---|---|
| 1 | ABL1 | 2HYY | **imatinib** (`STI`) | DFG-out | holds: imatinib is the canonical type II binder and wedges DFG out |
| 2 | ABL1 | 2GQG | apo (dasatinib stripped) | DFG-in | free to flip; the protonation arm acts here |
| 3 | ABL1 | 2GQG | apo, **DFG-Asp protonated** | DFG-in | the spec's headline kinase feature, kept (below) |
| 4 | β2AR | 2RH1 | **(S)-carazolol** | inactive | holds: an inverse agonist stabilising the inactive state |
| 5 | β2AR | GPCRdb refined | apo | active | drifts: Dror 2011 found active relaxes without G protein |

**Which state gets the ligand, and why that way round.** For each system the
ligand goes on the state that ligand is known to stabilise, and the apo run takes
the opposite state. That gives one trajectory expected to hold and one expected to
move, which is the comparison asked for. The reverse pairing was available and is
weaker: dasatinib is a type I binder whose hold on DFG-in is less decisive than
imatinib's on DFG-out, and an agonist-bound active β2AR still relaxes without a G
protein, so neither would have been a clean control.

**The protonation arm is kept**, as run 3. The instruction was about ligands and
apo partners, not about protonation, and section 5.2 calls the DFG-Asp arm "the
headline FlexAppeal feature for kinases": the same protein, one proton apart,
different motif behaviour. It is a variant of run 2 rather than a new system, so
it costs one more ABL1 row at about 16 hours. Flagged to Marc rather than assumed;
dropping it returns the matrix to four runs.

**Higher-resolution seeds exist and were not taken.** KLIFS lists 5hu9 at 1.53 Å
and 3qrj at 1.82 Å for DFG-out against 2HYY's 2.4 Å. Section 4.2 names 2GQG and
2HYY as preferred and says to accept them if KLIFS labels them as expected, which
it does: 2GQG is DFG-in/αC-in, 2HYY is DFG-out/αC-in, both at 2.4 Å. Following the
spec's named seeds also keeps continuity with BoltzMaker, which already uses 2GQG.

**What this gains.** Two of the three entries in the Runs library now carry a
ligand, so they get the binding-site panels and a grade alongside the motif and
state analysis. The library becomes three graded runs rather than one graded and
two ungraded.

**The complication, logged before it bites.** A GPCRdb refined model has its
fusions removed and carries no ligand, so run 4 cannot simply use one. Carazolol
has to be transplanted: superpose 2RH1 onto the refined inactive model on the 7TM
core, copy the ligand across, then check for clashes before parameterising. gemmi
does the superposition and `sse_comparison` already does apo-onto-holo alignment
on non-binding-site residues. If the transplant will not sit cleanly, the fallback
is to build the inactive start from 2RH1 itself with the T4 lysozyme fusion
deleted, which is rung 3 of section 4.3's ladder and is flagged in the app as
"start model built locally".

*Reverse:* set `ligand_runs = false` and restore the apo pairs from sections 4.2
and 4.3.

---

## D008 - Two pairs, and imatinib is the ABL1 ligand

**2026-10-10, phase 1.** Marc: "Let's just do 2 pairs for each, choose which is
mostly likely to hold for abl." This amends **D007** by dropping the protonation
arm, leaving one pair per system and four runs in total.

| # | System | Start | Ligand | State | Expectation |
|---|---|---|---|---|---|
| 1 | ABL1 | 2HYY | **imatinib** (`STI`) | DFG-out | holds |
| 2 | ABL1 | 2GQG | apo (dasatinib stripped) | DFG-in | free to flip |
| 3 | β2AR | 2RH1 | **(S)-carazolol** | inactive | holds |
| 4 | β2AR | GPCRdb refined | apo | active | drifts |

**Why imatinib rather than dasatinib.** Imatinib is a type II inhibitor: it binds
the allosteric pocket that only exists once the DFG motif has flipped out, so the
drug physically occupies the space the motif would have to move back through. The
DFG-out state is held by steric occlusion rather than by affinity alone.
Dasatinib is type I, binding the ATP site with DFG in, and DFG-in is the state the
kinase reaches more readily in any case, so a dasatinib run would be a weaker test
of whether a ligand holds a state: it would be hard to tell a ligand holding the
conformation from the conformation simply being stable.

That matters because the pair is a comparison. Run 1 should hold and run 2 should
be free to move, and the clearer the hold, the more a change in run 2 means
something.

**The protonation arm is dropped**, which was the point of raising it in D007.
Section 5.2 calls it the headline kinase feature and it is not in this build. It
is a variant of run 2 rather than new machinery, so it can be added later as a
fifth row for about 16 hours of compute with no change to the pipeline.

**This brings the projection inside the ceiling.** Four runs, 400 ns, about **92
GPU-hours or 3.8 days against the 96 hour default**, so G1 does not open on the
estimate. The headroom is four hours, which is thin: at 70 ns/day rather than 80
for the membrane rows the total reaches 100 hours and the gate opens. Phase 3
measures it.

*Reverse:* add the protonation arm back as run 5, or restore D007's five-run
matrix.

---

## D009 - Tier 2 is 50 ns a run, not 100

**2026-10-10, phase 1.** Ruled by Marc. Four runs of 50 ns, 200 ns in total,
about **46 GPU-hours or 1.9 days**.

**The reason is not the compute saving, although there is one.** Neither length
reaches the transitions these two systems are known for. The spec says so itself
about the receptor: the β2AR active-to-inactive relaxation Dror 2011 measured
happened on the microsecond scale, and section 9 already records its absence at
100 ns as an expected limitation rather than a failure. The ABL1 DFG flip is the
same: Shan 2009, the paper the protonation arm came from, puts it on microseconds
too. 100 ns is twenty times short of the headline event and 50 ns is forty times
short. Both sit in the same regime, where what is observed is motif fluctuation
around the starting state rather than the transition, and doubling the run does
not cross that gap.

**What decided it was the single replica.** D006 and D008 left one trajectory per
row, so a run that goes wrong costs the whole row and there is no second replica
to carry it. At 100 ns the budget is 92 hours against a 96 hour ceiling, and one
β2AR redo puts the build at 122 hours and opens G1. At 50 ns the budget is 46
hours and every run could be redone once and still come in at 61. Headroom goes
from four hours to fifty.

**What it costs, stated rather than glossed.** P5's window is "the final 20 ns",
which is a fifth of a 100 ns run and two fifths of a 50 ns one, so it is less
independent of the start. And if a flip did happen spontaneously, a 100 ns run is
twice as likely to catch it, which is a small difference in an already small
probability. 500 frames at 100 ps is the same frame count the committed lysozyme
example packs, so no panel is short of data.

*Reverse:* one number in the run matrix. Re-running longer needs no change to the
pipeline, and the analysis keys on `run_id`, so a longer re-run appears as another
row rather than replacing anything.
