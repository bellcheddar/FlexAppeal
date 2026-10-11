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

---

## D010 - Two deviations forced by the hardware and by the data

**2026-10-10, phase 1 and 3.** Both are departures from the spec, both measured
rather than assumed, and both are wrapped rather than edited into the vendored
copy.

### Single precision on the GPU, because Apple's OpenCL has no other

Section 2.2 expects `OpenCL` at `mixed` precision on Apple Silicon and section
5.2 says `never single for production`. Apple's OpenCL implementation is single
precision only, and asking for mixed fails context creation with "No compatible
OpenCL platform is available", which names the platform rather than the property
that caused it. FlexAppeal has met this before: the measured cost of taking the
message at face value and falling back to the CPU was **29.5 ns/day against 206**.

Each candidate is now tried for real, in order, and the first that builds a
context wins: CUDA mixed, OpenCL mixed, OpenCL single, CPU. On this machine that
lands on OpenCL single, which is what the committed 100 ns benzene example also
ran on. The energy-conservation check section 5.2 asks for is run on the Tier 1
smoke trajectory rather than taken on trust.

*Reverse:* delete the `single` entry from the attempt list and the build falls to
the CPU, at roughly a seventh of the speed.

### The DFG call is read from the nearer anchor, not from an absolute threshold

Positive control P1 failed at **11/17, 64.7%**, against a 90% bar. Every one of
the six mismatches had the same shape: KLIFS `out`, ours `other`, never a
confusion between the two states.

The vendored classifier measures the DFG-Phe ring against the αC-Glu+4 Cα and the
catalytic Lys Cα, then asks whether each distance falls the right side of 11.0 Å.
Measured across ABL1's references:

| | ring to αC | ring to Lys | difference |
|---|---|---|---|
| DFG-in (5) | 5.6 to 6.8 Å | 12.9 to 14.0 Å | **-6.1 to -7.9** |
| DFG-out (12) | 13.1 to 16.3 Å | 9.1 to 11.3 Å | **+2.4 to +5.8** |

The two classes separate by **8.5 Å** on the difference, with nothing in between.
The absolute threshold fails because six DFG-out structures sit at 11.1 to 11.3 Å
from the Lys: a tenth of an ångström past a line drawn through the middle of a
dense cluster. 2HYY, this build's own DFG-out start structure, is one of them.

So the state is read from which anchor the ring is nearer, which places no
threshold where the data is dense. **P1 goes to 17/17.** A margin of 1.0 Å around
equidistant still returns "other", placed where the data is empty rather than
where it is dense: that band matters for trajectory frames rather than
references, because a DFG mid-flip genuinely is neither state and a bare sign test
would snap every frame to one side and hide the transition.

The distances are stored beside every call, so the label can be checked against
the geometry and the absolute-threshold answer stays recoverable without re-running
anything.

**Candidate upstream change for BoltzMaker.** The same 11.0 Å threshold is in its
own `compare-sse` and will mislabel the same six structures. Wrapped here in
`flexappeal/refstates/kinase.py`, never edited there (spec 4.1).

*Reverse:* call `sse_comparison.metrics.classify_dfg_state` directly; P1 returns
to 64.7% and fails.

---

## D011 - The β2AR starts, and the transplant D007 planned for is not needed

**2026-10-10, phase 2.** The fallback ladder in section 4.3 was walked rather
than assumed, and it stopped on the top rung for both states.

| Run | Entry | Resolution | Refined model | Ligand |
|---|---|---|---|---|
| apo, active | **4LDE** | 2.79 Å | yes | agonist stripped |
| ligand, inactive | **5D5A** | 2.48 Å | yes | **(S)-carazolol, in place** |

**2RH1 has no refined model**, which the spec predicted and this build confirmed
by asking. Probing the best entries of each state by resolution found refined
models for 4LDE and 6MXT among the active, and 6PS2, 5D5A and 6PS3 among the
inactive. The URL pattern is `/structure/refined/<PDB>_full/download_pdb`, it
returns a zip rather than a PDB, and a `_noloops` variant sits beside it that
drops ICL3 entirely.

**D007's complication does not exist.** It recorded that a refined model carries
no ligand and that carazolol would have to be transplanted from 2RH1 by
superposition, with a clash check and a documented fallback. Refined models do
carry their ligands, and 5D5A's is the exact inverse agonist that was wanted.
So there is no transplant, no clash risk and no fallback rung: the ligand is
where the crystallographers put it. 6PS2 at 2.40 Å is marginally better resolved
but carries an antagonist rather than an inverse agonist, and the inverse agonist
is the stronger stabiliser of the state the run is meant to hold.

**ICL3 is fully modelled, 24 of 24 residues, in both.** Section 4.3 says over 20
means keep them, flag the loop as "modelled, flexible" and exclude it from every
RMSD. Recorded in each start's sidecar as `icl3_excluded_from_rmsd`.

**Orientation.** Neither refined model arrives in the membrane frame: they were
50 Å apart in z and tilted differently from each other, and `addMembrane` needs
the bilayer normal along z with the centre at z = 0. OPM has no entry for 4LDE or
5D5A, so the frame is borrowed from the receptor rather than the entry, which
section 5.1 allows: the bilayer frame belongs to the fold. Superposed on 225 TM
Cα onto OPM's 2RH1, whose dummy atoms put the headgroup planes at z = ±15.9 Å.
Loops are excluded from the fit because ICL3 is fully modelled here and is the
most mobile part of the receptor, so fitting on it would tilt the helices to
satisfy a loop.

Both pass section 5.1's validation: aromatic belt |z| median 12.0 and 12.6 Å
against a 10 to 17 window, intracellular basics at -37.5 and -27.9 Å against a
-13 ceiling.

The fit RMSDs are themselves a check worth keeping: the inactive model lands on
the inactive OPM reference at **0.20 Å** and the active at **2.51 Å**. That 2.3 Å
gap is the activation difference this build exists to measure, showing up before
any simulation has run.

---

## D012: the membrane systems are 190 k atoms, not the 80 to 110 k expected

Section 5.1 expects "roughly 80 to 110 k atoms" and prescribes a response to more
than 150 k: reduce `minimumPadding` to 1.2 nm and log it. That is what happens,
and both systems still land well above the window:

| system | atoms | lipids | waters | box, nm |
|---|---|---|---|---|
| β2AR apo (active) | 194,626 | 476 | 41,788 | 12.74 × 13.37 × 11.77 |
| β2AR carazolol (inactive) | 188,352 | 476 | 39,686 | 12.74 × 13.37 × 11.39 |

Two causes, and neither is fixable by padding.

**Height.** ICL3 is fully modelled in the GPCRdb refined models and hangs down
to z = -61 Å, 35 Å below the next lowest element. Everything else spans -42 to
+32 Å. The loop alone sets the box height, so the box is 11.8 nm rather than the
9.8 nm the receptor would otherwise need.

**Width.** `addMembrane` tiles whole copies of a pre-equilibrated patch, so the
lateral box is an integer number of patches. The receptor is 5.9 nm across and
one patch is about 6.4, so anything but a negligible padding needs two patches
each way and the lateral box is 12.7 × 13.4 nm at any padding worth using.

Kept as built rather than trimmed. Section 4.3 says to keep what the refined
model contains and flag the loop, which is what `icl3_excluded_from_rmsd` already
does, and a receptor whose ICL3 is deleted to save water is a different construct
from the one the annotation describes. The cost is compute, which is measured
rather than extrapolated: see the Tier 1 figures in BUILD_LOG.md.

## D013: the MD queue is strictly sequential

There is one GPU. The β2AR smoke test and the ABL1 production run shared it for
fifteen minutes, and the smoke test got 8% of a core and never reached its first
frame while ABL1 slowed from 213 to 205 ns/day. Section 2.3's concurrency is
between MD and *analysis*, not between two MD runs, and the system building that
runs alongside is CPU work for the same reason.

`tools/queue_md.py` therefore runs one row at a time, skips a row that already
has `run.json`, and can be killed and restarted because `run_md.py` resumes from
its own checkpoint. Four rows at 50 ns is most of a day and the machine sleeps,
so resumability is the difference between losing a nanosecond and losing a row.

Tier 1 is run before the long runs rather than after, which cost less than one
nanosecond of ABL1: the membrane protocol in section 5.2 is substantially
different from the water one, with a membrane barostat and a six-step restraint
release, and none of it had ever been stepped. Finding a fault in it after five
hours of ABL1 would have been five hours spent to learn the same thing.

## D014: `addMembrane` believed a placeholder unit cell, and produced NaN

Both membrane builds died at `addMembrane`'s inflation relaxation with "Particle
coordinate is NaN" on its very first step. Worth recording because the first two
hypotheses were wrong and the third was invisible from the error.

Not a clash: the assembly going in minimised cleanly from 150,000 to 8,286
kJ/mol, and the sub-1.6 Å contacts a clash check reported were peptide bonds.

Not precision: Apple's OpenCL is single-precision only and `addMembrane` creates
its own context without specifying a platform, which was a plausible cause and a
real property of this machine. It fails identically on the CPU platform at double
precision.

The cause is that gemmi writes a placeholder `CRYST1` of 1 1 1 for a `P 1`
structure with no cell. `addMembrane` takes the box height as the larger of the
patch height and the existing cell plus padding, so with a cell present at all it
took 7 nm from the patch rather than sizing from the receptor, which is 9.4 nm
tall. The receptor wrapped through the periodic boundary onto itself. Clearing the
cell makes `addMembrane` size z from the protein, which is what it does when no
cell is present.

The relaxation platform is still pinned to the CPU, now for the reason in D013
rather than for precision: it leaves the GPU free for production.

## D015: the R-spine contact threshold is 7.7 A, not 6.0

Section 5.4 calls the regulatory spine `assembled` if every consecutive pair of
side-chain centroids is within 6.0 A. On ABL1's reference set the pair that
decides the call, DFG-Phe to HRD-His, separates the states completely: active
5.33 to 5.89 A, inactive 9.53 to 11.08 A, nothing in a 3.6 A gap between them.

So 6.0 A is right for crystal structures and sits a tenth of an angstrom above
the top of the active range. That left no room for a side chain that is
breathing, and trajectory frames of the apo DFG-in system read 6.03 to 6.69 A,
clearly on the active side of the gap, and were all called `broken`. The score
then carried a standing `conflict` on the DFG motif for all 43 frames: geometry
active at s = +0.755, rule inactive.

Moved to 7.7 A, in the middle of the empty region rather than at the edge of the
active cluster. This is the same correction D010 made to the DFG criterion and
for the same reason: a crystal structure is one minimised conformation and an MD
frame at 300 K is not, so a threshold calibrated on crystals needs to sit where
crystals are absent.

Checked both ways. The references still classify perfectly at the new threshold,
5 of 5 active `assembled` and 14 of 14 inactive `broken`, and the trajectory's
standing conflict disappears. The published 6.0 A is kept as
`R_SPINE_CONTACT_CRYSTAL_A` for structure-only comparisons.

**Re-derived, with two numbers above corrected.** This paragraph was written
from a measurement that nothing in the repository could reproduce, which is the
fault 4g describes for the manifests. `tools/build_references.py` now recomputes
it on every kinase verify. The separation holds and the conclusion is unchanged,
but the figures need amending:

| quantity | stated above | re-derived |
|---|---|---|
| active DFG-Phe to HRD-His | 5.33 to 5.89 A | 5.33 to **5.96** A |
| deciding pair, active | not stated | 5.33 to 5.96 A |
| deciding pair, inactive | 9.53 to 11.08 A | **10.35 to 12.92** A |
| empty region | 3.6 A | **4.39** A, from 5.96 to 10.35 |
| 6.0 A clears the active cluster by | "a tenth of an angstrom" | **0.04 A** |

The ranges above were the DFG-Phe to HRD-His pair alone. The call is made on
whichever pair is *widest*, since one pair over the threshold makes the state
`broken`, so that is what the control reports and what the empty region is now
measured on. Correcting the active bound to 5.96 narrows the published
threshold's margin from a tenth of an angstrom to four hundredths, which
strengthens the case for moving it rather than weakening it.

**And one reference does turn on the choice.** 2G1T, the Src-like inactive
start, has a deciding pair of 6.14 A: `assembled` at 7.7 and `broken` at the
published 6.0. D015 counted 19 references and left the Src-like ensemble out of
that count, so this was not visible before. It is the clearest single
illustration of the point, because 2G1T is DFG-in with the alphaC helix out and
its spine stack is intact, so `assembled` is the right answer and 6.0 gives the
wrong one.

## D016: the kinase manifest holds one row per deposition, not per chain

KLIFS's `structures_list` returns a row per chain and per alternate
conformation. For ABL1 that is 61 rows over 20 depositions: 2G1T appears four
times, 2HZ4 three times. The manifest keeps the first row for each deposition
and `tools/build_references.py --verify` reproduces the committed file exactly
on that rule, including the six entries whose chain is not A.

The alternative would count every chain, giving 37 usable references instead of
17. It was tried, and it passes: P1 37/37 and P2 36/37. It is still the wrong
denominator. Chains of one asymmetric unit are copies in a single crystal, so
they are not independent tests of whether our label agrees with KLIFS's, and
counting them gives a four-chain deposition four times the weight of a
single-chain one while quoting an n that overstates what was measured.

Kept as `--all-chains`, which reports and refuses to write the manifest, because
it answers a genuinely different question: whether a label is stable within one
crystal. It found that 2HZ4 chain A reads the catalytic Lys to alphaC Glu
distance at 4.20 A against a 4.0 A cutoff and so calls alphaC `out`, where
chains B and C read 2.71 and 2.74 A and call it `in`. A 0.2 A threshold miss in
one copy, and P2 would still pass at 16/17 had that copy been listed first.

This is deliberately not applied to the receptor manifest, which is keyed on
`(pdb, chain)` because GPCRdb lists a preferred chain per entry and the 116 rows
are already one per structure.

## D017: one barostat helper, and checkpoints that say which phase they are

A resumed run had no barostat. The force was added inside the equilibration
branches, and a resume skips those, so production continued at constant volume.
OpenMM restores context parameters whether or not a force consumes them, so the
resumed context reported `MonteCarloPressure` from a System with no barostat in
it and the load succeeded without complaint.

**Both call sites now go through `barostat(membrane, temperature_k)`.** They had
each constructed their own, which is how the resume path came to have none, and
it is also where the membrane choice lives: `MonteCarloMembraneBarostat` at zero
surface tension, XY isotropic, Z free, so the bilayer sets its own area. A plain
isotropic barostat would scale z with xy and squeeze it.

**Checkpoints now record their phase in `checkpoint_phase.json`.** The membrane
ladder saves after every rung, and those checkpoints carry a non-zero clock and
a non-zero `k_restraint`, so they are not production checkpoints even though
they are the same file. An absent marker means production, which keeps every
checkpoint written before this existed resumable, and an equilibration marker
makes the run repeat its equilibration rather than resume into production with
the protein still held.

Repeating equilibration costs 10 ns on a membrane row, about four hours. That is
the right trade against a 50 ns row that is restrained, short and labelled
complete, because the second failure is one nobody would catch by looking.

The alternative was resuming the ladder at the recorded rung. Rejected for now:
it needs the heating stages and the barostat to be reconstructed in exactly the
state the rung expects, and a resume path that is itself hard to verify is not
an improvement on one that starts again from a known point.

## D018: the ABL1 rows are dephosphorylated, and that is now written down

`abl1_apo_start.pdb` carries PTR 393, the activation-loop phosphotyrosine.
PDBFixer's `replaceNonstandardResidues` rewrites it to TYR, so both ABL1 rows
simulate the dephosphorylated kinase domain started from a phosphorylated
crystal structure. Deposited residue 393 reads TYR in the built system and the
system contains no PTR at all.

**Kept, not re-run.** Amber14 ships no phosphotyrosine, so phosphorylating it
means bringing in a separate parameter set, and the four rows are already most
of a 96-hour budget with abl1_apo 4 ns in. Simulating the dephosphorylated
domain is also the common choice in the ABL1 literature, so this is a defensible
state rather than a broken one.

**What was wrong was the silence.** Nothing recorded the substitution, in a
build whose whole subject is which conformational state a kinase is in and where
Tyr393 phosphorylation is one of the main things that decides it. Both builders
now log every replacement, and FINDINGS 4j states what it means for the results:
the ABL1 rows measure a dephosphorylated kinase, P4 is unaffected because it
reads crystal geometry, and a drift towards inactive in the apo DFG-in row now
has a physical explanation that has to be offered before any other.

The spec does not mention phosphorylation, so this is not a departure from it.
It is a choice the spec left open, and the build had made it without noticing.

## D019: `run_md.py` imports the ligand recipe rather than keeping its own

A water row with a bound ligand needs a SMIRNOFF template built from the same
molecule object as the topology. Handing the generator the CCD SDF fails with
"No template found for residue 262 (STI)" even though both hold 68 atoms,
because the generator matches a chemical graph and a PDB-derived ligand has
bonds inferred from distance.

`tools/build_membrane.py` had already solved this for carazolol. `run_md.py` now
imports `ligand_molecule` from it instead of carrying a second implementation,
and follows the same sequence: replace nonstandard residues, take the ligand out
before the hydrogens go on, add hydrogens to the protein, then merge the ligand
back from the object the generator was given.

`--ligand` therefore takes a CCD residue code, not an SDF path. The old
signature could not express what is needed, since the chemistry has to come from
the dictionary and the coordinates from the deposition.

The queue now carries a `ligand` field per row, and a test reads each start
model and asserts that every het residue is either that row's declared ligand or
one PDBFixer rewrites. Without it `abl1_imatinib` would have crashed at build
time twenty hours into the queue.

## D020: the reuse path caps its reference ensembles, the production path does not

`flexappeal refstates annotate` is section 6.5's G-REUSE hook, and on a
well-studied kinase it tried to download every structure the catalogue offered:
CDK2 returns 337 that KLIFS calls active, and the first corrected run was still
fetching after fourteen minutes.

**Capped at 40 per state for `annotate`, best first.** By KLIFS's own
`quality_score` descending for a kinase and by GPCRdb's `resolution` ascending
for a receptor, since neither catalogue offers the other's column. An entry with
no value for either sorts last rather than first, which needed its own fallback
on the ascending side: a missing value read as minus infinity would have been
the best resolution there is.

The justification is what the score is. A motif's state comes from its distance
to the *nearest* active and nearest inactive reference, so beyond a few dozen
structures an extra one only wins if it is nearer than everything already there,
and the cost is linear. Forty is above both worked examples' largest state (45
inactive for ABL1 is the one exception, and the production path is uncapped so
it keeps all 45).

**The production path passes no limit.** That is deliberate and it is the whole
reason the cap is an option rather than a constant: every committed number,
including P1, P2, P3, P4 and both manifests, is computed over the full
ensembles. Checked after the change: P4 still reads +0.305, -0.368, +0.593 and
-0.934, and ABL1 still resolves 12 active and 45 inactive.

`--max-references 0` turns the cap off for a user who wants the full set.

## D021: a thin reference ensemble is reported, not silently used

CDK2's reuse run kept 40 active references and the core alignment dropped 39 of
them, leaving a state score computed against one structure. Without a word, that
output is indistinguishable from one computed against forty-five.

Any state whose ensemble falls below five usable references now adds a note to
the analysis payload and the log, saying how many it rests on and that the
scores are indicative rather than measured. Five is the smallest number that can
show a spread rather than a pair of points, and both worked examples clear it
several times over, so it only fires on an uploaded protein whose catalogue
coverage is thin.

The alternative was refusing to score at all below some count. Rejected: a
one-reference comparison is still the most informative thing available for that
protein, and a refusal would throw away a usable answer rather than qualify it.
The spec's own rule on this (0.1 rule 3) is to degrade gracefully and flag, not
to stop.
