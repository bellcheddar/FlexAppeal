# Reference states: what the build found

Results from `FLEXAPPEAL_REFSTATES_BUILD_SPEC.md`. Every number here was measured
during the build rather than quoted, and the measurement is named beside it so it
can be re-run. Where a result contradicts the spec's expectation, the spec's
expectation is stated too.

Units are angstrom and nanoseconds throughout.

---

## 1. Results that hold before any simulation ran

The reference ensembles and the state calls made on them are evidence in their
own right, and three of them came out before the first production frame.

### 1.1 The DFG two-distance criterion has the wrong threshold

The published Dunbrack criterion, as implemented in BoltzMaker's
`sse_comparison`, measures the DFG-Phe ring to two anchors and asks which side
of 11.0 Å each distance falls. Against ABL1's own KLIFS reference set that rule
agrees with KLIFS on **11 of 17** structures.

All six disagreements are KLIFS `out` called `other`, and all six sit between
11.1 and 11.3 Å from the catalytic Lys: a tenth of an angstrom past a line drawn
at 11.0. 2HYY and 3PYY are among them, and both are plainly DFG-out.

The classes do not overlap. On the difference between the two distances:

| class | αC anchor | Lys anchor | difference |
|---|---|---|---|
| DFG-in | 6.6 to 6.8 | 12.9 to 14.0 | −6.1 to −7.9 |
| DFG-out | 13.1 to 16.3 | 9.1 to 11.3 | +2.4 to +5.8 |

An **8.5 Å gap** separates them, with nothing inside it. Reading the state from
which anchor the ring is nearer needs no threshold where the data is dense, and
agreement goes to **17 of 17**. A margin is kept, but placed in the empty region
rather than the populated one: within 1.0 Å of equidistant the answer is `other`,
which matters for trajectory frames rather than references, because a DFG
mid-flip genuinely is neither state and a bare sign test would hide the
transition.

![ABL1 reference state map](data/figures/abl1_map.png)

The figure is the argument. Both dashed lines are the published 11.0 Å cut,
which is a two-distance criterion and so applies to each axis. On the αC axis it
falls in an empty six-angstrom gap and separates the classes cleanly. On the Lys
axis it cuts straight through the DFG-out cluster: six red points sit between
11.1 and 11.35 Å, just above the line, and those are exactly the six structures
the rule calls `other`. The shaded path is `abl1_apo`, which stays inside the
DFG-in cluster throughout.

Logged as D010. The same threshold is in BoltzMaker's own `compare-sse` and will
mislabel the same structures there, so this is a candidate upstream change.

### 1.2 The β2AR ionic lock carries no categorical information

Measured on the two start models, between guanidinium and carboxylate heavy
atoms rather than between one named atom of each symmetric group:

| | TM3–TM6 (3.50–6.34) | ionic lock | Y5.58–Y7.53 OH |
|---|---|---|---|
| active, 4LDE | 14.68 | 15.47, broken | 4.81 |
| inactive, 5D5A | 8.35 | 9.95, broken | 11.26 |

The cleft and the tyrosine bridge both separate the states cleanly and in the
expected direction. The ionic lock reads **broken in both**.

That is the real β2AR result rather than a fault. The R3.50–E6.30 lock is known
not to form in the inactive crystal forms of this receptor. The consequence for
this build is specific: for β2AR the lock's categorical rule corroborates
nothing, so the state label falls back to geometry alone, while the lock's
*distance* still separates the states by 5.5 Å and is worth plotting. A scorer
that treated the categorical label as evidence would be agreeing with itself.

`R3.50` resolves to `R131` from the GPCRdb residue table, which is section 5.5's
own worked example, so the generic numbering is wired up correctly.

### 1.2b The refined models earn their place, measurably

16 of the 116 GPCRdb references cannot have their TM3–TM6 cleft measured at all,
because the residue at 3.50 or 6.34 in the deposition is not the one the
catalogue names. These are engineered constructs, and the identity check refuses
them rather than measuring a distance between the wrong two residues.

**4LDE is one of the 16**, and 4LDE is the entry the active start model comes
from. Its *refined* model, which is what section 4.3's ladder actually uses,
measures fine: the refinement reverts the thermostabilising mutations and
rebuilds the missing segments, so the residue at 3.50 is the one GPCRdb says it
is. The choice of rung was made on the spec's instruction rather than on
evidence, and this is the evidence: on the raw deposition the headline
activation measure is unavailable.

### 1.3 The orientation fit already shows the activation difference

Both start models were superposed onto OPM's 2RH1, an inactive-state entry, on
TM Cα only. The inactive model lands at **0.20 Å** and the active at **2.51 Å**.
That 2.3 Å gap is the difference this build exists to measure, visible before
any simulation.

### 1.5 The R-spine threshold has the same fault as the DFG one

Section 5.4 calls the regulatory spine `assembled` if every consecutive pair of
side-chain centroids is within 6.0 Å. Measured on ABL1's own reference set, the
pair that carries the call is DFG-Phe to HRD-His, and it separates the states
completely:

| class | DFG-Phe to HRD-His centroid | inside 6.0 Å |
|---|---|---|
| active | 5.33 to 5.89, median 5.70 | 12 of 12 |
| inactive | 9.53 to 11.08, median 9.86 | 0 of 45 |

So 6.0 Å is a correct threshold for crystal structures. It is also a tenth of an
angstrom above the top of the active range, which leaves no room for a side chain
that is breathing, and that is what broke it: frames of the apo DFG-in trajectory
read **6.03 to 6.69 Å**, plainly on the active side of a 3.6 Å gap, and every one
of them was called `broken`. The whole-protein score then reported a standing
`conflict` on the DFG motif, where the geometry said active at s = +0.755 and the
rule said inactive, for 43 frames out of 43.

The threshold is moved into the empty region between the clusters, at 7.7 Å,
which is the same correction D010 made to the DFG criterion and for the same
reason: a crystal structure is one minimised conformation and an MD frame at 300
K is not. Checked both ways afterwards, the references still classify perfectly
(5 of 5 active `assembled`, 14 of 14 inactive `broken`) and the trajectory's
standing conflict disappears.

The crystal value is kept beside it as `R_SPINE_CONTACT_CRYSTAL_A` so a
structure-only comparison can still use the published number.

### 1.4 Positive controls

| control | result |
|---|---|
| P1 DFG label against KLIFS | 17 / 17 after D010 |
| P2 αC label against KLIFS | 17 / 17 |
| P3 TM3–TM6 state separation, leave-one-out | 100 / 100 at an 11.39 Å threshold |
| R-spine label against KLIFS state, after recalibration | 5 / 5 active, 14 / 14 inactive |
| P4 each start structure scores as its own state at t = 0 | 4 / 4 |

#### P4 in full

| row | expected | whole-protein s | verdict |
|---|---|---|---|
| abl1_apo, DFG-in from 2GQG | active | **+0.305** | pass |
| abl1_imatinib, DFG-out from 2HYY | inactive | **−0.368** | pass |
| b2ar_apo, active refined from 4LDE | active | **+0.593** | pass |
| b2ar_carazolol, inactive refined from 5D5A | inactive | **−0.934** | pass |

Measured on the start structures themselves, with no simulation in between: the
start model is both topology and single frame. A scorer that cannot place its
own start model has nothing to say about a trajectory that began there.

**The kinase margins are thin and that is worth stating.** Both ABL1 rows clear
the ±0.3 band by less than 0.07, where the receptor rows clear it by 0.29 and
0.63. The reason is interpretable rather than worrying: ABL1's active and
inactive references differ mainly at the DFG motif and the activation loop,
while the hinge, the catalytic loop and the glycine-rich loop sit at nearly
identical distances from both ensembles. Those motifs contribute real residue
weight and almost no signal, so the residue-weighted mean is diluted towards
zero. The per-motif breakdown is where the kinase signal actually lives: xDFG
reads +0.755 on the apo start while the hinge reads −0.078.

That is an argument for reading the motif table rather than the headline number
for kinases, and the explorer shows both side by side for that reason. It is
also a limitation to carry into any future weighting decision: a whole-protein
mean over motifs of unequal discriminating power understates a real state call.


---

## 2. Measurement choices the data forced

### 2.1 The kinase superposition core is the KLIFS scaffold, not αE/αF/αH

Section 5.4 names the C-lobe helices αE, αF and αH. Only αE is inside the KLIFS
pocket, so taking the spec literally gave a core of **five residues**, which is
not a rigid-body fit. With it, motif RMSDs to the inactive references came out at
17 to 22 Å.

The core is now the KLIFS pocket scaffold, every pocket position outside the
mobile and functional regions, which is the choice BoltzMaker's own
`alignment.stable_reference_positions` makes. This is a deviation from the letter
of the spec, so it is measured rather than asserted: **43 residues superposing
across 19 references at 0.456 Å mean and 0.856 Å worst**. A core that moved with
activation would show up there as a large spread.

The receptor core needed no such deviation. Section 5.5's x.38 to x.58 span of
TM1 to TM4 gives **82 residues across 100 references at 0.829 Å mean, 1.289 Å
worst**.

### 2.2 KLIFS's `a.l` region is not the activation loop

KLIFS's own `a.l` is two pocket positions, because the pocket alignment stops
where the loop leaves the pocket. Using it as the activation loop gives a
five-residue span, and every A-loop RMSD measured on it would be an RMSD of the
DFG motif under another name. The loop is therefore taken from the DFG aspartate
to the APE motif found in the chain's own sequence, which gives 26 residues for
this construct.

### 2.3 Rule direction belongs to the metric, not to the word

The two vocabularies collide. `broken` means inactive for the R-spine, whose four
side chains stack in the active kinase, and active for the GPCR ionic lock, which
holds the inactive receptor shut. A single table keyed on the word alone labelled
every broken lock inactive and then agreed with the geometry for the wrong
reason.

---

### 2.4 A motif's residue list has to be agreed before it can be compared

An RMSD between two structures over different residue subsets is not a
comparison, so a motif is measured on exactly the residues both sides resolve.
Applied naively that rule deletes the motifs that matter: TM6 is 38 residues
through the most mobile helix of the receptor, no deposition resolves all of
them, so every reference failed the check and the motif went unscored. The
kinase activation loop went the same way. Both are the motif that carries the
state of its protein.

Intersecting across every reference is worse. On 116 GPCRdb entries it leaves
every motif empty, because a handful of partial depositions decide the answer
for all of them.

What works is to count coverage over the references that contain the motif at
all, and keep the residues 90 per cent of those resolve. Helix 8 is the case
that forces the distinction: absent or disordered in most β2AR structures, it
keeps no residues measured against all of them and 11 of 13 measured against
the ones that have it. The ensemble size travels with each motif, so a reader
can see that H8 is scored against fewer structures than TM6 is.

| motif | residues kept | of | references with it |
|---|---|---|---|
| TM5 | 35 | 42 | 107 |
| TM6 | 31 | 38 | 107 |
| TM7 | 25 | 25 | 107 |
| H8 | 11 | 13 | fewer |
| NPxxY | 5 | 5 | 107 |

**The same operation is wrong for the kinase**, which is the more useful half of
this. ABL1's references are split between the ABL1a and ABL1b numbering, 19
residues apart, so intersecting residue *numbers* compares residue 290 in one
structure against 309 in another: the intersection collapsed and five of six
motifs were dropped. The receptor's references all share one numbering, which is
the only reason it is safe there. Each kinase reference therefore keeps its own
numbering, and the activation loop is compared on the same count of residues
from the DFG aspartate that both spans share.

### 2.5 The receptor path's first end-to-end result is a control

Run on the built active system, every motif reads active, and the nearest
references are the ones they should be:

| motif | s | nearest active | nearest inactive |
|---|---|---|---|
| NPxxY | +0.839 | 3P0G, 0.27 Å | 8W1V, 3.14 Å |
| TM6 | +0.778 | 9LW5, 0.54 Å | 2RH1, 4.33 Å |
| TM7 | +0.742 | 3P0G, 0.33 Å | 6PS0, 2.25 Å |
| H8 | +0.555 | 3P0G, 0.35 Å | 6PS0, 1.23 Å |
| TM5 | +0.298 | 9BUY, 1.21 Å | 3PDS, 2.23 Å |

Whole-protein score +0.593, active, no conflicts. The system was built from
4LDE, an active-state refined model, so active is the answer it should give;
3P0G is an agonist-bound active entry and 2RH1 is the canonical inactive one.

The inactive system is the other half of the control, and it is the stronger
one. Built from 5D5A, it scores **−0.934, inactive**, every motif inactive, with
5D5A itself as the nearest inactive reference at **0.00 Å**, which is what it
must be, since that is the structure the system was built from. The active
distances are 2 to 4 Å. Secondary structure agrees independently: every TM helix
comes out 92 to 97 per cent helical and NPxxY 100 per cent.

So the receptor pipeline reads active for an active start and inactive for an
inactive one, with no conflicts either way.

### 2.6 A criterion that cannot tell the states apart should not be asked to

TM6 originally reported `conflict` on every frame of the inactive run, with the
geometry at s = −1.000 against a nearest inactive reference of 0.00 Å. The cause
was attaching the ionic lock as TM6's corroborating rule, which section 1.2
had already shown reads `broken` in both β2AR reference states. Mapping `broken`
to active then guaranteed a conflict on exactly the frames the measurement had
got right.

TM6 is now scored on geometry alone for this receptor, and the lock's distance,
which does separate the states by 5.5 Å, is reported beside it as a measurement
rather than used as a vote. This is D015 one step further: D015 moved a
threshold that was in the wrong place, and this removes a criterion that has no
discriminating power for this protein at all. A permanent conflict on a correct
measurement is worse than no corroboration, because it trains a reader to
ignore the flag that exists to be noticed.

## 3. Faults worth recording because the output looked right

Three bugs in this build produced plausible numbers. They are listed because
the plausibility is the point: in each case nothing in the output indicated a
problem.

**References were each in their own crystallographic frame.** Scoring superposes
a trajectory frame onto a core and compares it against every reference, but the
references come from separate mmCIFs in unrelated frames, and only the one that
happened to be the alignment target was comparable. Motif RMSDs to the inactive
set read 17 to 22 Å while the target reference read 1 to 3 Å, and the
whole-protein score was a confident **0.877 "active"** for every frame with a
believable per-motif spread. After aligning every reference onto one frame on its
core, the same frame reads **0.173 "intermediate"**, every motif RMSD falls
between 0.16 and 4.11 Å, and the structure of the answer is interpretable: the
hinge is 0.18 Å from active and 0.16 Å from inactive because the hinge does not
distinguish the states, and xDFG is the motif that does, at +0.755.

**Trajectory residues had lost their deposited numbering.** OpenMM renumbers from
1, so the ABL1 construct's residues 233 to 498 arrive as 1 to 262 while every
annotation addresses deposited numbers. A naive offset would also have been
wrong: the construct spans 266 residues and contains 262, so trajectory residue
100 is deposited 333, not 332. The mapping now checks the one-letter sequence at
every position and refuses rather than returning a map that is wrong in the
middle.

**A restraint measured a wrapped coordinate.** Section 5.2 asks for lipid
phosphorus to be held in z while the bilayer relaxes, written as
`k*(z - z0)^2`. OpenMM wraps coordinates into the periodic box; the lower
leaflet sits at negative z and wraps by one box vector, so the restraint saw an
11.77 nm deviation on 238 atoms that had not moved. At 4184 kJ/mol/nm² that is
about 138,000,000 kJ/mol. Written as `periodicdistance(0, 0, z, 0, 0, z0)²` the
same system reads −2,213,385 with both restraints at exactly zero. The protein
restraint needs the same treatment for the same reason: written as a plain
harmonic it reported 5,045,427,552 kJ/mol on atoms sitting exactly on their
reference positions.

**And a bad comparison read as a result.** That +136,000,000 was attributed to
single precision for several hours, on the strength of a CPU run reading
−2,060,000 for what was described as the same system. The two runs did not
differ in precision: one had no restraints and the other had both. Nothing
measured in this build ever showed the two platforms disagreeing, and once the
real cause was fixed they agreed exactly. Four diagnostics went after the wrong
cause on the strength of it. The lesson is narrow and worth keeping: a
difference between two runs is only evidence about the variable that actually
differed.

**Membrane quality measured across the periodic boundary.** The bilayer QC
worked on the built system and was wrong on the trajectory, which is the worst
place for a difference like that to live, because the built system is the one
you check by hand. A trajectory stores wrapped coordinates and under a membrane
barostat the bilayer drifts until it straddles the box edge. The phosphorus z
values are then bimodal at *both* ends of the box with nothing in between, so
the median used to split the leaflets lands in the empty middle:

| | measured | actually |
|---|---|---|
| P–P thickness | 78.01 Å | 37.55 Å |
| waters in the hydrophobic core | 8,622 | 4 |

The 78 Å is the long way round the box between two leaflets 37.6 Å apart, and
the water count is bulk water measured against a midplane that is not where the
midplane is. Both numbers are the kind that would have gone into the deliverable
unchallenged: 78 Å is wrong but not absurd for a membrane, and a four-figure
core water count reads as a leaky bilayer rather than as an artefact.

Fixed by centring the bilayer on a circular mean before anything is measured,
and by taking both the thickness and the water separation as periodic distances.
The arithmetic mean is what caused it: the mean of two clusters at opposite ends
of a periodic axis is the empty middle.

**`addMembrane` believed a placeholder unit cell.** gemmi writes a `CRYST1` of
1 1 1 for a `P 1` structure with no cell. `addMembrane` takes the box height as
the larger of the lipid patch and the existing cell plus padding, so with a cell
present at all it took 7 nm from the patch rather than sizing from the receptor,
which is 9.4 nm tall. The receptor wrapped through the periodic boundary onto
itself. The error named a NaN coordinate and nothing about a box.

---

## 4. Deviations from the spec, with reasons

| Section | Spec | What was done | Why |
|---|---|---|---|
| 5.1 | 80 to 110 k atoms | 188 to 195 k | ICL3 is modelled hanging 35 Å below everything else and sets the box height; the lateral box is two lipid patches at any usable padding. Section 5.1's prescribed response to over 150 k, a padding reduction to 1.2 nm, is applied. D012. |
| 5.2 | mixed precision, never single for production | OpenCL single | Apple's OpenCL is single-precision only. The measured alternative was the CPU platform at 29.5 against 206 ns/day. Energy drift is measured rather than assumed. D010. |
| 5.2 | 4 fs step throughout | first 100 ps of heating at 1 fs | A freshly minimised bilayer does not survive 4 fs at full temperature. The 500 ps of NVT heating is kept whole. |
| 5.3 | 3 replicas, 100 ns, 4 rows | 1 replica, 50 ns, 4 rows | D009, agreed with Marc. |
| 5.4 | core is αE, αF, αH | KLIFS pocket scaffold | See 2.1. Verified at 0.456 Å. |

---

## 4b. Page QC, and the eleven failures it found

Section 10's Playwright and axe-core checks were initially skipped: the
environment had no `playwright` and the spec said never to run
`playwright install`. Marc pointed out that the constraint was not his, so both
were installed and the checks run properly. They were worth running.

**Final state: 25 tests pass.** Eight pages at 1440x900 and 390x844 with no
console errors, eight pages through axe-core with **zero serious or critical
violations**, and one test asserting the three Plotly panels actually draw
rather than leaving an empty div behind.

**Eleven findings, every one pre-existing or mine, none visible from the server
side.**

| Finding | Impact | Where |
|---|---|---|
| `aria-required-children` | critical | `role="tablist"` on a radio group |
| `color-contrast` x4 | serious | navigation hints, footer, grade meter, terminal clips |
| `definition-list` / `dlitem` | serious | `dt`/`dd` outside a `dl`, then a `p` inside one |
| `link-name` | serious | Mol*'s own logo link |

### The pattern worth keeping

Four of the contrast failures were the same mistake in four places: `opacity` on
small text. `--md-text-muted` was already 4.26 against white, under the 4.5 AA
asks for, and fading it made it worse:

| element | opacity | computed | ratio |
|---|---|---|---|
| `.md-tab-hint` | 0.8 | `#8996a9` | 3.00 |
| `.md-footer-meta` | 0.8 | `#798799` | 3.41 |
| `.md-grade-score em` | 0.65 | `#909cab` | 2.59 |

Each was fixed by deleting the opacity and letting the colour token say what it
means, and the token itself was darkened to `#5a6b80`, which is 5.46 on white
and 5.08 on the alt background. **Opacity on text is the finding**: it makes the
rendered colour something no token states, so no amount of reading the
stylesheet reveals it.

### A fix that was the wrong shape

The first pass at the terminal clips patched the two colours axe had named. That
was wrong, and for an instructive reason: axe samples nodes, so its report is
evidence that a colour fails and not evidence that the others pass. Three more
failing colours were still there afterwards, including `#800080` at **1.91**.
`tools/fix_clip_contrast.py` now audits every colour in every clip against the
threshold and exits non-zero if one is missed. Only colour changes; the recorded
text is untouched.

### What the browser caught that a server could not

The `link-name` violation on Mol*'s logo was fixed, then came back: the library
builds its viewport controls after the structure-loading promise resolves, so
labelling it once found no logo. Only a real browser shows that, and only after
the page has settled. The retry is bounded rather than an observer left running.

## 4c. The QC checklist

Section 10's list, with where each one stands. Nothing here is marked done on
the strength of having been attempted.

| Item | State |
|---|---|
| Numbering maps, KLIFS position to residue | done. 3PYY resolves at the ABL1b +19 offset with nothing typed |
| Numbering maps, generic number to residue | done. `R3.50` resolves to `R131`, which is section 5.5's own worked example |
| Metrics on synthetic geometries | done. A salt bridge at 3.5 Å measures 3.5 Å; a helix rotated 0, 15, 30 and 45 degrees reports each angle |
| State score symmetry | done. +1 on the active reference, −1 on the inactive, exactly 0 halfway, for any separation |
| Membrane builder smoke test | done, and on a real run rather than 5 ps: the Tier 1 row is 1 ns of production |
| Resume from checkpoint | done. A reloaded checkpoint continues the clock, and the resume branch is asserted to come before the equilibration branch |
| P1 to P6 | P1 to P4 pass. P5 and P6 need the Tier 2 rows |
| No metric all-NaN for a motif | done, and it is why TM6 and the activation loop were caught: both were silently unscored rather than NaN, which is the same fault wearing a different mask |
| Playwright at both viewports | done. 8 pages at 1440×900 and 390×844, no console errors |
| axe-core, zero serious or critical | done. 8 pages, zero. Eleven violations found and fixed getting there |
| Performance targets | done. First paint 396 ms median against 2 s; bundle 0.4 MB against 150 MB |
| Units everywhere | angstrom, nanoseconds, degrees, Å². Threshold lines labelled on the state map and the timeline |
| Energy drift within limit | done. Within noise at 3σ on both runs; the 2 fs fallback is unused |

Two things on that list were found to be wrong by writing the check rather than
by reading the code, and both were mine: the drift criterion at one sigma, and
the paint test's race. Neither would have been visible from inspection.

### 2.7 The receptor state map had only one axis

Section 5.6 asks for the GPCR map to be the TM3-TM6 cleft distance against the
NPxxY RMSD to inactive, which is the plane Dror et al. draw their state maps on.
Only the first of those was available: the cleft distance is in the phase 1
manifest and the NPxxY RMSD was nowhere, so the endpoint returned 100 points
with `y: None`.

That is a failure that reports success. The endpoint returns 200, the panel
renders, and nothing is drawn in it, which is the same shape as the empty-panel
problem found by the browser QC except that no status code reveals it.

The second axis is now computed from the reference coordinates during the
analysis, since they are already aligned into one frame, rather than being added
as another manifest field to keep in step. A reference is not compared with
itself, which matters: without that exclusion every inactive structure would
sit at zero and the map would be two lines rather than two clouds.

Both axes separate the states cleanly:

| state | n | cleft distance | NPxxY RMSD to inactive |
|---|---|---|---|
| active | 71 | 15.84 Å | 2.79 Å |
| inactive | 22 | 8.25 Å | 0.15 Å |

A point missing either axis is left out rather than sent as null, and a map with
nothing left returns 404 so the panel is hidden rather than blank.

### 2.8 One failure shape, three times

Three defects in this build were the same thing wearing different clothes: an
endpoint returned 200 with data nothing could draw.

| Where | What was returned | What was drawn |
|---|---|---|
| state map panel | a 404 the page had already rendered a panel for | an empty box, plus a console error |
| state map data | 100 points with every `y` null | nothing, on a 200 |
| motif series | 8 metric rows and no time axis | an invisible line, on a 200 |

None of them shows in a status code, none raises, and none is visible from the
server side without asking what the numbers mean. The third was found by
constructing the case rather than by waiting for it: a motif that resolves its
inactive distance and not its active one is **legal**, because the rule
throughout this build is that an unmeasured key is omitted rather than written
as zero, so every reader has to cope with it. The series was keyed on the
`d_act_a` row's arrival, so when that row never came the time axis was empty
while the metrics were full.

The invariant is now asserted directly rather than instance by instance: a 200
from any of these endpoints must carry arrays that are non-empty and the same
length as each other, and an endpoint with nothing to show must 404 so the panel
is hidden. A panel is allowed to be absent. It is not allowed to be present and
blank.

## 4d. What was outstanding, found by asking

Marc asked whether everything was ready for the results landing and said he
found it hard to believe nothing was. He was right. Checking rather than
asserting turned up three things, one of them a bug that would have corrupted a
whole row silently.

**`traj.load` selected `protein or resname CAU`.** CAU is carazolol. The ABL1
imatinib row's ligand is STI, so that row would have been loaded without its
ligand, analysed as an apo system, and written a viewer bundle with no ligand in
it, all without raising. The selection is now an exclusion: protein plus
anything that is not solvent, counter-ion or bilayer. Written that way because
the inclusion version is exactly the wrong shape for a package whose claim is
that it handles a protein it has never seen. Verified on all three systems: STI
kept for imatinib, CAU for carazolol, nothing extra for the apo runs.

**P5 and P6 did not exist.** Both are section 9 acceptance criteria and both
read the explorer database, so neither could be written against data that had
not arrived; neither had been written at all. They now exist, report `pending`
rather than nothing when there is no data, and are exercised against fixtures so
that when the rows land the only new thing is the data. P5's adaptation is
recorded rather than reinterpreted: section 9 asks for 2 of 3 replicas and D009
dropped the replicas to one, so the rule is per run, because one replica cannot
support a majority.

**And I had deleted four functions.** Rewriting `traj.load` replaced everything
from that function to the end of the file, taking `deposited_numbering`,
`renumber`, `scaffold_core` and `core_agreement` with it. All 56 tests still
passed, because none of them imports those four; the loss surfaced only on
running the real tool. There is now a test that reads what every script imports
from the package and asserts each name resolves, and it was confirmed to catch
the deletion by making one.

**A resume would have duplicated frames and shifted every later timestamp.**
Checkpoints were written every nanosecond and frames every 100 picoseconds, so
the checkpoint sat up to ten frames behind the trajectory. A resume restarts
from the checkpoint and the DCD reporter appends, so those frames get written a
second time. Because frame time here is the frame index times the interval,
duplicates do not merely repeat a point: they shift the timestamp of every frame
after them.

`abl1_apo` had exactly this waiting for it, 43 frames against a checkpoint at
4.0 ns, and it is the row the queue will resume. The membrane rows take 20 hours
each on a machine that sleeps, so resumes are expected rather than exceptional
and the error would have accumulated.

Two changes. The checkpoint interval now equals the frame interval, so the two
stay in step. And a resume truncates the trajectory and the state table back to
the frames the checkpoint actually covers, which makes it exact rather than
nearly right. Verified on the real `abl1_apo` trajectory: 43 frames to 40, the
kept frames bit-identical to the originals, box vectors preserved, and a no-op
when the checkpoint is not behind. The truncation reads only the frames it
keeps rather than loading the file, which matters at a gigabyte.

One concern turned out to be unfounded and is recorded because checking it was
the point: the final membrane pass loads a 50 ns solvated trajectory in one go,
which is 1.16 GB against 128 GB of RAM. Left alone rather than "fixed".

## 4e. Eight metrics that were not implemented, and one that was a parameter choice

An audit of sections 5.4 and 5.5 against what existed found eight named metrics
missing: the alphaC salt bridge, alphaC displacement, the P-loop, the hinge
backbone flags, TM6 rotation, the PIF connector, helix 8's angle and the TM
kinks. All are now implemented, and checking them against the start models
found a methodological fault worth more than the metrics themselves.

### The kink angle was the window, not the helix

The first implementation fitted a principal axis straight to the Calpha either
side of each conserved proline. It gave **P6.50 at 62.7 degrees for the active
structure against 50.2 for the inactive one**, a 12.5 degree difference in the
direction the literature describes. It would have been written up.

It was wrong. Sweeping the only free parameter, how many residues to fit,
produced:

| residues per half | P6.50 | P7.50 |
|---|---|---|
| 3 | 24.3 | 23.0 |
| 4 | 50.2 | 36.0 |
| 5 | 29.2 | 14.4 |
| 7 | 25.7 | 8.3 |
| 10 | 33.2 | 10.2 |

The answer was the window. Calpha positions spiral around the helix axis, so a
principal axis through a short segment is dominated by where in the turn it
starts and ends rather than by the axis. Fixed by averaging over one full turn
before fitting, which removes the spiral, and that is asserted against an ideal
helix built to the standard geometry: turn-averaging holds the axis on z across
every usable window while the raw fit does not, by more than a factor of ten.

**And the honest reading changed with it.** With the window swept and the spread
reported:

| | P6.50 | P7.50 |
|---|---|---|
| active | 30.8 ± 3.5 | 7.4 ± **17.2** |
| inactive | 34.1 ± 5.2 | 14.3 ± 5.3 |

The spreads overlap, so **neither kink distinguishes these two structures**, and
P7.50's active value is not determined at all since its spread exceeds twice its
value. That is a negative result where an hour earlier there was a confident
12.5 degree activation signal pointing the way the literature says it should.

This is the most dangerous shape a measurement error takes: wrong, and in
agreement with expectation. A metric now reports `resolved: false` when the
windows disagree by more than half the value, rather than letting a median stand
alone.

### What the new metrics do separate

| measure | active start | inactive start |
|---|---|---|
| alphaC salt bridge | 3.34 Å, in | 3.41 Å, in |
| TM6 axis rotation, one against the other | **25.4°** with a **7.56 Å** centroid shift | |
| PIF connector RMSD, one against the other | **3.29 Å** over 14 shared side-chain atoms | |
| H8 angle to the membrane plane | 20.4° | 16.4° |

Those two figures are larger than the first ones measured, 20.9° and 3.13 Å,
and the difference is the core superposition. Without it the whole-protein
misalignment was partly cancelling the helix motion, so the uncorrected numbers
*understated* the signal. The fit has to be onto the structure being compared
against, which cost a round of debugging: passing the transform the scoring uses
instead, which maps a frame into the reference ensemble's frame rather than the
start model's, gave a 78 Å centroid shift and an 82 Å connector RMSD. Both
obviously wrong, which was luck. The invariant that is not a matter of luck is
now asserted: a structure compared with itself reports exactly zero.

25.4° of axis rotation with a 7.6 Å centroid shift for the 6.28 to 6.40 span is
the TM6 outward swing itself, and it is consistent with the roughly 14 Å the
cytoplasmic tip is described as moving.

The alphaC salt bridge reading `in` for both is correct rather than a failure:
imatinib-bound ABL1 is the classic DFG-out, alphaC-in type II conformation, and
19 of 20 KLIFS references for this kinase are alphaC-in.

## 4f. Four functions that were implemented, tested, and called by nothing

`smooth`, `transitions`, `sse.agreement` and `sse.reference_profile` all existed
with tests and no caller. The SSE comparison against the reference ensembles is
an explicit section 5.6 requirement, so that was a gap in the deliverable rather
than unused convenience code.

Wired up, it is an independent line of evidence and it agrees with the geometry.
On the Tier 1 receptor run:

| motif | active | inactive | margin | lean |
|---|---|---|---|---|
| TM7 | 0.712 | 0.561 | 0.151 | active |
| TM5 | 0.982 | 0.832 | 0.150 | active |
| TM6 | 0.934 | 0.868 | 0.066 | active |
| NPxxY | 0.739 | 0.774 | 0.035 | **neither** |
| H8 | 0.864 | 0.888 | 0.025 | **neither** |

**Scored across every residue, it gave two different structures the same
answer.** Running the ligand-bearing receptor path for the first time exposed
it. The active-start and inactive-start systems, built from 4LDE and 5D5A, read:

| motif | active run | inactive run |
|---|---|---|
| TM6 | 0.934 / 0.868 → active | 0.934 / 0.868 → active |
| H8 | 0.864 / 0.888 → neither | 0.864 / 0.888 → neither |
| TM5 | 0.982 / 0.832 → active | 0.953 / 0.813 → active |

TM6 and H8 bit-identical for two different proteins in two different states.
The reason is simple once seen: a transmembrane helix is helical in both states,
so its SSE string is the same in each and so is its agreement. The lean was
reporting that the active reference set is more uniformly helical than the
inactive one, which is a fact about the catalogue rather than about the frame.

The comparison is now scored only over residues where the two reference sets
disagree about the dominant class, with at least three of them required. That
line is drawn from the data: TM5 has two discriminating residues and gave the
identical answer for both structures, TM7 has seven and tells them apart
correctly in both directions.

| run | built from | SSE call |
|---|---|---|
| β2AR apo | active 4LDE | TM7 = **active**, 7 residues |
| β2AR carazolol | inactive 5D5A | TM7 = **inactive**, 7 residues |
| ABL1 apo | DFG-in 2GQG | A-loop = **active**, 14 residues |

One call per run instead of five, every one correct, and on the motifs that
genuinely change secondary structure on activation. The test asserts the
property the first version lacked: it is not enough that a frame matching the
active set reads active, the frame matching the inactive set has to read
inactive.

**The margin threshold is still needed on top of that.** Secondary structure is far less
discriminating than geometry: most residues of a transmembrane helix are helical
in both states, so every agreement here is high and the differences are small.
The first version reported H8 as leaning inactive on a margin of 0.024, which is
reading noise and would have sat on the page as a finding. A lean is now only
called above a 0.05 margin, chosen from the observed spread rather than picked:
the discriminating motifs separate by 0.07 to 0.15 and the ones carrying no
signal by 0.02 to 0.04, so the line goes between them.

The three TM helices leaning active corroborates the geometric score of +0.277
for the same run, from a measurement that shares none of its machinery.

### The kinase result is the better one

Two independent measurements agree on *which motifs carry the signal*, which is
a stronger statement than either agreeing on a state.

| motif | SSE active | SSE inactive | margin | SSE lean | geometric s |
|---|---|---|---|---|---|
| A-loop | 0.285 | 0.129 | 0.155 | active | +0.136 |
| xDFG | 0.600 | 0.482 | 0.118 | active | **+0.755** |
| c.l | 0.950 | 1.000 | 0.050 | inactive | −0.087 |
| alphaC | 0.945 | 0.903 | 0.043 | neither | +0.305 |
| g.l | 0.408 | 0.423 | 0.015 | neither | +0.019 |
| hinge | 0.667 | 0.667 | **0.000** | neither | −0.078 |

The motifs that discriminate geometrically, xDFG and the activation loop, also
discriminate by secondary structure. The hinge discriminates by neither: its
geometric distances are 0.18 Å from active and 0.16 Å from inactive, and its SSE
agreements are identical to three decimal places. Two measurements with nothing
in common reaching the same conclusion about which parts of this kinase hold
state information is the most useful thing in this section.

Both A-loop agreements are *low*, 0.285 and 0.129, which is not a fault: the
activation loop varies between structures of the same state, so no reference set
agrees with itself closely. It still separates them, which is what the margin
measures.

### A cache that is a requirement, not an optimisation

Running DSSP over the 116 receptor references takes about nine minutes, and
section 2.3's streaming analysis re-runs every half hour while production is
going. The labels are properties of deposited structures and never change, so
they are cached: nine minutes becomes seven seconds. The key covers the residues
asked for as well as the structure codes, so changing a motif's span invalidates
the cache rather than silently reusing labels for different residues.

### And a lost edit, worth recording as a process failure

Three replacements were applied to one file by a single script whose last
assertion failed. Python raised before the write, so all three were lost
together, not just the one that failed. The run afterwards therefore used the
old code and reported no agreement at all, and I initially mis-diagnosed that as
a logic fault in `motif_agreement`, which worked perfectly in isolation when
tested. The fix was to reapply and then verify each part had landed, rather than
trust that an assertion passing meant the file had changed.

## 4g. Three numbers nobody could re-derive

An audit for unreferenced public functions turned up six. Three were dead and
are now deleted. The other three were the interesting half:
`refs.prepare_kinase_reference`, `refs.write_manifest` and `gpcr.fit_threshold`
were all reachable only by reading them, because the code that drove them was a
one-off in a session and was never committed.

What they had produced was not scratch output. `data/manifests/refs_abl1.jsonl`
and `refs_b2ar.jsonl` feed the reference browser and both state maps, P1 and P2
are counted off the kinase manifest, and P3's 11.39 angstrom threshold sits in
`tools/ingest.py` as a constant. So four of the spec's own acceptance numbers
rested on artefacts that could not be regenerated. That is a reproducibility
gap and not a tidiness one: a number nobody can re-derive is a number nobody
can check, including me.

`tools/build_references.py` closes it. `--verify` rebuilds everything in memory
and diffs it per record rather than per file, because the order a catalogue
returns rows in is not stable and a file diff would report a shuffle as a
change.

Both halves now reproduce exactly:

| manifest | rebuilt | on disk | changed | controls |
|---|---|---|---|---|
| `refs_abl1.jsonl` | 20 | 20 | 0 | P1 17/17, P2 17/17 |
| `refs_b2ar.jsonl` | 116 | 116 | 0 | P3 100/100 at 11.39 A |

and the script checks that `GPCR_CLEFT_THRESHOLD_A` still matches the fit rather
than leaving the two to drift apart silently.

### The first rebuild disagreed, and it was the script that was wrong

It reported 42 kinase rows against the 20 on disk, with P1 at 37/37 and P2 at
36/37. Both numbers pass their marks, which is what made this worth chasing
rather than accepting: a rebuild that passes is not the same as a rebuild that
reproduces.

KLIFS lists one row per chain and per alternate conformation, so ABL1 comes
back as 61 rows over 20 depositions, with 2G1T alone appearing four times. My
script deduplicated on `(pdb, chain)` where the original had deduplicated on the
deposition. Per-PDB first occurrence reproduces all 20 chains exactly,
including the six that are not chain A, so the original intent was unambiguous
once tested.

The per-entry choice is also the defensible one. Chains of one asymmetric unit
are copies in one crystal, not independent checks of whether our label agrees
with KLIFS's, so counting them weights a four-chain deposition four times and
quotes an n larger than what was actually measured. The wider set is kept as
`--all-chains`, which reports rather than writes, because it answers a different
question: whether our labels are stable *within* a crystal.

### What the per-chain check then found

One chain in 37 disagrees, and it is informative:

| 2HZ4 chain | Lys NZ to Glu carboxylate | our aC | KLIFS |
|---|---|---|---|
| A | 4.20 A | out | in |
| B | 2.71 A | in | in |
| C | 2.74 A | in | in |

The cutoff is 4.0 angstrom. Chain A misses it by 0.20 angstrom while B and C sit
at 2.7 with the bridge plainly intact, so this is a side chain modelled at the
edge of a threshold in one copy, not a systematic error in the call.

It is still worth stating plainly that the manifest uses chain B here, and that
first occurrence picked it rather than any judgement of mine. Had the order come
back the other way, P2 would read 16/17, or 94 per cent, against its 85 per cent
mark. The control passes either way, which is the useful part: the result does
not depend on which copy of this crystal the catalogue happened to list first.

A test now asserts that no deposition appears twice in the kinase manifest. The
invariant needs guarding precisely because breaking it is silent: the per-chain
rebuild passed both controls and only the record count gave it away.

### And they were never committed at all

Checked before pushing, and it corrects the paragraphs above. `data/manifests/`
is in `.gitignore`, so neither manifest is tracked: they are local files, and
the word "committed" in my own account of this was wrong.

That makes the gap worse rather than smaller. A fresh clone has no manifests at
all, so before `tools/build_references.py` existed the reference browser and
both state maps had no way to be populated on any machine but this one, and
P3's threshold in `ingest.py` was a number with no reachable derivation. The
fix is the same fix and it was the right one; the description of what it fixed
was too generous.

What is tracked under `data/` is five files: `references.json`, which the About
page reads, and four `*_start.json` carrying each row's provenance, which hold a
PDB code, a chain, a span and an expected state and no coordinates. Everything
else, the trajectories, the caches, the built systems, the manifests and the
explorer database, is ignored and travels by rsync or is rebuilt.

### The same gap, in a second place

Having found it once, the obvious question is which other committed numbers are
fitted rather than chosen. Going through the module constants, most are policy
or published contact distances: the 4.0 angstrom salt bridge, the 3.5 angstrom
sodium coordination, the decisive bands. Two were fitted to data. One was P3's
11.39 angstrom cleft threshold, now re-derived above. The other was
`R_SPINE_CONTACT_A`, D015's 7.7 angstrom, and it had the same problem: a
measurement quoted in a decision record with nothing able to reproduce it.

The kinase verify now recomputes it. The separation holds and D015's conclusion
is unchanged, but re-deriving it corrected two of its figures and found a case
it had missed:

| quantity | D015 stated | re-derived |
|---|---|---|
| active cluster, deciding pair | 5.33 to 5.89 A | 5.33 to 5.96 A |
| inactive cluster, deciding pair | 9.53 to 11.08 A | 10.35 to 12.92 A |
| empty region between them | 3.6 A | 4.39 A, from 5.96 to 10.35 |
| margin left by the published 6.0 A | a tenth of an angstrom | 0.04 A |

D015 had quoted the DFG-Phe to HRD-His pair, but the call is made on whichever
consecutive pair is widest, because one pair over the threshold is enough to
make the state `broken`. Measured that way the active bound is 5.96 rather than
5.89, so the published 6.0 angstrom clears the top of the active cluster by four
hundredths of an angstrom. That is a stronger argument for moving the threshold
than the one D015 actually made.

It also surfaced the reference that turns on the choice. 2G1T, the Src-like
inactive start, sits at 6.14 angstrom: `assembled` at 7.7 and `broken` at 6.0.
D015 counted 19 references and had left the Src-like ensemble out of that count,
so the one structure where the two thresholds disagree was invisible. It happens
to be the clearest illustration of the point, because 2G1T is DFG-in with the
alphaC helix out and its spine stack is genuinely intact, so `assembled` is
correct and the published threshold gives the wrong answer.

Two notes on what this does and does not establish. The control checks that the
threshold still sits in an empty region and that every crystal classifies
correctly, and it fails the run if either stops being true, so a catalogue
change that moved the clusters would be caught. It does not validate 7.7 against
anything other than these 20 depositions, and the trajectory evidence for it is
still the 43 frames D015 cites.

A stale docstring went with it: `r_spine_state` still said "within 6.0 A" after
D015 moved the threshold, so the code described a behaviour it no longer had.

## 4h. What was still untested with the results about to land

Asked once more whether everything is ready for the Tier 2 trajectories, the
useful move was to stop asserting it and go and run the paths that had never
run. Two had not.

**The receptor ingest path had never executed.** The database held one system,
one run and 20 reference structures, all of them ABL1. The 116 beta2AR
references, the five receptor motifs and the whole GPCRdb annotation path were
written, tested against fixtures and never once run against a real analysis.
Run against the existing receptor analysis on a scratch database it works: 5
motifs, 100 frame rows, 116 references annotated and 100 of 100 agreeing with
GPCRdb, with the second map axis merged onto the same records. Nothing needed
fixing, but "it is implemented" and "it has run" were not the same claim, and
only one of them was true.

**The bilayer QC had never written a row.** `membrane` was empty. The pass runs
at ingest rather than analysis time, which is why it is invisible in the
analysis JSON, and it needs `solvated.pdb` beside the trajectory. The Tier 1
smoke run does not have one, because `--system-dir` put the system files
elsewhere, so pointing the pass at the smoke trajectory with the production
topology was the way to exercise it. On 1 ns:

| quantity | measured | window |
|---|---|---|
| lipids, per leaflet | 476, 238 / 238 | balanced |
| area per lipid | 68.9 A^2, sd 0.23 | 60 to 68 expected, 55 to 75 acceptable |
| thickness | 37.62 A, sd 0.11 | 36 to 41 expected |
| core waters | 10.1 | was 8,622 before the periodic fix |
| G6 | closed | the bilayer does not invalidate the result |

Area per lipid is 0.9 A^2 above the literature window and inside the gate's, so
the chunk reads `outside expected` while G6 stays closed. I first read that pair
as a contradiction. It is not: there are deliberately two windows, a tight one
from the literature and a wider one for the gate, and a value between them is
exactly the case they exist to distinguish. Worth recording because the mistake
was mine and the design was right.

The production run directory does have its `solvated.pdb`, so this pass will
find its topology when `b2ar_apo` lands rather than silently returning nothing.
That was the specific failure I went looking for and it is not there.

**A loop that could never stop, found by reading its own guard.** The streaming
analysis re-analyses and re-ingests every row every half hour, and then runs one
final pass with the bilayer QC once production ends. Its condition was
`while pgrep -f "tools/run_md.py"`. That can never go false. The script's own
text contains the pattern, so the shell that launched it carried
`tools/run_md.py` in its command line and matched: the loop was kept alive by
its own parent, and the condition was true because the loop was running.

Checked rather than reasoned about, by listing every process the guard matched:
three, of which one was PID 78508, the loop's own grandparent.

The symptom would have been the worst available. The streaming passes work, so
the explorer would have filled in correctly and nothing would have looked
wrong. Only the pass *after* the loop, the one that computes the bilayer QC and
G6 on the finished trajectory, would never have run, and the `membrane` table
would have stayed empty with no error anywhere to say why.

Replaced with two conditions on state rather than on a pattern that matches us:
the queue's own `queue done` marker, and a process check that excludes the
script's ancestors by PID. Both paths tested before the swap, including the case
where the marker is present and the case where no log exists at all.

**What the fix immediately produced.** Restarting the loop populated the
receptor half of the database for the first time: 2 systems rather than 1, and
136 reference structures rather than 20. Everything above about the receptor
ingest path having never run was fixed by the act of running it.

**P6 is no longer a control that exists only in fixtures.** Both P5 and P6 were
implemented with tests for the passing, failing and empty cases, and neither had
ever seen a real trajectory. Feeding the smoke run's bilayer through ingest into
a scratch database and running the controls command gives P6 a real verdict:
`1 chunks, APL 68.9 A^2 (68.9 to 68.9), 37.6 A  PASS`.

Worth being exact about why that passes when the chunk itself reads `outside
expected`. Section 5.8 says "expect roughly 60 to 68 A^2 ... outside 55 to 75
A^2 is G6", so the advisory window and the gate are two different things in the
spec and the code encodes both. P6 is the gate, and 68.9 is inside it.

That leaves **P5 as the only control still genuinely waiting**, and not for want
of code: it needs an inactive-start run, which means `b2ar_carazolol` or
`abl1_imatinib`, third and fourth in the queue. P4 is 4 of 4 on the start
structures, and the controls command reports P5 as `pending` rather than
omitting it, because a control nobody mentions reads as a control that passed.

**And a false alarm I raised myself.** Running `pytest -q` from the repository
root reported 17 failures. They are all in `vendor/boltzmaker/tests`, which
shells out to `vendor/boltzmaker/.venv/bin/python3`, an interpreter that does not
exist inside this pixi environment. `pixi run test` was already scoped to
`tests/` and reports 871 passing; the bare command was not. `pytest.ini` now
scopes collection so both agree.

Two process points from that, both about the same habit. The run I did before it
piped pytest into `tail`, so the exit code I read was `tail`'s and it was 0 while
pytest was returning 1: a suite I reported as passing had 17 failures in it.
And the count in `GATE_OPEN.md` has now been wrong twice, at 833 and 837,
because I wrote it from memory and then from arithmetic instead of from a run.
It is 871, measured.

## 4i. A resumed run lost its barostat, and nothing said so

This is the most consequential thing the readiness audit found, and it was
about to matter: `abl1_apo` is queued to resume from 4.0 ns.

`run_md.py` builds its System, creates the Context, and then branches. If a
checkpoint exists it loads it; otherwise it equilibrates. The barostat was added
inside the equilibration branches, which is the branch a resume skips, so a
resumed run never got one.

What makes it a silent fault rather than a loud one is how OpenMM behaves.
Checked rather than assumed, by loading the real `abl1_apo` checkpoint into a
freshly built Context:

    context built with 5 forces, OpenCL
    loadCheckpoint SUCCEEDED, time = 4.0000 ns
    params: {'MonteCarloPressure': 1.0, 'MonteCarloTemperature': 300.0}

The load succeeds, and the restored parameters include the barostat's own, from
a System that has no barostat in it. OpenMM restores context parameters whether
or not a force consumes them. So a resumed run would report `MonteCarloPressure`
while producing at constant volume, and the only symptom would have been a
density that stopped responding, in a column nobody reads row by row.

For `abl1_apo` that would have been 46 of its 50 ns as NVT. For a membrane row
it is worse, because the bilayer could not then relax its area and the area per
lipid that P6 checks would be whatever the box was fixed at.

Three things were wrong and all three are fixed:

1. **The barostat is now added on the resume path**, before the Context is
   built, because a Context cannot gain a force afterwards without being
   rebuilt. Both equilibration branches and the resume now call one `barostat()`
   helper: they had each built their own, which is how they drifted apart, and
   the comment above the Simulation already claimed "the barostat goes on
   before the Simulation is built" while no code did it.
2. **An equilibration checkpoint is no longer mistaken for a production one.**
   The membrane ladder saves after every rung, so a run killed part way down it
   left a checkpoint carrying a non-zero clock and a non-zero `k_restraint`.
   Resuming that as production would have started with the protein still held,
   counted the equilibration time against the 50 ns target and written a short
   row. `checkpoint_phase.json` now records which phase each checkpoint belongs
   to, absent meaning production so older checkpoints still resume, and an
   equilibration one makes the run repeat the equilibration rather than guess.
3. **The water branch saved its checkpoint before resetting the clock**, so a
   kill in that window left a checkpoint holding the equilibration time as
   though it were production progress. The reset now comes first.

The resume also sets `k_restraint` and `k_lipid` to zero explicitly. The last
ladder rung is already zero so a finished equilibration needs it, but whether
production is unrestrained should not depend on the ladder happening to end
there.

Verified after the change on the real checkpoint: six forces including
`MonteCarloBarostat`, loaded at 4.0000 ns, barostat present in the resumed
context. The running `b2ar_apo` is untouched, because its process holds the old
code in memory and the next row is the one that picks this up.

**The general shape, since this is the third time this session.** All three
faults were a precondition quietly not met: a loop guard that matched its own
parent, a bilayer pass that could not find its topology, and a branch that
skipped the force it needed. None of them raised, none of them logged, and each
produced output that looked entirely normal. The ones left to find are the ones
that still look normal.

## 4j. The imatinib row would have failed, and ABL1 is dephosphorylated

Two findings from checking the three rows that have not run yet, rather than the
one that is running.

### abl1_imatinib had no ligand and would have crashed at build time

`queue_md.py` had no `ligand` field at all. `abl1_imatinib_start.pdb` contains
STI, imatinib, 37 heavy atoms, and amber14 has no template for a drug, so the
row would have failed the moment it tried to build, after the three rows ahead
of it had spent their GPU hours. It is fourth in the queue, so this was about
twenty hours from happening.

Proved rather than assumed: after PDBFixer, STI is the only unmatched residue in
the topology, at 68 atoms once hydrogens are on.

**Passing the CCD definition is not enough**, which is worth recording because
it was my first attempt and it fails in a way that reads as something else:

    No template found for residue 262 (STI). The set of atoms is similar to
    CHL1, but is missing 15 H atoms and has 2 C atoms and 7 N atoms too many.

CHL1 is cholesterol. The SDF and the topology both hold 68 atoms, but the
generator matches on a chemical graph and a PDB-derived ligand has its bonds
inferred from distances, so the graphs differ. `tools/build_membrane.py` had
already worked this out for carazolol: take the chemistry from the CCD, the
coordinates from the deposition, add the hydrogens, and build topology and
template from one molecule. `run_md.py` now imports that same function rather
than keeping a second copy, takes the ligand out before PDBFixer adds hydrogens
and merges it back from the object the generator was given. The imatinib system
builds: 48,185 atoms, formal charge zero.

**And one guard that was wrong in an instructive way.** I added a refusal using
`getUnmatchedResidues`, and it rejected a setup that worked, because that method
reports against the static templates only and does not consult a registered
template generator. It now runs only when nothing is generating templates, which
is the actual hazard: a ligand in the start model and no `--ligand` given.

A test now reads each row's start model and asserts that any het residue is
either declared as that row's ligand or is one PDBFixer rewrites, taking the
second list from PDBFixer's own table rather than hard-coding it.

### ABL1's activation-loop phosphotyrosine is silently replaced

That test immediately found a second undeclared residue, and this one is not a
packaging detail. `abl1_apo_start.pdb` carries **PTR 393**, phosphotyrosine, and
PDBFixer reports and performs `PTR 393 -> TYR`. The built system contains no
PTR, and deposited residue 393 reads `TYR` in it.

Tyr393 is the ABL1 activation-loop tyrosine, and its phosphorylation is one of
the main things that holds the kinase active. So `abl1_apo`, already 4 ns in, is
a simulation of **dephosphorylated** ABL1 started from a phosphorylated crystal
structure.

This is a defensible choice and it is the usual one: amber14 ships no
phosphotyrosine, and simulating the dephosphorylated kinase domain is what most
ABL1 work does. It is not defensible to do it without saying so, and nothing
anywhere recorded it. Both builders now log every replacement, so the line
`replacing nonstandard residues: PTR 393 -> TYR` appears in the build log.

**What it means for the results**, stated plainly rather than left for a reader
to work out:

- The ABL1 rows measure a dephosphorylated kinase domain. Any statement about
  the activation loop or the DFG state holding is a statement about that, not
  about the phosphorylated protein the structure came from.
- P4 is unaffected: it scores the start structure's geometry at t = 0, and the
  geometry is the crystal's whatever the residue is called.
- It gives the apo DFG-in row a reason to drift towards inactive that is
  physical rather than an artefact, since the phosphate that would hold the
  loop is absent. If that drift appears, this is the first thing to say about
  it, and the honest reading is a limitation rather than a result.
- The spec does not mention phosphorylation, so this is not a departure from
  it. It is a choice the spec left open and the build made without noticing.

Not re-run. Adding phosphotyrosine means a parameter set amber14 does not ship,
and the four rows are already most of the compute budget. Recorded as a
limitation, which is what section 0.1 asks for when a choice is forced.

## 4k. Measured at the size a finished row will actually be

Every performance number in section 4 had been taken on a partial run: 43
frames of the kinase and 26 of the receptor. A finished 50 ns row at 100 ps is
500 frames, so none of those numbers answered the question that matters. Tiled
the receptor trajectory to 500 frames, which repeats coordinates and so measures
cost rather than science, and ran the real analysis and the real endpoints
against it.

**The analysis**, on 194,032 atoms and 500 frames: 19.8 s wall, DSSP over 322
residues, and a viewer bundle of 50 frames at 1.4 MB against section 6.3's
150 MB ceiling. The analysis JSON comes to 1.9 MB.

**The endpoints**, against a database holding the full 5,000 motif rows:

| route | time | payload |
|---|---|---|
| `/refstates/run/b2ar_apo` | 33.7 ms | 27.9 KB |
| `.../frames.json` | 17.1 ms | 97.5 KB |
| `.../dssp.json` | 7.5 ms | 161.5 KB |
| `/refstates/map/b2ar.json?run=` | 8.2 ms | 31.4 KB |
| `/refstates/references/b2ar` | 1.3 ms | 39.0 KB |
| `/refstates/` | 6.5 ms | 5.6 KB |

Against a two second budget, so there is nothing to do here. A run page plus its
three fetches is about 290 KB.

**And one measurement I nearly reported without having made.** The first pass
timed `frames.json` at 12 KB and 8 ms and I almost wrote that down as the
500-frame figure. It reads the SQLite table, not the analysis file, and the
database still held the real run's ten frames, so the number was the old one
measured twice. Ingesting the tiled analysis first gives 97.5 KB and 17.1 ms.
Pointing a test at the big input is not the same as the code under test reading
it, and the only thing that caught it was the payload looking too small.

## 4l. G-REUSE annotated every kinase against ABL1

Section 1 lists six goals and G-REUSE is one of them: the module has to work for
any kinase or class A GPCR a user uploads, not only these two examples. Section
6.5 makes it concrete as `flexappeal refstates annotate` with a topology, a
trajectory, a family and a protein.

The entry point existed, took all four, and was never once run on a third
protein. Run on CDK2, with a two-frame trajectory made from the 1HCK crystal so
no MD was needed:

    family detected as kinase
    KLIFS kinase 392, 43 structures
    final frame: whole score +0.001 (intermediate), 0 conflicts

Kinase 392 with 43 structures is **ABL1**. `--protein CDK2` was accepted by the
parser and then dropped: `annotate` called `kinase_setup(start, reference_pdb)`
and never passed it, so the setup took its own default. CDK2 was scored against
ABL1's active and inactive ensembles, and the answer came back as a clean
`intermediate` with no conflicts, which is exactly what it would look like if it
had worked.

The receptor half had the same hole, `gpcr_setup` defaulting to `adrb2_human`,
so any uploaded class A receptor would have been scored against beta2AR.

**Both halves now thread the protein through**, with the default chosen per
family rather than shared, because "ABL1" is not a receptor and one default
cannot be right for both.

**And the fix exposed the second half of the same bug.** With the lookup
corrected, CDK2 then failed on `--reference-pdb`, which defaulted to `2GQG`, an
ABL1 entry, and raised a bare `KeyError: '2GQG'` from a dictionary of CDK2's
structures. A user would have had no way to read that as "the default belongs to
a different protein". It now picks one of the kinase's own structures when none
is given and says which, and refuses with a message naming valid codes when the
one given is not that kinase's.

With both fixed, CDK2 resolves properly:

    KLIFS CDK2 is kinase 198, 457 structures
    no --reference-pdb given; resolving the pocket from 1AQ1
    active: 337 references
    inactive: 2 references

So the hook works on a protein this build has never touched. Two things follow
from the numbers, and both needed fixing rather than noting.

**337 downloads where a dozen would do.** The first corrected run took over
fourteen minutes and was still fetching when I stopped it, because CDK2 has 337
structures KLIFS calls active and the ensembles were quality-filtered but
uncapped. A state score is the distance to the *nearest* reference, so the three
hundredth structure costs a download and changes almost nothing. `annotate` now
caps at 40 per state, best first by the catalogue's own quality for KLIFS and by
resolution for GPCRdb, and the same run takes 90 seconds. The cap is the reuse
path's alone: the production path passes no limit, so every committed result is
still computed over all of ABL1's 12 and 45 and beta2AR's 72 and 28, and P4
still reads +0.305, -0.368, +0.593 and -0.934 exactly as before.

**And the ensemble that survived was one structure.** Of the 40 active
references kept, 39 were dropped at the core alignment, which needs a matching
core length and does not get one across a structurally heterogeneous set:

    active: 40 references, the best of 337 by quality
    inactive: 2 references
    aligned 1 references onto 1F5Q: core residual mean 0.0 A, max 0.0 A, 28 dropped
    the active ensemble has only 1 reference structure, so every active
    distance here rests on 1; treat the scores as indicative rather than measured

That last line is new. How many structures a score rests on is part of the
score, and without it a comparison against one reference reads exactly like a
comparison against forty-five. It is written into the analysis payload as well
as the log, so anything downstream can see it.

The honest summary of G-REUSE after this: the plumbing works on an unseen
kinase, resolves its own catalogue entry, caps its own downloads and says when
its own ensemble is too thin to trust. For CDK2 specifically the answer is not
yet worth much, because the core-based alignment keeps one reference in forty.
Improving that means a core chosen per target rather than from one reference
structure, which is a larger change than this pass, and it is now a stated
limitation with a number attached rather than a silent one.

Every existing caller passes these positionally and still gets the right
default, so the four production rows are unaffected.

## 4m. The kinase scores rest on five active structures, not twelve

Chasing why CDK2 kept one reference in forty led back to the worked example, and
to a number I had been quoting wrongly all build.

Every log line and several paragraphs of mine say ABL1 has 12 active and 45
inactive references. Those are **catalogue rows**. KLIFS lists one row per chain
and per alternate conformation, and the motif ensemble is keyed on the
deposition, so what the scores actually compare against is:

| state | catalogue rows | structures | usable after motif resolution |
|---|---|---|---|
| active | 12 | 5 | 5 |
| inactive | 45 | 14 | 14 |

The five are 2F4J, 2GQG, 2HZ4, 2V7A and 4TWP. The fourteen inactive are
2E2B through 6NPV. Eighteen of the inactive rows are additionally recorded as
partial for a missing activation loop, which is a real limitation of those
structures rather than a fault here.

Nothing is lost on the way: the ensemble builder records a reason for every row
it does not keep, and it recorded none, because the rows collapse onto
depositions rather than being dropped. D016 already made exactly this choice for
the manifest and for the same reason, that chains of one crystal are not
independent. What was wrong was only the reporting, which counted the flattering
number.

**Why it matters enough to correct.** An active state score is the distance to
the nearest of five crystal structures. That is a legitimate measurement and it
is what the reference-ensemble method is, but it is a much smaller basis than
"12 references" suggests, and a reader deciding how much to trust a motif score
needs the real figure. The receptor is in a different position: GPCRdb lists one
row per structure, so beta2AR's 72 active and 28 inactive are 72 and 28.

`analyse_run.py` now logs both: `active: 12 catalogue rows over 5 structures`.

## 4n. The membrane panel had never been rendered with a membrane in it

Every streaming pass through the night printed `membrane 0 rows`, fourteen times
in a row, on a run that is a POPC bilayer. That turned out not to be a fault:
`tools/stream_analysis.sh` passes `--skip-membrane`, because the bilayer pass
needs the full solvated trajectory while the state analysis needs protein and
ligand only, and the one is cheap to run every half hour where the other is not.
The final pass, after production stops, runs it.

The fault was that the final pass had never been run on this trajectory, so
nothing about the bilayer panel had ever been exercised against real data. Three
things were wrong with it, and all three would have appeared for the first time
on the finished row.

**The tilt column was a dash, and always would have been.** `analyse_chunk`
takes `tm_ca_indices` as an optional argument and measures the receptor tilt
only if it is given. Nothing gave it. `tools/ingest.py` calls
`membrane.analyse_chunk(piece, start_ns=..., end_ns=...)` and always had, so
every chunk this build could ever have written had `tilt_deg = None`, and
section 5.8 asks for receptor tilt by name. An optional metric that the caller
has to remember to switch on is a metric that does not get measured, so it is
now derived inside `analyse_chunk` when the caller names nothing, and the
argument survives only as an override.

Deriving it needs the transmembrane bundle rather than the whole receptor, and
the bilayer locates that better than a residue table does: a Calpha between the
two phosphate planes is membrane-spanning by construction, so the set follows
the system that was built rather than a published segment range that may not
match the construct. On `b2ar_apo` that selection finds 201 of the 322 protein
Calpha, and they fall into exactly seven contiguous runs of 27 to 30 residues.
That is the 7TM bundle, recovered from the geometry with no catalogue lookup.

**Then the measurement itself was wrong.** With the indices finally passed in,
the first number out was 25.1 degrees of tilt, which is far too much for a
receptor that was built square in the membrane. `receptor_tilt` took the
principal axis of all the membrane-spanning Calpha at once, and its docstring
asserted that for a 7TM receptor this is the bundle axis. It is not. Measured on
this system the slab Calpha span 40 by 45 by 37 angstrom with singular values
159, 139 and 109: a first-to-second ratio of 1.15. A cloud that round has no
principal direction, so the axis returned is set by noise, and the tilt it
reports is noise dressed as a measurement.

The spec had already written down the right definition. Section 4.3's fallback
orientation uses "the mean direction of the seven TM helix axes", and each helix
on its own is a rod whose axis is well determined: the same structure gives
first-to-second ratios of 4.8 to 7.8 per helix. Averaging those seven axes, each
folded towards +z first because consecutive TM helices run in opposite
directions and would otherwise cancel, gives 6.7 degrees for the as-built system
and 13.5 degrees averaged over the first 7.6 ns of production. A receptor that
settles from 7 to 13 degrees during equilibration is what equilibration does.

The broken function was correct for one helix, which is what both of its tests
fed it, and that is why they passed. The error was in applying it to a bundle.
So it is now two functions with different arguments rather than one: `helix_axis`
and `axis_tilt_deg` for a single rod, `bundle_tilt` for a receptor. Nothing can
reach for the wrong one by accident, because they do not take the same inputs.
The regression test builds a seven-rod bundle with a tilt known by construction
and asserts the two apart: the per-helix mean recovers the angle to within half
a degree, and the test also asserts that the whole-cloud axis is more than 20
degrees out, so if the fixture ever stops demonstrating the fault it was built
for, it says so rather than passing quietly.

**And the panel showed numbers with no verdict.** Section 5.8 asks for the panel
"so a user can see the bilayer behaved", and five numbers against two windows
held in the reader's head is not that. The verdicts exist: `MembraneChunk` has
`apl_verdict` and `thickness_verdict`, and `as_row` emits both. `put_membrane`
drops them, which is correct, because section 7 gives the `membrane` table six
columns and no verdict. A verdict is a pure function of the number and the
window, so the panel recomputes it through `MembraneChunk` itself rather than
the schema growing a column the spec does not have, and the windows keep one
definition instead of a second copy in the view.

Three states, not two. An unbadged row cannot be told apart from one with
nothing measured, so "expected" carries a quiet badge too, while a filled badge
stays reserved for the exception. The first real chunk came out at 68.3
angstrom squared per lipid, which is 0.3 above the narrow window and well inside
the acceptable one, so the panel's first live reading is an amber "outside
expected" next to an "expected" thickness of 38.3 angstrom. That is the panel
doing its job on its first real data.

## 5. Production runs

Tier 2 started at 17:04 on 10 October, in section 5.3's priority order: β2AR
active-start first, then the ABL1 rows, then β2AR inactive. Four rows at 50 ns,
about 59 GPU-hours. This section fills in as they land, through the streaming
analysis in `tools/stream_analysis.sh`.

### 5.0 The figures

Regenerated by `tools/make_figures.py` from the explorer database, so they are
reproducible rather than drawn once by hand. They are committed because a figure
that has to be regenerated to be looked at is a figure nobody looks at.

![ABL1 apo per-motif state score](data/figures/abl1_apo_motifs.png)

![beta2AR apo per-motif state score](data/figures/b2ar_apo_motifs.png)

![beta2AR reference state map](data/figures/b2ar_map.png)

The receptor map is the one to read carefully. The active references cluster at
a TM3-TM6 distance of 15.3 to 16.5 Å and the inactive ones at 8.0 to 8.5 Å, and
`b2ar_apo` sits in neither: it runs between 11 and 14 Å for the whole of
production so far, having left the active cluster during equilibration. The
cleft threshold is 11.39 Å, so the frames still classify as active, but with
little margin.

That is the apo system doing what an apo system should. 4LDE was solved with an
agonist and a nanobody, and this one has neither, so the cytoplasmic cleft
closes. The NPxxY axis is the useful counterpart: it stays at 2.6 to 3.8 Å from
the nearest inactive reference, which is active-like, so the receptor has not
moved towards the inactive state so much as relaxed out of an agonist-bound
one. Two motifs disagreeing in this way is more informative than a single
whole-protein number, and is the reason section 2 kept the scores per motif.

**And it is a working argument for not classifying on a threshold.** Over the
first 8.5 ns the TM3-TM6 distance runs from **11.17 to 15.32 Å**, mean 13.02,
and one frame is already below the 11.39 Å cleft threshold that section 1.2b
fitted. A classifier that read the state off that one distance would therefore
call a frame of this run inactive, eight nanoseconds into a receptor that is
active-like on every other measure: NPxxY sits at +0.63, TM6 at +0.49, TM7 at
+0.48. The ensemble score does not, because it asks how far the motif is from
the nearest active and the nearest inactive structure rather than which side of
a line it falls. The threshold earns its place classifying crystal structures,
where the two clusters are 7 Å apart and nothing sits between them. A
trajectory does sit between them, which is the whole reason it is interesting,
and that is where a single cut stops being informative.

### 5.1 Measured speeds

| row | atoms | ns/day |
|---|---|---|
| ABL1 in water | 46,721 | 213 |
| β2AR in POPC | 194,032 | 59.8 |

The membrane figure is with the `MonteCarloMembraneBarostat` running. Without
it, during the NVT heating stages, the same system does about 86 ns/day, so the
barostat costs roughly 30 per cent. That difference is the reason the first
compute projection in `GATE_OPEN.md` was too optimistic and was corrected.

### 5.2 The bilayer behaves

From the Tier 1 run, after 0.2 ns of equilibration:

| measure | value | expected |
|---|---|---|
| area per lipid | 68.4 Å², from 71.6 as built | 60 to 68 |
| P–P thickness | 37.6 Å | 36 to 41 |
| leaflet balance | 238 / 238 | even |
| waters in the hydrophobic core | 4 | a handful |
| temperature | 303.5 to 304.3 K | 303.15 |
| density | 1.017 g/mL | about 1.02 |

Area per lipid is 0.4 Å² outside the expected window and still falling from
where the builder left it, which is what an equilibrating bilayer should do:
`addMembrane` removes lipids near the protein and leaves the rest with more area
each. G6 stays closed. The full per-chunk series over the 50 ns rows will say
whether it settles inside the window.

### 5.3 The 4 fs step holds, so the 2 fs fallback is not needed

Section 5.2 sets a 4 fs step with hydrogen mass repartitioning and names a
fallback: 2 fs without HMR if the energy drift exceeds 0.01 kT/ns/dof.
Measured by `tools/check_drift.py`:

| run | span | mean T | setpoint | T slope | energy drift |
|---|---|---|---|---|---|
| abl1_apo | 4.2 ns | 300.77 K | 300 | −0.36 K/ns | within noise |
| β2AR Tier 1 | 0.9 ns | 303.68 K | 303.15 | −0.41 K/ns | within noise |

Both hold their setpoint and neither energy slope is resolvable from its own
uncertainty, so the fallback stays unused.

**The number is reported carefully, because it is easy to over-read.** These
runs use a Langevin middle integrator and, in production, a barostat; both
exchange energy with a bath by design, so the total energy is not a conserved
quantity and a non-zero slope is not a fault on its own. What the threshold
usefully catches in a thermostatted run is a *systematic* trend alongside a
temperature that has left its setpoint. The tool therefore reports the slope,
its uncertainty and the temperature together, and prints "within noise" rather
than quoting a figure the data cannot resolve. Quoting 0.0003 kT/ns/dof from ten
frames would be a false precision dressed as a pass.

"Resolved" means three standard errors, not one, and the difference matters. A
one-sigma criterion calls about a third of pure-noise draws a real drift, so a
gate built on it would flag a third of perfectly good runs and quickly be
ignored. Three sigma brings that to roughly 0.3 per cent. This was found by a
test rather than by reasoning: the first seed tried gave a 1.7 sigma
fluctuation on noise and read as a genuine drift, and the test now measures the
false-positive rate over sixty draws instead of trusting one.

The degrees of freedom come from the serialised system rather than from 3N:
HBonds and rigid water remove 167,040 constraints from the membrane system's
194,032 particles, so using 3N would understate the per-degree-of-freedom drift
threefold.

### 5.4 Section 6.3's performance targets, measured

| target | budget | measured |
|---|---|---|
| explorer first paint | < 2 s | **176 to 456 ms**, median 396 |
| frame scrub | < 100 ms | within budget on 20 relayouts |
| per-system web bundle | ≤ 150 MB | **0.4 MB** for the partial run, about 1 MB for a finished row |

First paint is measured from the browser's own `first-contentful-paint` entry
rather than by timing the navigation call, because the target is about when a
reader sees something and not when the last asset arrives. The viewer and the
plots load after first paint by design, so a slow Mol\* does not count against
it.

The bundle is far inside its budget because the viewer coordinates are protein
and ligand only at a 1 ns stride. The solvated trajectory a membrane row
produces is about a gigabyte, and none of it is drawn.

The paint test was intermittently failing when first written, and the cause was
the test rather than the page: it read the paint entry immediately after
`DOMContentLoaded`, which on a page that paints in 400 ms is sometimes before
the entry exists. It now waits for the entry, which does not inflate the
measurement because the entry carries its own `startTime`.

### 5.5 One preliminary observation, on 1 ns

Worth stating because it is the first real dynamics, and worth hedging for the
same reason. The Tier 1 β2AR run is 1 ns of the active-start system, and over it
the whole-protein score falls from **+0.593** on the static start structure to
**+0.277**:

| motif | s at 0.9 ns | nearest active | nearest inactive |
|---|---|---|---|
| NPxxY | +0.589 | 3SN6, 0.66 Å | 6PS0, 2.57 Å |
| TM6 | +0.357 | 8GE8, 1.88 Å | 3PDS, 3.96 Å |
| TM7 | +0.303 | 7DHI, 1.07 Å | 6PS4, 2.00 Å |
| H8 | +0.188 | 8GED, 0.99 Å | 6PS4, 1.44 Å |
| TM5 | +0.171 | 9BUY, 1.47 Å | 3PDS, 2.08 Å |

The direction is what an apo receptor built from an agonist-bound active
structure should do: 4LDE was solved with its agonist, this system has none, and
the motifs move off the crystal geometry. It is still clearly active-like.

**What this is not.** One nanosecond, one replica, and no equilibration between
the start structure and the first frame of the comparison. Most of a 0.32 drop
on that timescale is relaxation away from a crystal pose rather than progress
towards the inactive state, and the inactive distances barely move. Nothing here
distinguishes the two, and on the Dror timescale nothing at 50 ns will either.
The 50 ns rows will show whether the drift continues or settles, which is a
different and answerable question.

## 4o. The deploy would have shipped a gigabyte of trajectory

G3 has not been run, so this was found by reading rather than by failing.
`deploy/deploy.sh` rsyncs the project root with excludes for the virtual
environments, the scratch directories and the secrets, and no exclude for
`data/`. That was harmless while `data/` was manifests and a small database. It
is not harmless now: `data/md/` is 565 MB with one row a sixth of the way
through its 50 ns, and four finished rows are several gigabytes of trajectory,
checkpoint and system XML. `data/cache/` and `data/refs/` add another 92 MB of
structures downloaded so a rebuild does not hit KLIFS and GPCRdb again. The web
app reads none of the three.

Measured on a dry run, the deploy sends **1.08 GB today and 193 MB with the
excludes**, and the difference grows with every hour of production.

The tidier fix would have been `--filter=':- .gitignore'`, which would make the
deploy send what the repository considers source. It would also have been
wrong, and quietly so. `data/refstates.db` is gitignored, and it is exactly
what the explorer reads. The app treats an absent database as an empty section
rather than an error, deliberately, because the feature is additive and a
deployment without it should not 500. So the filter would have produced a
droplet that came up perfectly, with the reference-states tab showing nothing,
and nothing anywhere saying why. The explicit list is less elegant and says
what it means.

A test asserts both halves: that the three working directories are excluded,
and that the three runtime paths are not. Verified against a dry run rather
than by reading the flags, after a first check that found every file missing
because `rsync -n` without `-v` prints no file list at all, so the list being
grepped was empty. That is the same shape as reading `tail`'s exit code instead
of pytest's, which section 4h records: the check ran, the check said nothing,
and nothing is easy to read as a pass.


## 4p. Half the run matrix had an analysis path nobody had ever run

`abl1_apo` and `b2ar_apo` have been analysed every half hour for a day.
`abl1_imatinib` and `b2ar_carazolol` are built, queued, and have never produced
a frame, so nothing downstream of the simulation had ever been run on either.
They are not the same shape as the two that have: one is a kinase with a ligand
in the binding site, the other is the only inactive-start receptor, and between
them they are the only rows that can satisfy P5. Finding out whether they
analyse by waiting for the queue means finding out after about 24 GPU-hours.

So they were rehearsed instead. A twelve-frame trajectory was written from each
built system's own coordinates with a little noise on each frame, in the
scratch directory rather than in `data/md/`, where the streaming loop would
have ingested it into the live database. That tests the plumbing and not the
physics, which is the question being asked: can these rows be analysed at all.

Both can, first time, with no conflicts:

| row | whole score | motifs |
|---|---|---|
| `abl1_imatinib` | **−0.281**, inactive | xDFG −0.867, A-loop −0.379, all six resolved |
| `b2ar_carazolol` | **−0.662**, inactive | all five inactive, TM6 −0.821, NPxxY −0.758 |

Those are the right answers. `abl1_imatinib` is imatinib-bound and DFG-out, and
its xDFG motif comes out 4.27 Å from the nearest active structure and 0.31 Å
from the nearest inactive one. `b2ar_carazolol` starts from 5D5A and sits on top
of it, 0.35 to 0.41 Å away on four of five motifs.

**And P5 produced output for the first time.** It has read `pending: no
inactive-start run analysed yet` for the whole build, because the only two rows
it applies to are last in the queue. Against the rehearsal database:

    abl1_imatinib     12 frames, mean s -0.295, 100% below zero  PASS
    b2ar_carazolol    12 frames, mean s -0.674, 100% below zero  PASS

The numbers are from a fabricated trajectory and mean nothing about the
science. What they establish is that the control reports rather than staying
silent, which is the thing that could not be known until a row with
`start_state` of `DFG-out` or `inactive` reached the database.

P6 on the receptor gives an area per lipid of 71.6 Å², outside the expected 60
to 68 and inside the acceptable 55 to 75, so it passes and is flagged. That is
a useful cross-check rather than a result: 71.6 is exactly the as-built figure
section 5.2 records before equilibration pulls it to 68.4, and the rehearsal
trajectory is the built system. The measurement agreeing with a number taken
independently a day earlier is evidence the fabricated frames are being read
correctly.


## 4q. A row that died would have shown as running until somebody looked

Section 0's fallback ladder ends with degrading gracefully and recording
`status = "failed:<reason>"`. For references that happens: `refs.py` writes
`failed:no_motifs`, `failed:missing_motifs` and `failed:<ExceptionName>`. For
runs, nothing wrote it at all.

`run_md.py` writes `run.json` when it finishes and not otherwise.
`tools/ingest.py` reads an absent `run.json` as a row still in flight, which is
right while a row is in flight and wrong once it has stopped:

```python
else:
    # A run still in flight is recorded as such rather than left out
    fields["status"] = "running"
```

So a row that died at three in the morning would be ingested as `running`, and
the explorer would show it as `running` for as long as anyone cared to look.
The partial trajectory beside it would make that look entirely plausible: the
streaming analysis would keep re-reading the same frames and reporting the same
motif counts every half hour, because mdtraj reads whatever is on disk and
nothing on disk says the writer has gone.

The record is now written by the queue rather than by the runner. The queue is
what survives the failure, which is the point: a runner killed by the OOM
reaper or by a segfault inside the platform layer gets no chance to write
anything, and those are the failures most worth recording. The reason is taken
from the last line of the row log that looks like an error, so a NaN
coordinate is recorded as `failed:Particle coordinate is NaN` rather than as
`failed:exit 1`. A row whose runner did write its own `run.json` keeps it,
because the runner's record is the better one.

Found by asking what the explorer would show if a row failed overnight, which
is a question with an obvious answer and no test behind it. There is now a test
for the recording and a second for the rendering, because a failed row also has
to render rather than 500 or disappear.


## 4r. Section 6.2's viewer asks for four things and had one

The run page's Mol* viewer is specified in one line of section 6.2: "cartoon
coloured by motif tokens; the GPCR shows the membrane as transparent phosphate
spheres (toggle); frame scrubber; reference overlay toggle (nearest active and
nearest inactive reference superposed on the core, ghosted)."

The frame scrubber is there. The cartoon colouring is not, and the page says so
in its own read-more, which is the right way to be missing something:

> Colouring the cartoon by motif is not wired up yet; the table is the honest
> version of it rather than a half-working theme.

The other two were missing with nothing anywhere saying so, which is the part
worth recording. A deliberate omission that is disclosed is a decision; an
undisclosed one is just absent.

**The membrane toggle is now built**, because the reason it was not turned out
not to survive checking. The viewer bundle drops the solvent and the bilayer on
load and explains why: a 190,000-atom membrane system is mostly water and none
of it is drawn. That is right for the water. It is not right for the
phosphorus. A POPC bilayer around β2AR is 476 lipids, so one atom each is 476
atoms against roughly 2,600 for the protein, and measured rather than
estimated the pair comes to **38 KB of PDB and 26 KB of XTC** beside a 0.6 MB
protein bundle. The gigabyte the comment is defending against is the water, and
the phosphorus was never part of it.

It is written as its own structure rather than merged into the viewer bundle.
Merging would put non-protein atoms into the topology every downstream
selection runs against, and that path is currently producing correct numbers
for a live run; a second structure the viewer loads beside the first cannot
affect any measurement. It is loaded on demand and deleted again rather than
loaded and hidden, which also avoids Mol*'s visibility API, the part of its
surface most likely to move between versions: `loadTrajectory` returns the
model's state ref, and deleting that subtree is enough.

Verified in a real browser rather than by reading the code. Both files serve
200, the load resolves, the delete commits, the control reports what it did,
and there are no console or page errors on either transition. What could not be
verified here is that the spheres are actually drawn: this Chromium has no
WebGL, so Mol* renders its "WebGL does not seem to be available" notice and a
screenshot proves nothing either way. The structure reaching the plugin's state
is as far as this environment can take it, and that is stated rather than
dressed up as a visual check.

**The reference overlay is still missing**, and is now written down as missing
rather than merely absent. It needs the nearest active and nearest inactive
reference superposed on the core and ghosted, which is more than a file to
serve: the superposition exists in the analysis, but the overlay needs those
coordinates written out per frame and a second representation style, and doing
it badly would put a structure on the screen that is not where the analysis
says it is. The honest position is the one the cartoon colouring already takes.


---

## 4s. Three columns of section 7's schema were never written

`systems` is specified in section 7 as

```sql
CREATE TABLE systems (system_id TEXT PRIMARY KEY, family TEXT, uniprot TEXT,
                      name TEXT, environment TEXT, n_atoms INTEGER,
                      forcefield TEXT, notes TEXT);
```

and the store creates all eight columns. `ingest_analysis` wrote four of them.
`uniprot`, `n_atoms` and `forcefield` were NULL for both systems and always had
been, because nothing passed them.

The consequence was visible on the page and read as a design choice. Section
6.2 asks the Examples cards for a thumbnail, an atom count, a simulated time
and a generated finding. The card template already had the atom count in it,
behind `{% if card.system.n_atoms %}`, so it rendered nothing and looked like a
card that was never meant to show one. The markup was right and the data was
absent, which is the harder of the two to notice.

All three are now written, and none of them is guessed:

| system | uniprot | n_atoms | forcefield |
|---|---|---|---|
| `abl1` | P00519 | 46,721 | not recorded for this build |
| `b2ar` | P07550 | 194,032 | amber14-all.xml, amber14/tip3p.xml |

The accessions come from the catalogues the build already queries, not from
memory: KLIFS `kinase_information?kinase_ID=392` returns `"uniprot": "P00519"`
and GPCRdb `protein/adrb2_human/` returns `"accession": "P07550"`. They are
pinned in `tools/ingest.py` rather than fetched on every ingest, because an
accession does not change and the ingest should not need the network, with the
endpoint named beside each so the value can be rechecked.

The atom counts are measured. For the receptor they come from `build.json`,
which `tools/build_membrane.py` writes; the water builder writes no manifest,
so ABL1's is counted from the solvated topology, where counting coordinate
records is exact even above 99,999 atoms because the record name does not
change. Both agree with section 5.1's independently recorded figures, 46,721
and 194,032, which is a cross-check rather than a coincidence.

ABL1's force field stays null, because this build has no record of it that is
not a guess. `run_md.py` now writes the list it actually loaded into
`run.json`, so the next finished water row fills the column by itself. A null
rather than a plausible string is the point: section 0 says to omit a key that
cannot be produced instead of writing a stand-in.

**And that change nearly cost a run.** The force field list is chosen inside
`build` and recorded inside `main`, which are different scopes, so
`"forcefield": forcefield_files` parsed cleanly and would have raised
`NameError` at the one moment it costs most: after the production loop had
finished and before anything was written down. A seventeen-hour row would have
completed and recorded nothing, and `record_failure` from section 4q would then
have marked a successful run as failed. It is a module constant now, and there
is a test that reads `main` with `ast` and asserts every name it uses is one it
binds, a module-level name, or a builtin. The test was checked against the bug
rather than only against the fix: reintroducing the local makes it fail.

### The card wanted a thumbnail too

Section 6.2 asks each Examples card for "a thumbnail, atom count, simulated
time, and a one-line finding generated from the data". The finding was there
and is the hard one. The atom count was markup waiting for data, above. The
thumbnail was missing outright.

There is no GPU in this pipeline, so a cartoon render is not available, and a
stock picture of a kinase would be the opposite of "filled from SQLite, never
written by hand". So the thumbnail is a line drawing of the coordinates that
were actually built: the Calpha trace in order, and for the receptor the 476
lipid phosphorus atoms as well, projected looking along the membrane. 32 and 36
KB. β2AR comes out immediately recognisable, seven helices crossing a bilayer
with the intracellular loops hanging below it.

The first version was wrong in a way worth recording, because it looked fine.
The receptor appeared to sit mostly above the upper leaflet, half out of the
membrane. These are wrapped coordinates, and under a membrane barostat the
bilayer drifts until it straddles the periodic boundary, so the two leaflets
had been drawn at opposite ends of the box. That is the same fault section 2
records `centre_bilayer` being written for, and `centre_bilayer` already
returns its shift so that a caller can move other atoms into the same frame.
The protein is now shifted with the lipids, and any Calpha still more than half
a box from the bilayer centre is moved the short way round, because that is
what periodicity means. A picture is as capable of being confidently wrong as a
number, and harder to check.

## 4t. Checking the About page against the nine things section 6.2 names

Section 6.2's About item asks for a "software and data table (OpenMM, PDBFixer,
MDTraj, MDAnalysis, gemmi, KLIFS, GPCRdb, OPM, BoltzMaker's `sse_comparison`
with commit SHA), references verified at build time". That is a checklist, so
it was checked rather than assumed. Seven of the nine were there, the
references do verify through Crossref at build time, and the vendored SHA
`da939b1be54a` is on the page as asked.

**PDBFixer was missing**, and it is not a small omission. PDBFixer is what
`tools/run_md.py`, `tools/build_membrane.py` and `tools/verify_references.py`
use to add missing atoms and replace nonstandard residues, and replacing
nonstandard residues is the whole of D018: it is the reason both ABL1 rows
simulate an unphosphorylated kinase. The About page already carries that
limitation and did not name the tool that causes it. PDBFixer has no paper of
its own, being part of the OpenMM project, so it is named on the row that
already cites OpenMM rather than repeating the same reference: "Simulation
engine, and structure preparation through PDBFixer".

**MDAnalysis is not there because it is not used.** It appears once in this
repository, in `flexappeal/runtime/pixi.toml.j2`, which is the environment for
the run bundle a user downloads and executes on their own machine, not for the
reference-states module. Nothing in `flexappeal/refstates/` or in `tools/`
imports it: the trajectory work is mdtraj and the structure reading is gemmi.
Section 6.2 names it in an expected list, and the honest response to a list
naming a library this build does not use is to not list it as though it were
used. Recorded here instead, which is the difference between a gap and a
decision.


## 4u. The motif table was missing the column my own decision relied on

Section 6.2 specifies the motif table as "per motif, start value, final 10 %
mean, nearest reference PDB at end, `s`, label". The page had five of the six.
The nearest reference was not there.

That matters more than a missing column usually would, because D023 had just
been written, declining to build the viewer's reference overlay on the grounds
that

> The motif table below gives the nearest reference and the distance to it, per
> motif, which is the same comparison in numbers.

It did not. The justification for not building one thing rested on another
thing that was also absent, and I wrote it without checking. A decision is only
as good as the claim it is resting on, and this one was resting on a column
that had never existed.

The data was there the whole time: every frame in the analysis payload carries
`nearest_active`, `nearest_inactive`, `d_act_a` and `d_inact_a` per motif.
`motif_frames` has no column for it, and should not: section 7 does not give it
one, and the nearest reference is a property of a frame rather than of a run.
So the view reads it from the payload it already opens for the viewer bundle
and merges it into the summary rows, which is the same shape as the membrane
verdicts in section 4n.

Both references are shown, not the nearer one. The score is a comparison of two
distances, and naming only the winner would hide the case the comparison exists
for. β2AR's H8 at the end of the run reads 0.63 Å from 8GE1 and 0.85 Å from
6PS0, which is what `intermediate` means; NPxxY reads 0.48 Å from 3P0G against
3.30 Å from 8W1V, which is what `active` means. One column, and the label stops
being something you have to take on trust.

### And the timelines are not the chart the spec describes

Section 6.2 asks for "one small panel per motif metric, shaded bands for the
active and inactive reference ranges (5th to 95th percentile of each
ensemble)". What is on the page is one panel for all motifs, plotting the state
score `s` against time, with bands at ±0.2 marking where the label changes.

That is a defensible chart and arguably the better one for this build, because
`s` is the single comparable quantity across motifs whose metrics are in
different units, and the ±0.2 bands are where the label actually changes. But
it is a different chart from the one specified, and until now nothing said so.
It is recorded here rather than rebuilt: the per-metric panels with ensemble
percentile bands would need the ensembles' per-metric distributions carried to
the page, which is a larger change than this pass, and the information a reader
would get from them is in the reference state map and the new nearest-reference
column.


---

## 6. The positive controls

Section 9 names six. Four are decided, one is pending on data that is still
being generated, and one is pending on the final analysis pass. Everything here
is `flexappeal refstates controls` output rather than a transcription, and the
command is the thing to re-run rather than this table to trust.

| control | what it asks | result |
|---|---|---|
| **P1** | Reference structures run through the trajectory code path reproduce their own KLIFS DFG labels | **17 of 17** |
| **P2** | The same, for the αC-helix labels | **17 of 17** |
| **P3** | The β2AR cleft threshold separates the catalogue's states, leave-one-out | **100 of 100**, needs 90% |
| **P4** | Each start structure scores as the state it was built from at t = 0 | **4 of 4** |
| **P5** | An inactive-start run stays inactive over its final 20 ns | **pending** |
| **P6** | Every membrane chunk inside section 5.8's windows | **pending** |

P4 in full, since it is the one that tests the whole scoring chain end to end
on real structures rather than on references the chain was fitted to:

| row | expected | s at t = 0 | label | verdict |
|---|---|---|---|---|
| `abl1_apo` | active | +0.305 | active | PASS |
| `abl1_imatinib` | inactive | −0.368 | inactive | PASS |
| `b2ar_apo` | active | +0.593 | active | PASS |
| `b2ar_carazolol` | inactive | −0.934 | inactive | PASS |

The band is ±0.3, and it is deliberately wide: a near-zero score on a structure
that was built from a known state is a failure of the scoring, not a marginal
pass, so the band has to exclude it.

**P5 is pending on the queue, not on the code.** It reads the explorer database
and reports on runs whose `start_state` is `DFG-out` or `inactive`, which are
`abl1_imatinib` and `b2ar_carazolol`. Both are last in section 5.3's priority
order, so neither has produced yet. The code path has tests against synthetic
runs, and the start states in `tools/ingest.py` match the two the control looks
for, so it will report rather than stay silent when they land.

**P6 is pending on the final analysis pass.** The streaming passes run with
`--skip-membrane`, so no membrane chunk has been written to the live database
and P6 has had nothing to read all build. Rehearsed against a scratch database
on 7.6 ns of `b2ar_apo`, the pass produces a chunk at 68.3 Å² per lipid and
38.3 Å thickness with G6 closed, which is section 4n.

---

## 7. Limitations

Stated here rather than left for a reader to infer. None of these is a fault to
be fixed later; each is a boundary on what the results can be asked to support.

**Timescale.** Each row is 50 ns. The β2AR active-to-inactive relaxation these
systems are built to illustrate occurred on the microsecond scale in Dror 2011,
roughly twenty times longer. Nothing at 50 ns can show that transition, and its
absence is an expected limitation rather than a negative result. What 50 ns can
show is whether a start state is stable, which is what P5 asks and what the
state scores answer. Section 9 says this in the same terms, and section 5's
β2AR drift is hedged on exactly this ground: most of an early fall in the score
is relaxation away from a crystal pose rather than progress towards the other
state, and nothing here distinguishes the two.

**One replica per row.** Section 9's P5 asks for 2 of 3 replicas. D009 reduced
that to one replica for compute, so the criterion became "the inactive-start run
keeps its score below zero". One replica cannot support a majority rule, and
reinterpreting a 2-of-3 rule as a 1-of-1 pass would be claiming a majority from
a sample of one. The adaptation is recorded rather than quietly absorbed, and no
result here should be read as reproducible across replicas, because none was
tested that way.

**POPC only, and no palmitoylation.** The β2AR bilayer is pure POPC. A real
plasma membrane has cholesterol, mixed chain lengths and charged headgroups, and
cholesterol in particular is known to sit in β2AR grooves and change the
receptor's dynamics. Palmitoylation at Cys341 anchors helix 8 to the inner
leaflet in the real receptor and is not modelled here, so helix 8's angle and
its score are measured on a receptor that is less constrained at that end than
the real one. The H8 numbers should be read with that in mind.

**Apo core runs.** Two of the four rows are apo. `abl1_apo` starts from 2GQG and
`b2ar_apo` from 4LDE, both of which were solved with something bound that these
systems do not have, so the first part of each run is the receptor relaxing away
from a ligand-bound pose. For β2AR this is visible: the TM3-TM6 distance leaves
the active cluster within the first nanoseconds. That is the apo system behaving
correctly rather than the method failing, but it means an apo row's early frames
are not a sample of the state its start structure was labelled with.

**ABL1 is dephosphorylated.** PDBFixer's `replaceNonstandardResidues` converts
the phosphotyrosine at 393 to plain tyrosine, so both ABL1 rows simulate the
unphosphorylated kinase. This is D018 and it is logged at build time now rather
than happening silently, but it remains a difference between these runs and the
crystal structures they are scored against.

**The kinase reference set is small.** ABL1's scores rest on 5 active and 14
inactive structures, not the 12 and 45 catalogue rows that KLIFS lists, because
KLIFS lists one row per chain and the ensemble is keyed on the deposition. A
kinase motif score is a distance to the nearest of five crystal structures.
β2AR is in a better position at 72 and 28, because GPCRdb lists one row per
structure. This is section 4m.

**The reuse path is capped, and one target shows it.** Annotating an unseen
protein caps each reference ensemble at 40 structures, best first. For CDK2 the
core-based alignment then keeps one reference in forty, so a CDK2 score is a
distance to a single structure and the analysis says so in its own output.
Improving it needs a superposition core chosen per target rather than taken from
one reference, which was larger than this pass. Section 4l has the numbers.
