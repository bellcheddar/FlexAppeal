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

### Re-measured after the explorer grew

Sections 4n to 4x added a bilayer toggle, two consensus columns and five motif
bands to the ribbon, a badge strip that updates on every slider move, and a
thumbnail to each card. Three of section 6.3's budgets touch that work, so they
were measured again rather than assumed to still hold. The numbers are in
section 5.4 beside the originals. Two things there are worth being precise
about rather than rounding into a pass.

The scrub figure is the synchronous handler cost, which is the part this work
added: setting the slider, mapping the viewer frame onto an analysis frame and
rewriting five badges. It is not the whole visible update, because Mol\* swaps
its model asynchronously and this does not wait for it. What can be said is
that the strip added nothing measurable to a budget of 100 ms.

The bundle projection is not linear, and the first version of it was. The
topology PDB is a fixed size whatever the run length, and only the coordinates
and the analysis JSON grow with frames, so projecting the whole directory
linearly from 4.3 ns overstated the kinase by most of its headroom. Split into
the parts that scale and the parts that do not, both systems land at about 3.5
MB. The largest growing piece is not the trajectory but the analysis JSON, at
545 KB for 13 ns.

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
its own ensemble is too thin to trust. For CDK2 specifically the answer was not
yet worth much, because the core-based alignment kept one reference in
twenty-nine.

### Resolved, and it was not what this section guessed

I wrote above that fixing it "means a core chosen per target rather than from
one reference structure, which is a larger change than this pass". That was a
guess, and measuring it found something smaller and worse.

The cores are not the problem. Across CDK2's twenty-nine references the KLIFS
scaffold gives 42 residues to twenty-eight of them and 41 to one. The problem
is which one `align_to_common_frame` picked as its target: the first in sorted
order. That was **1F5Q**, the odd one. Every other reference is dropped when
its core length differs from the target's, so the single shortest reference in
the set took the whole ensemble with it, and the survivor was the target
itself. Hence a "core residual mean 0.0 A, max 0.0 A" in the log, which is what
a set of one fitted onto itself reads as, and which should have been the tell.

The target is now the first, in sorted order, whose core is the **most common
length** in the set. On CDK2:

| | before | after |
|---|---|---|
| target | 1F5Q, 41 residues | 1FVV, 42 residues |
| aligned | **1** of 29 | **28** of 29 |
| core residual | 0.0 Å, which is one structure onto itself | mean 0.473 Å, max 1.425 Å |
| core agreement | not measurable | 0.49 Å mean, 1.425 Å worst |

The score on a two-frame 1HCK trajectory now comes out at **+0.364, active, 0
conflicts**, which is the right answer for an ATP-bound CDK2, resting on 28
structures instead of one. The thin-ensemble warning still fires, and now says
something true and narrower: CDK2's *inactive* side has one reference, because
the catalogue has two inactive entries for this kinase. That is a limit of
KLIFS and not of this code.

**Nothing in production moves.** Choosing by modal length cannot change a set
that is already uniform, and every ABL1 and β2AR ensemble here is: ABL1 aligns
19 of 19 onto 2F4J at 0.568 Å mean and 0.989 Å max, β2AR 100 onto 3P0G at 0.829
and 1.289, and P4 still reads +0.305, −0.368, +0.593, −0.934. Those are the
same numbers to three decimals, which is the check that mattered before
touching shared code with four rows in flight.

The thorough fix, intersecting the cores, is still not available and the reason
is worth keeping: these are residue numbers, and ABL1's references are split
between the ABL1a and ABL1b numbering nineteen residues apart, so the
intersection across its nineteen references is **three residues**. It would
have to be taken over KLIFS positions and resolved back per structure. The
position-level intersection is 38 of 42 for CDK2 and would be the better core;
it is a change to how cores are built rather than to how one is chosen, and it
is not needed to get twenty-eight references out of twenty-nine.

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


## 4v. The DSSP ribbon was missing the half that makes it a comparison

Section 6.2 asks for a "DSSP ribbon: residues x time heatmap with motif bands
marked, plus active/inactive reference consensus rows above it". The heatmap
was there. The reference consensus was not, which leaves a chart that says what
the trajectory's secondary structure is and nothing about whether that is what
an active receptor looks like. The comparison is the point.

The data had been computed all along. `report.sse_profiles` holds, per state,
per motif, per position within that motif, a count of DSSP codes across the
whole reference ensemble: `{"H": 70, "T": 2, "E": 1}` and so on. Two steps were
missing. Pick the modal code, and resolve the motif position to a residue
number through `motif_residues`, which is the same index space. Both are cheap,
and neither had been done.

It is drawn as two columns to the left of t = 0 rather than as a separate
chart. A separate chart would mean the browser keeping two residue axes in
step, and they would come apart the moment either scrolled. Three things had to
be got right for the columns to be readable rather than merely present:

- **Width.** A heatmap takes each cell's width from the gaps between
  neighbouring x values, so at one frame step apart the consensus columns came
  out as 0.1 ns slivers against an axis tens of nanoseconds long. The width is
  now a fraction of the run, which makes it constant in pixels whether the run
  is 10 ns or 50, and the tick labels always have room. A fixed multiple of the
  frame step does not: it shrinks as the run grows.
- **Separation.** An empty column sits between the consensus and the first
  frame, with a rule drawn through it, so that two kinds of column that mean
  different things do not touch. Without it the consensus reads as the opening
  frames of the run.
- **Honest gaps.** Only the motif residues have a consensus, because the
  ensembles are built per motif. The other 220 residues of the construct are
  null and draw as blank. Filling them with anything would be inventing a
  reference the build does not have.

Each cell reports how many of the references agree, because a residue where 70
of 72 agree and one where 38 do are not the same claim. The test asserts that
a 50/50 split still reports, with `agreement` of 0.5, rather than being
discarded as indecisive: a reference set that disagrees with itself about a
residue is a result about that residue.


## 4w. The state badge strip, and a shape nobody could see

Section 6.2's explorer list ends with "State badge strip: one badge per motif
for the current frame". There was no strip. The motif table gives every motif's
label at the end of the run, which is a different thing: the strip is supposed
to answer "what is this frame" while the scrubber moves.

It is built now, and it is rendered server-side at frame 0 rather than drawn by
JavaScript on arrival, so it reads correctly with no JavaScript and does not
pop in. That needed `start_label` in the summary, which was not there: the
summary carried `start_s` and `final_label` but no label for the first frame.
Rendering the final label instead would have put the strip in disagreement with
the viewer next to it from the moment the page loaded, and in a way a reader
has no way to notice.

The slider indexes the viewer's strided frames and the labels are per analysis
frame, so one has to be mapped onto the other. That mapping already existed in
the state map's click handler, in the other direction.

### The shape was computed, written down, and never drawn

`web.py` had `STATE_SHAPES = {"active": "circle", "inactive": "square",
"intermediate": "diamond", "conflict": "triangle"}`, passed into every run page
render, and used once to write `data-shape="diamond"` onto a span. Nothing read
it. There was no `.md-state` rule anywhere in `brand.css`, which is the only
stylesheet: the state labels were unstyled text and the attribute was decoration
on decoration.

Section 6.1 asks that "every state is also distinguishable by shape or label,
never by colour alone". That was satisfied, by the label, which is why nothing
ever failed and why this survived a page QC pass and an axe-core run. The
requirement was being met by the half of it that happened to work.

The marker is now drawn, and from the state class rather than from the
attribute. The shape follows from the state, so a separate map is a second
place for one fact to live, which is the thing this build keeps removing: the
map and the attribute are gone, and a test asserts they stay gone. The badges
now carry a shape, a colour and the word, which is more than the spec asks for
and is what makes a row of five of them readable at a glance while the slider
is moving.


## 4x. A shared colour palette that resolved to nothing

The DSSP ribbon's last clause is "with motif bands marked". There were none:
the endpoint sent a `motifs` object, the chart fetched it, and nothing drew it,
so the ribbon's y axis ran from residue 23 to 345 with no way to tell which
rows are TM6.

Drawing the bands needed colours, which is how this turned into something
larger. Section 7's `motifs` table has a `token` column, `store.MOTIF_TOKENS`
has been filling it since the schema was written, and the comment above it
explains the point:

> Section 6.1's motif colour tokens. Shared by every view, so they live with
> the data rather than in a stylesheet: the explorer reads a motif's token from
> the same row it reads the motif's residues from, and a motif cannot then be
> one colour on the timeline and another on the viewer.

None of that was true. `--motif-1` to `--motif-6` were stored on every motif
row and **were not defined anywhere in `brand.css`**, which is the only
stylesheet this app has. Nothing read them either. A `var(--motif-1)` would
have resolved to nothing, and `getComputedStyle` returns an empty string for a
property that does not exist, which is the quiet kind of wrong.

**And the map did not match the motifs.** I first wrote that it had been
written against a stale planned list. That is wrong, and checking section 6.1
rather than assuming is what caught it: the spec fixes the tokens in a table,
and `MOTIF_TOKENS` was that table, exactly. What diverged is the other side.
The spec assigns colours to the R-spine, the ionic lock, the PIF connector and
the toggle switch, and this build measures all four but as pocket states,
micro-switch metrics and corroborating rules rather than as motifs with
reference ensembles, so none of them has a motif row or anything to colour.
Meanwhile `c.l`, `TM5` and `TM7` are scored as motifs and the table does not
name them:

| family | measured, with no token | token, but not a motif |
|---|---|---|
| kinase | `c.l` | `r_spine` |
| receptor | `TM5`, `TM7` | `ionic_lock`, `PIF`, `toggle` |

Three of the eleven motifs this build scores had no colour, and four tokens
were spent on things that are measured but are not motifs: `r_spine` as a
pocket state, the ionic lock as section 1.2's TM3-TM6 distance, the PIF
connector with its own number in section 4e at 3.29 Å.

**So this is a deviation from an explicit table, and it is kept minimal.**
Every assignment section 6.1 makes that applies to a motif this build renders
is honoured, and only the tokens left free by a motif the build does not render
are reused. For the kinase that is tokens 1 to 5 unchanged, with `c.l` taking 6
from the R-spine. For the receptor it is 1, 3 and 6 unchanged, with `TM5`
taking 2 from the ionic lock and `TM7` taking 4 from the PIF connector, and 5
left unused because nothing replaces the toggle switch.

My first version moved helix 8 from token 6 to token 5 as well, for no reason
beyond tidiness, which is a change to the spec's table that buys nothing. It is
back on 6. Reading the clause rather than working from the shape of the problem
is what caught it, and it is a fair warning about the rest of this section: the
case for changing something is weaker than it looks once the document turns out
to have been specific.

All three are fixed together, because they are one thing. The variables are
defined, in Okabe-Ito, which is the usual colour-blind-safe set and is what
section 6.1 asks for. The map is now exactly the motifs that become rows, with
a test that asserts it both ways: no measured motif without a colour, and no
colour for something that is not measured. And both charts read the token
rather than taking what the plotting library offers next, so TM6 is the same
blue on the ribbon's band and on the state timeline, which is the thing the
comment promised and nothing delivered.

The bands themselves are laid out in lanes, because the ranges overlap: NPxxY
is residues 322 to 326 and sits inside TM7's 304 to 328, so a band drawn over
another band hides it. Each motif takes the first lane that is free at its
start residue. A motif shorter than twelve residues gets its name written
beside the bar rather than turned on its side inside it, because "NPxxY" is
longer than a five-residue bar is tall.

One thing worth stating plainly: Plotly's autorange looks at traces and not at
shapes, so the first version drew all five bands outside the axis and clipped
every one of them. The range is set explicitly now. A chart that silently
declines to draw what it was given is the same failure mode as the rest of this
section.

`tools/make_figures.py` had the same problem in a fourth place: it coloured
each motif by its index in a list of its own, so a motif's colour in a FINDINGS
figure had nothing to do with its colour on the page. It reads the tokens now,
resolving them against `brand.css` itself rather than against a copy of the
values, so the stylesheet is the one definition and a figure cannot drift from
the site it illustrates.

### And the same table names four more tokens

Section 6.1's table has a last row: `--state-active`, `--state-inactive`,
`--state-intermediate` and `--state-conflict`, "state badges, map colours".
Those were not defined either, and the first version of the state markers in
section 4w put the four hex values straight into the CSS rules. That is the
same fault as the motif tokens, written one paragraph later in the same
stylesheet, while I was in the middle of writing up the motif tokens as a
fault. They are defined now and the markers reference them, and
`tools/make_figures.py` reads them for the reference state map, which is the
other place the spec names them for.

The state map on the page had the same gap from the other direction: it set
point shapes, correctly, and no colours at all, so Plotly chose from its own
cycle. The figure version had just been changed to read the tokens, so the two
maps would have shown the same references in different colours. Both read the
tokens now.

Reading the clause is the lesson of this whole section. Five separate colour
decisions here were made by guessing at what the build needed, and every one of
them would have been made correctly by looking at a table that had been sitting
in the spec the whole time.

---

## 4y. Section 2 names five files and none of them existed

Section 2 is the hardware section and it is specific. `tools/hwprobe.py` writes
`config/hardware.toml`; section 2.2 derives `config/tuning.toml` from it;
measured speeds go in `config/benchmarks.toml`; and section 2.3's analysis pool
polls `data/manifests/md.jsonl`. There was no `config/` directory at all, no
`hwprobe.py`, and `data/manifests/` held the two reference manifests and
nothing else.

**The settings were right anyway**, which is why this survived so long.
`tools/run_md.py` works the platform out at run time rather than reading it
from a file, and does it better than the spec's table asks: section 2.2 expects
OpenCL at mixed precision, Apple's OpenCL is single precision only, and the
runner tries each candidate against the real system and drops to what the
device supports rather than failing. One GPU worker, analysis on the CPU
alongside, 4 fs with HMR. All of it correct, none of it written down.

What was missing is the record. Section 2.2 ends "Every compute budget in
Section 5.6 is derived from these measurements, never assumed", and a budget
cannot be derived from a file that does not exist or checked against one. The
probe now writes all three:

| file | what it holds |
|---|---|
| `config/hardware.toml` | Apple M2 Ultra, 16 performance and 8 efficiency cores, 128 GB, 1,369 GB free, a 60-core GPU, and the three OpenMM platforms this machine has with the speed each claims |
| `config/tuning.toml` | `md_platform` OpenCL, `md_precision` single, `cpu_workers` 14, `io_workers` 8, `concurrent_md` 1, `memory_headroom_gb` 24 |
| `config/benchmarks.toml` | measured ns/day per row, with the finished ones marked |

Three things worth saying about how it was built rather than what it writes.

The platform is **read, not re-benchmarked**. Section 2.2 asks for a 20 ps
benchmark of each platform, and `run_md.py` already tries each candidate
against the real system and records what it kept, so the probe reads
`run.json`. Running the benchmark again would mean competing with the queue for
the one GPU this machine has, to re-derive a number the queue has already
measured on the actual system rather than a 20 ps stand-in.

`md_precision` records `single`, not the `mixed` the table asks for. That is
the D010 deviation, and writing `mixed` into a file that claims to describe
what ran would be worse than the deviation itself.

And a key that cannot be measured is **omitted**, not defaulted. A machine with
no `hw.perflevel0.physicalcpu` gets no `cpu_workers`, rather than a 2 that
reads like a measurement. That is section 0's rule and it matters most in a
file whose only purpose is to be trusted later by someone who was not here.

The test for this caught me writing the same `relative_to` crash that section
4o records from `tools/make_figures.py`, where `--out` outside the project
kills the log line. The same bug, written again three hours later, in the same
shape of function. It now has the same helper, with a comment saying so.

### The fifth file, and why it starts with two reconstructed lines

`data/manifests/md.jsonl` is the other one. Section 0 asks for a manifest per
stage and section 2.3 names this one specifically. The queue now appends a line
when a row starts, finishes or fails.

Append-only, deliberately. It is a history and not a view of what is true now:
a row that started and then failed is two lines, and collapsing them to one
would throw away the thing most worth knowing. `tools/ingest.py` reads the run
directories for current state; this sits beside it as the record of what
happened. The write is wrapped so that a queue which cannot write its own notes
still runs the simulations, because losing the record is bad and losing a day
of GPU time to an unwritable directory would be worse.

The manifest starts with two lines marked `reconstructed`, taken from
`production.log` and a finished `run.json`. The queue process running this
campaign loaded its code before this existed, and Python does not reload a
running module, so the rows it is working through now will never append. The
honest options were to leave the file empty until the next campaign or to write
down what can be recovered. It writes it down and flags it, because a line
worked out afterwards is not the same kind of fact as one written as it
happened, and a log that blurs the two is worse than one that starts late.

Restarting the queue to pick up the change was the other option and it was not
taken: `b2ar_apo` is thirteen nanoseconds into a fifty nanosecond row, it would
resume from its checkpoint, and interrupting a good run to improve its
bookkeeping is the wrong trade.

---

## 4z. Two dependencies added on an estimate nobody measured

Section 3's package table is primary and fallback in two columns, and D004
recorded which of the spec's packages were added, which were left out, and why.
It is a careful entry. Three of the five it added are not used by anything.

Found by grepping for every name in the table rather than reading D004 and
believing it, which is the only reason it came up at all.

**`mdanalysis`**, added for "membrane analysis and helix axes". The membrane
analysis is `flexappeal/refstates/membrane.py`, written against mdtraj. The
helix axes are `bundle_tilt` in the same file, which fits each transmembrane
helix with SVD and averages the seven, and that is section 3's own named
fallback for HELANAL: "hand-rolled helix axis fit (Kahn-style local axis)". It
was written in section 4n to fix the tilt measurement, without my noticing the
spec had named it as the fallback for a library already in the environment. The
choice is fine. Not recording it was not.

**`pyarrow` and `duckdb`**, added for "the per-residue DSSP store, which is too
wide for SQLite and too large for JSON". No such store exists. The per-residue
DSSP went into the analysis JSON, which is what `web.py` serves the ribbon
from, and measured on `b2ar_apo` at 135 frames it is **46 KB** inside a 575 KB
payload.

That is the part worth recording. Two dependencies were added to a lock that
also serves the deployed application, on a size estimate that was out by orders
of magnitude, and the estimate was never checked against a file. It would have
taken one `len(json.dumps(...))`. The same build measured the analysis payload,
the endpoint timings, the web bundle and the per-chunk bilayer cost, carefully,
and did not measure the one number that justified a dependency.

They are left installed. Removing three packages from a working lock to save
nothing is the trade D004 itself declined when it refused BoltzMaker's pins.
What is corrected is the reason recorded beside them.

---

## 4aa. A clean page called broken by a screenshot

Three of the sixteen page tests failed on a run while the GPU was carrying a
50 ns membrane row, and all sixteen passed in isolation a minute later. The
failure was not what the test is for. Its assertion is that the page loads
with no console errors, and that passed; what failed was capturing the
screenshot afterwards, with the browser protocol's "Unable to capture
screenshot" under load.

The test already knew about this shape of problem. It falls back from a
full-page capture to a viewport one, with a comment saying why:

> A viewport shot is still the artefact section 10 asks for, and failing the
> console-error assertion over a screenshot would be the tail wagging the dog.

The fallback had no fallback. When the viewport capture failed too, the
exception propagated and a clean page was reported as broken, which is the
exact thing that comment argued against, one level further down.

The capture is best effort now: both attempts are tried, a failure prints what
the browser said, and the test passes or fails on the page rather than on the
artefact. A missing screenshot is worth saying out loud and is not evidence of
anything being wrong with the page.

Worth recording because of what it nearly cost. I had committed on the previous
run and read "3 failed" afterwards, and the first plausible reading of three
page tests failing is that something I had just changed broke three pages. It
had not; the machine was busy. A test that fails for a reason unrelated to its
assertion spends someone's attention every time it does, and the cheapest
moment to fix that is the first time it happens.

## 4ab. Section 12's two files, one written and one not

Section 12 gives the contents of `.claude/settings.json` and a project
`CLAUDE.md`. Neither existed. One is written now and one is not, and D024 has
the reasoning: `CLAUDE.md` is standing instructions and is written, with the
two lines the user overrode during the build changed and marked;
`.claude/settings.json` is a permissions policy whose deny list includes
`Bash(git push:*)`, which the user has explicitly overridden, so writing it
either breaks the authorised workflow or means choosing what a future session
may do on this machine without being asked.

## 4ac. The fifth test in section 10's list

Section 10 names five pytest items. Four existed: the numbering maps against
BoltzMaker's fixtures, the helix rotated by a known angle, the salt bridge at
3.5 angstrom, the state-score symmetry, and the resume. The fifth, "membrane
builder smoke test on a 5 ps run", did not.

Section 5.3's Tier 1 smoke covered the same ground two hundred times over, at
1 ns on the same system, and `run.json` has the result. What was missing is the
version that lives in the suite and runs again whenever anyone changes the
builder.

It runs against the system the membrane builder actually produced rather than
a fixture, because the builder's output is the thing being tested. That needs
the GPU, which has been carrying a 50 ns row for most of this build, so it
skips with a stated reason while a run is in flight rather than competing for
the device or spending ten minutes on 187,848 atoms on the CPU.

**A test that has only ever skipped is not a test**, which is most of what this
section has been about, so the body was proved separately: the same code on the
CPU platform with ten steps instead of 2,500, to show the XML deserialises, the
topology matches the positions, a context builds and stepping does not raise.
It loads 187,848 particles, which is `build.json`'s figure for that system to
the atom, and comes out at −2,353,259 kJ/mol with every coordinate finite. All
three of the test's assertions hold on that run.

## 4ad. Phase 0's acceptance criteria, checked at last

Section 8's phase 0 has two testable criteria and neither had been run as
written.

The first is literal: `pixi run python -c "import openmm, mdtraj, MDAnalysis,
gemmi"`. It passes, which is slightly funny given section 4z: MDAnalysis is in
that import line and nothing else in the build uses it.

The second is "BoltzMaker's own `sse_comparison` tests pass from the vendored
copy". **They do, 50 of 50.** I had never run them. What I had run was the
whole vendored suite, which reports 17 failures because parts of it shell out
to `vendor/boltzmaker/.venv/bin/python3`, an interpreter that does not exist
inside this pixi environment. `pytest.ini` excludes the vendored tree for
exactly that reason, and the comment in it, as written, reads as though
BoltzMaker's tests are unusable here. The ones this build depends on are not:
they pass completely, and the file now says so with the command to repeat it.

The distinction matters because `sse_comparison` is the thing section 4.1 says
to reuse and never reimplement. "Its tests fail in our environment" and "its
tests pass, and a different part of its suite needs an interpreter we do not
have" are different claims about a dependency the whole build rests on, and
only the second one is true.

## 4ae. The resumed-run fault, demonstrated instead of asserted

Section 4i records a resumed run losing its barostat and nothing saying so. The
fix was an ordering one: the force has to be on the System before the Context
exists, and a resume skips the branch that adds it. The test for it reads the
source and says so in its own docstring, because the real path needs a built
membrane system and a GPU.

What that test could not show is **why** getting the order wrong is invisible,
which is the entire reason the fault survived long enough to produce four and a
half nanoseconds at constant volume. `abl1_apo` resumes from its checkpoint
tonight, unattended, so "the source is in the right order" felt like thin cover
for a path that has never been executed since it was changed.

It turns out eight particles are enough. A checkpoint written from a system
with a barostat loads **without complaint** into a system that has none, and
the context then reports the barostat's own parameters:

    destination forces:      ['NonbondedForce']
    destination parameters:  {'MonteCarloPressure': 1.0,
                              'MonteCarloTemperature': 300.0}

So a resumed run prints the pressure it is supposedly holding while having
nothing to hold it. That is now a test on the Reference platform that runs in
milliseconds, alongside the same sequence done correctly, which keeps the force
and the parameter together.

The source-order test stays. The two are different claims: one is that
`run_md.py` is written correctly today, the other is that OpenMM really does
restore a parameter whose force is gone. The second has an assertion attached
saying that if OpenMM ever stops doing it, the fault becomes self-announcing
and the test should be revisited.

---

## 4af. Rehearsing the pass that nobody will be awake for

Every piece of the streaming loop's final pass had been checked on its own.
The sequence had not, and it is the part that runs unattended: when the queue
writes `queue done`, the loop breaks and runs the bilayer analysis, the
hardware probe, the figures, the verification and the controls, back to back,
with nobody watching. So it was run, with the loop stopped so nothing competed.

It completes, exit 0, and two things happened for the first time.

**P6 passes on the live database.** It has read `pending: no membrane chunks
ingested yet` for the whole build, because the streaming passes use
`--skip-membrane` and only the final pass does the bilayer. Now:

    P6: every membrane chunk inside section 5.8's windows
        b2ar_apo     2 chunks, APL 67.7 A^2 (67.3 to 68.0), 38.6 A  PASS

That is section 9's sixth control, green, on real production data rather than
on the rehearsal database of section 4n. Five of six are now decided and only
P5 is waiting, on a row that has not started.

**And the bilayer figure drew for the first time**, which is worth more than
the control line because it shows the direction. Area per lipid reads 68.0 then
67.2 A^2 across the two chunks, with the receptor tilted 14.9 degrees. Section
5.2 measured 68.4 on the Tier 1 run and said it was "still falling from where
the builder left it, which is what an equilibrating bilayer should do". Ten
nanoseconds later it is 67.2, inside the expected window rather than 0.4 above
it. The prediction was right and now has three points behind it.

### The verification step was announcing a failure on every success

The rehearsal found one fault, in the step added an hour earlier:

    [04:04:08] every stored number agrees with the analysis file it came from
    verification reported a disagreement, see above

Both lines, on the same run. The block was

```bash
pixi run python ... report --verify \
  | grep -E "agrees|disagree|..." \
  || echo "verification reported a disagreement, see above"
```

and `||` on a pipeline takes the **pipeline's** status, which is grep's. `log`
writes to stderr, so the "agrees" line never reached the pipe at all; grep
matched nothing in the report table, exited 1, and the alarm fired. Every time.

The same construct fails the other way round too, which is the worse half. On a
genuine disagreement `report --verify` exits 1, but grep would match the word
"disagree" and exit 0, so the `||` would not fire and the failure would print
nothing. An alarm that cries wolf on every success and stays silent on the one
failure it exists for.

It now tests the command's own status, with both streams captured to a file,
and both branches were run rather than reasoned about: against the real
database it prints the agreement line and nothing else, and with one score
moved by hand in an analysis file it prints

    VERIFICATION FAILED; the first lines of data/md/verify.log:
        abl1_apo g.l.start_s: database -0.0204, analysis file 0.98765

naming the run, the motif, the field and both values. The analysis file was put
back afterwards.

That fault would have survived to the morning. It is only visible by running
the sequence, because every piece of it passes on its own.

### A checkpoint is bound to the platform that wrote it

Found while trying to answer a narrower question. `abl1_apo` resumes the moment
`b2ar_apo` finishes, and its checkpoint was written three minutes before its
trajectory stopped, so frames exist past it and the resume has to drop them.
Rather than infer how many from file timestamps, the checkpoint can be read:
deserialise `system.xml`, build a context, load it, ask for the time.

Not on the CPU, it turns out:

    OpenMMException: loadCheckpoint: Checkpoint was created with a different
    Platform: OpenCL

Which is a better finding than the number would have been. `run_md.platform()`
tries CUDA, then OpenCL at mixed precision, then OpenCL at single, then CPU,
and that fallback is right for starting a row. For resuming one it is a trap: a
GPU that has gone away between a row starting and its resume means the runner
picks CPU, and `loadCheckpoint` then throws. Unguarded that is an OpenMM
traceback at the top of a queue with 24 hours of work behind it, the row's
nanoseconds still sitting on disk, and nothing anywhere saying they are still
good.

It is guarded now, and deliberately not recovered from. Starting the row again
would throw away the progress in the checkpoint without being asked, and
section 0's ladder is about degrading with the reason recorded rather than
quietly redoing work. So it says which platform the checkpoint needs, says the
checkpoint itself is fine, and stops; the queue records the row as failed with
that line as the reason, which is section 4q's machinery doing its job.

The question that started this is still unanswered: reading the checkpoint
needs the OpenCL device that production is using, and a second context on the
one GPU is exactly what section 2.2's `concurrent_md: 1` exists to prevent. The
truncation it would have confirmed is covered by its own test, and the resume
will log the figure when it runs.

---

## 4ag. The resume ran unattended, because I killed the queue myself

Section 4af left one question open: where exactly would the `abl1_apo` resume
restart, given a checkpoint written three minutes before its trajectory
stopped. Reading the checkpoint needed the GPU, and the GPU was busy, so the
answer waited.

It stopped waiting at 04:19 on 11 October, when I stopped a background shell I
had read as a stale poller. It was not stale. Its command ended in `exec
caffeinate -dimsu pixi run python tools/queue_md.py`, so the shell I stopped
*was* the production queue, and `b2ar_apo` went from 17.2 ns of 50 to no
process at all. Nothing in the shell's one-line summary said so; the summary
showed the heredoc that wrote `tools/after_smoke.sh`, and the `exec` was the
last of twenty lines. The lesson is cheap to state and was not free to learn:
read what a shell actually runs before stopping it, not what its first line
suggests.

What it bought was the rehearsal. The resume path for a membrane row had never
executed. At 04:20 it did, unattended, and every piece of it behaved:

    [04:20:25] barostat re-added for the resumed run
    [04:20:26] platform: OpenCL {'Precision': 'single'}
    [04:20:26] resumed from checkpoint at 17.000 ns
    [04:20:26]   dropped 2 trajectory frame(s) and 2 state row(s) written past the checkpoint
    [04:20:26] producing 33.0 ns (8,249,999 steps)

So the open question has its answer, on a different row than expected and with
a real number rather than an inference: **17.000 ns, with two frames and two
state rows dropped.** The barostat re-add from section 4ac ran on a resumed
membrane row for the first time. The platform guard from section 4af did not
fire, which is the right outcome: the OpenCL device that wrote the checkpoint
was still there, and the guard exists for the case where it is not. Cost of the
whole episode was 0.2 ns of discarded sampling and two minutes of idle GPU.

### A row in flight runs the code it started with

The restart turned up something that would have mattered more if it had gone
the other way. `data/md/b2ar_apo/` had no `checkpoint_phase.json`, although
`save_checkpoint` writes one on every save. The reason is that the row started
at 16:07 and the marker was added to `tools/run_md.py` at 21:47, in commit
f225a33. Python had already loaded the module. **A fix committed while a row is
running does not reach that row**, and the queue's own 24-hour rows are long
enough for that gap to be ordinary rather than exotic.

Here it was harmless because `checkpoint_phase()` treats an absent marker as
production, which is what its docstring says it is for. Had the default been
`equil`, the restart would have thrown away 17 ns and repeated a 10 ns
equilibration, and the log would have looked entirely normal while doing it.

That left the clock question. A resume that mistook production time for total
time would have read 27.2 ns rather than 17.2 and stopped the row 10 ns short,
and the code that resets the clock was last touched in the same post-21:47
commit, so the running row could not be assumed to have it. The row's own
`state_data.csv` settled it without the GPU: step 4,300,000 at 4 fs is
17,200 ps exactly, so the clock had been reset at production start. Then the
resume read 17.000 ns and agreed.

### A row certified itself on the number it was aiming at

Checking the clock question raised a worse one. `run.json` recorded
`production_ns` as **the target the row was asked for**, and `queue_md.complete`
then checked that number against the same target:

    return float(done.get("production_ns", 0.0)) >= target_ns

Which cannot fail. A row that stopped at 40 ns of 50 still wrote 50 into
`production_ns`, and `complete` certified it, and `tools/ingest.py` ingested it,
and the explorer showed a 50 ns row. Run against the short case the old
expression returns `True`; the new one returns `False` and the queue resumes the
row to its target on the next pass.

So `run.json` now carries `produced_ns` as well, read off the context at the end
rather than assumed:

    produced_ns = simulation.context.getState().getTime().value_in_unit(
        unit.nanoseconds)

`complete` prefers it, falling back to the target for rows written before the
field existed, and a row more than 50 ps short records
`status = "short:40.00_of_50ns"` instead of `complete`. The tolerance covers
integer step counts rounding sub-picosecond and nothing larger. It is a literal
in both files, because the queue should not build an OpenMM system to read one
float, and a test parses both and asserts they agree.

The same exit had a second fault. The branch for a row already at its target

    log(f"already at {done_ns:.1f} ns of {target_ns}; nothing to do")
    return 0

announced itself and wrote nothing. A row resumed to exactly its target would
have had no `run.json` at all, and `complete` and `ingest.py` both read an
absent `run.json` as a row still in flight, so a finished row would have shown
as running for good. It writes one now, through the same writer as the normal
exit, omitting `ns_per_day` and `wall_hours` rather than zeroing them because a
row that did no stepping has no speed to report, which is section 0.3's rule.

Section 4ag's own lesson then applied to the fix within the hour. `b2ar_apo` is
in flight on the `run_md.py` that existed at 04:20, so it will write the old
`run.json` with no `produced_ns`, which is exactly the case the fallback was
written for. The live queue process is likewise the 04:20 `queue_md.py` and
carries the old `complete`; that one is harmless, because `complete` only
decides what to skip at the head of a row and each of the four runs once per
invocation. The next invocation gets both fixes.

## 4ah. A pre-flight check, because the questions kept being asked by hand

Deciding whether the restart in section 4ag was safe meant asking the same
questions of a row on disk, one `sed` at a time: is there a built system or will
it be built, is there a checkpoint, was it written during equilibration or
production, where has the production clock got to, how many frames ran on past
the checkpoint, which barostat does the resume add, and which platform does the
checkpoint need. Every one of those decides whether a row resumes or quietly
starts again, and the queue had no way to state any of them.

`tools/queue_md.py --preflight` states all of them and exits. It reads only: it
builds nothing, steps nothing and opens no OpenMM context, so it is safe beside
a production row, which is the only time anyone wants it. A row that cannot run
at all gives a non-zero exit, so it can gate a queue rather than merely inform
one.

On the four rows as they stand:

    b2ar_apo:       production checkpoint, at or below 17.4 ns of 50
                    174 frames, last written 0 min after the checkpoint
                    MonteCarloMembraneBarostat re-added
    abl1_apo:       production checkpoint, at or below 4.3 ns of 50
                    43 frames, last written 3 min after the checkpoint;
                    those frames are dropped on resume
                    MonteCarloBarostat re-added
    b2ar_carazolol: no checkpoint; starts from minimisation
    abl1_imatinib:  no checkpoint; starts from minimisation

Two things came out of running it that were worth knowing before the rows do.

**The two unstarted rows are already built.** Both carry `system.xml` and
`solvated.pdb`, and their ligands are in them: 45 CAU atom lines for carazolol
and 69 STI for imatinib. The `ligand` field's own comment in `queue_md.py`
records that `abl1_imatinib` would once have failed at build time after the
three rows ahead of it had spent their GPU hours. It will not: there is nothing
left to build, and the ligand that amber14 has no template for is already
parameterised and serialised.

**The frame-drop figure had to stop being arithmetic.** The first version
derived the resume point from the current checkpoint interval, and reported
abl1_apo at 4.2 ns when its checkpoint is at 4.0: that row was written before
the checkpoint and frame intervals were matched, so its checkpoint sits on the
old 1 ns boundary. Arithmetic from the current interval is wrong on exactly the
rows that still need resuming. What can actually be measured is how long the
trajectory went on being written after the checkpoint was, and that is what it
reports now, with the exact count left to the resume, which logs it. The
difference between the two versions is the difference between a check and a
restatement of an assumption.

## 4ai. The streaming loop would have stopped at the first row boundary

Section 4af's rehearsal exercised the final pass. It did not exercise the
condition that decides when the final pass runs, and that condition had a race
in it that the restart in section 4ag came within a minute of triggering.

`producing()` looked once:

    pgrep -f "tools/run_md.py" 2>/dev/null | grep -qvE "^(${mine})$"

The queue runs each row as a fresh subprocess. Between one row's runner exiting
and the next being spawned there is a moment with no runner alive and nothing
yet written to say the queue is done, and a single look inside that window
returns false. The loop then leaves, runs the final pass, and exits.

The consequence is the worst shape a fault can have here. Streaming would stop
silently after the first row. The final pass, which computes the bilayer QC and
decides G6, would run against a third of a trajectory and write its numbers as
though they were final. Nothing would re-run it, and nothing in the output would
say which row it had stopped at. The window is a second or two against a
half-hour interval, so it would have happened rarely, which is worse than
happening reliably: it would have happened on one row and not the others.

The restart at 04:19 opened the same window for about a minute. The loop
survived it only because it happened to be asleep.

Two changes. The check now names the queue as well as the runner, because the
queue is the process whose lifetime actually matches production and it lives
straight across row boundaries. And an absence has to persist before it counts:
twelve looks, five seconds apart, so a row boundary or a restart by hand passes
and a real stop costs one minute of delay. The PID exclusion still does the work
of not matching ourselves, which matters more now that the pattern names two
scripts and both appear in this file's own text.

Proved rather than asserted, which needed the script to be testable. It now
returns early when sourced, the way the Python tools guard `__main__`, so a test
can source it and call `producing` against a stubbed `pgrep` on PATH. The grace
is two environment variables so the test for a real stop does not take a minute
to decide nothing is there. Run against the gap case with the grace set to one
look, which is the old behaviour, `producing` returns false; with two, true.

## 4aj. Six classes the pages used and the stylesheet never defined

Giving the new `short:` status a colour meant finding the class the status
already used, `md-tag`. It is used in two templates. It is defined in none.

`brand.css` is the only stylesheet these pages load, and not one of the
reference-state templates carries a `<style>` block, so a class that is not in
that file is simply unstyled markup. Scanning the templates for `md-*` tokens
and differencing against the stylesheet found five more:

    md-stats      index cards and every run page: the figures block
    md-scrub      the trajectory scrubber under the viewer
    md-readmore   the "how to read this" disclosures, eight of them
    md-sub        a card's one-line subtitle
    md-viewer     the Mol* container

`md-stats` is the one that mattered. It is the `<dl>` of `<div><dt>/<dd>` pairs
carrying start state, length, speed and status at the top of every card and
every run page, and undefined it rendered as a browser-default description
list: each label on its own line with its value indented under it, where the
design is a row of labelled figures. `md-viewer` escaped visibly only because
the element carries `min-height:460px` as an inline style.

Nothing was ever going to catch this. The pages return 200, the console is
clean, the page tests pass, and the markup is valid. The only way to see it is
to compare the page against the design it was meant to have, which is to say by
looking, and the thing about a page that renders is that it does not ask to be
looked at.

So the difference is now a test. It reads every template, collects the `md-*`
tokens out of `class` attributes and asserts each one is a selector in
`brand.css`. Five exceptions are named and explained rather than skipped:
`md-reanalyse` is a JavaScript hook that `analysis.js` selects on, which is a
real use of a class with no style attached, and four others belong to pages
outside this work and are recorded instead of fixed in passing. A second test
asserts the allow-list has no stale entries, so an exception cannot outlive the
markup that earned it and quietly cover a later class of the same name.

The six are defined now, in the stylesheet's own idiom and from its tokens.
`md-tag-bad` is a status that is a caveat rather than a stage: `running` is a
stage and stays quiet, while `failed:` and `short:` say the numbers below
describe less run than the page claims, and they take the same red as the bad
verdict badge from one shared declaration rather than a second copy of the hex.

## 4ak. Tier 3, and the push that could not run

Section 11 asks for `deploy/runpod/` written and not executed, and section 5.3
adds what it should contain: job script, pixi lock, resume logic. Auditing the
deliverable list against the tree, I read a truncated `ls` and concluded the
directory was empty. It was not: `job.sh` and `sync.sh` were written on
10 October and are good. The lock requirement was already met too, since
`pixi.toml` carries `osx-arm64` and `linux-64` and the lock's linux build of
OpenMM 8.5.2 comes with `cuda-version-13.3` and `libcufft` behind it.

What was missing was the README and the environment example, and reading the
two scripts properly to write them turned up a fault worth the detour.

### The push was unreachable on the path that needed it

`job.sh` ran the queue, captured the status, and pushed the result to durable
storage:

    pixi run python tools/queue_md.py --only "$ROW" --production-ns "$NS"
    status=$?

The file opens `set -euo pipefail`. Under `set -e` a failing command exits the
script there and then, so `status=$?` was unreachable and the push below it
never ran. The push is how a pod's work survives the pod. So the one path where
it mattered, a row that crashed or a pod reclaimed mid-run, was the one path
where it did not happen, and what was lost was the checkpoint the row had just
spent hours earning.

Making the call part of a condition is what stops `set -e` firing:

    status=0
    pixi run python tools/queue_md.py ... --require-platform CUDA || status=$?

and the push is now outside any success branch, deliberately, because a failed
row's checkpoint is the most valuable thing on the disk. The test runs both
shapes rather than asserting about the line, since the fault is in how bash
behaves rather than in how the line reads.

### The fallback ladder is wrong on a rented GPU

`run_md.platform()` tries CUDA, then OpenCL mixed, then OpenCL single, then
CPU. That ladder is right on a machine you own: section 2.2 measured the last
step at 29.5 ns/day against 206, and a slow run costs time you already have.

On a pod billed by the hour it is a trap. A CUDA device that fails to
initialise drops a 1 µs row to a fiftieth of the speed and the pod goes on
billing, with one line in a log nobody is reading at the time to say so. So
`platform()` takes `require`, which turns the ladder off and refuses rather
than descending it, `run_md.py` and `queue_md.py` both take
`--require-platform`, and `job.sh` passes `CUDA`. Asking for CUDA on this Mac,
where it is not registered, exercises the refusal for real rather than through
a mock.

### What Tier 3 cannot do

A checkpoint cannot cross the platform that wrote it, which section 4af
established and which has a consequence here that is easy to miss: a
part-finished Tier 2 row on the Mac cannot be carried up to a pod and
continued. Its checkpoint is OpenCL and the pod is CUDA. Tier 3 rows start from
their own equilibration. Between pods it is not a problem, both being CUDA, and
that is exactly what `sync.sh` relies on.

There is no launcher in the directory and that is deliberate. Creating the pod
is a console action behind G1 and G2, and a script that creates billable
infrastructure is the wrong thing to leave in a public repository beside one
that says it has never been run.

## 4al. The one failure the round-trip verifier cannot see

Section 10's QC list reads as six items. Five have tests: the numbering maps,
the metric functions on synthetic geometries, the state-score symmetry on a
two-reference toy, the membrane builder on a 5 ps run, and the resume. The
sixth did not:

> any metric with all-NaN for a motif is a test failure

Nothing checked it. Worse, the check that looks as though it would is exactly
the one that cannot. `verify_summary` compares every stored number against the
analysis file it came from, and an all-NaN metric passes that comparison
perfectly, because the file says NaN too. Both sides agree, and they agree
about nothing.

What makes it worth a gate rather than a footnote is how it renders. An
all-NaN metric draws as an empty column and a flat line. On the page that
cannot be told apart from a quantity which genuinely did not move, and "did not
move" is the headline this build exists to report. A missing atom, a residue
the numbering failed to resolve, a motif whose span ran off the end of a chain:
each of them produces a result that looks like the finding.

So `verify_live_metrics` is its own check, against the database rather than a
fixture, for the same reason `verify_summary` is: four rows land unattended and
the question is whether what landed is measurable, not whether the code could
measure it. It reports three distinct failures.

**All-NaN.** Every frame of a metric with no value. Some NaN frames are not
reported, deliberately: a salt bridge whose atoms are absent from a handful of
frames has an honest reading in the rest, and complaining about that would make
this noise rather than a gate.

**A value with no score.** The distance is there and the state score derived
from it is not, so the motif has a measurement and no reading of it. On the
page the score column is empty while the distance column is full, which looks
like two bugs and is one.

**A motif that contributed nothing.** This one needed a second query. A motif
whose every frame failed to ingest leaves no rows for a `GROUP BY` to find, so
it passes a check built on groups while being as absent as one that is all NaN.
The group query cannot see what is not there, which is the recurring shape of
the faults in this section.

On the live database: 2,536 rows, no NULL in any of them, every motif metric
with a reading. The final pass echoes both verdicts now rather than only the
round trip's, because a check whose result never reaches the log is close to a
check that did not run, and that log is the only account of the pass anyone
reads.

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

**And it is not drifting.** Measured again at 14.6 ns over 147 frames, the
distance has gone nowhere: mean 13.02 Å over the whole run, 12.94 across the
first half and 13.09 across the second, 13.03 over the final five nanoseconds.
Still one frame below the threshold and it is an early one.

That is a more interesting result than a drift would have been. The cleft
closed during equilibration, when the agonist and the nanobody that 4LDE was
solved with were taken away, and then stopped. What the run shows so far is an
apo receptor sitting at a stable distance of its own, four angstrom below the
active crystals and five above the inactive ones, rather than travelling
between them. Section 7's timescale limitation still applies and this does not
contradict it: nothing here is the microsecond relaxation, and a plateau at 15
ns is not evidence that there is no slower motion. It does mean the early fall
is relaxation out of a bound pose and not the beginning of a transition, which
is what section 5.5 hedged about on the 1 ns data and can now be said with
fifteen times as much of it.

### 5.0b Per motif at 14.6 ns, with the whole-protein number hiding it

The first production observation that needs the per-motif breakdown to be
visible at all. β2AR, 147 frames, 14.6 ns of 50, so preliminary.

| motif | start | first half | second half | change |
|---|---|---|---|---|
| H8 | +0.474 | +0.297 | +0.120 | **−0.176** |
| NPxxY | +0.677 | +0.646 | +0.747 | **+0.101** |
| TM5 | −0.013 | −0.046 | −0.036 | +0.009 |
| TM6 | +0.534 | +0.476 | +0.437 | −0.040 |
| TM7 | +0.377 | +0.403 | +0.400 | −0.003 |
| **whole protein** | +0.319 | +0.278 | +0.256 | **−0.022** |

The whole-protein score has not moved. Underneath it, two motifs have moved in
opposite directions by roughly a tenth and two tenths, and the weighted mean of
a fall and a rise is neither. Section 2 argued for keeping the per-motif
breakdown beside the headline number rather than replacing it with it; this is
the first instance in real data where the headline alone would have said
nothing happened.

**What the two motifs are actually doing, which is less than it first looks.**

Split into thirds, H8 reads 0.328, 0.153, 0.143: a step down early and then a
plateau, not a continuing decline. Its scatter is 0.168 across the run, which
is almost the size of the whole change, so the honest description is a motif
that settled once and is noisy, not one that is drifting. The distances behind
it wander in both directions, 0.64 then 0.99 then 0.76 Å to the nearest active
reference. H8's helix fraction is 0.909 at the first frame, 0.909 at the last
and 0.909 on average, so whatever it is doing it is doing as an intact helix
rather than unwinding.

NPxxY is the cleaner signal of the two. Its thirds are 0.655, 0.672, 0.764 with
a scatter of 0.079, and it is driven by a real approach: 0.66, 0.66 then 0.42 Å
to the nearest active reference. The receptor's NPxxY motif is getting closer
to the active crystals, modestly, while the cleft below it stays put.

I had started writing this up as "H8 is drifting away from the active state",
which the thirds do not support and which the helix fraction argues against. It
is in here as a caution: a first-half against second-half comparison will call
a step change a drift every time, and the step was in the first third.

**So the halves and the scatter are now in the summary**, not computed by hand
at four in the morning for one row. `store.summarise` reports
`first_half_mean_s`, `second_half_mean_s`, `sd_s` and `n_scored_frames` beside
the endpoints it already had, which puts them where the page reads them too so
the two cannot derive the same number differently. `flexappeal refstates
report` prints them, and marks a change larger than the scatter:

    H8     start +0.474 -> final +0.257  halves +0.297 -> +0.120
                                         (-0.176, sd 0.168, larger than the scatter)
    TM6    start +0.534 -> final +0.528  halves +0.477 -> +0.436
                                         (-0.040, sd 0.061, within the scatter)

That comparison is against the standard deviation and is deliberately blunt. A
real test would compare the difference of two means against their standard
error, and the standard error of a mean over MD frames needs the
autocorrelation time, because frames 100 ps apart are not independent samples.
Against the standard deviation the bar is far too high for a genuine trend and
will never call noise a movement, which is the error worth avoiding in a line
somebody skims. The code says so rather than implying a test it does not run.

Three faults of my own while adding it. The first count read 294 frames for a
147-frame run, because `motif_frames` holds one row per metric and a motif has
two: the same rows-against-things error that section 4m records for the kinase
reference ensembles, made again four sections later. The first version of the
test asserted on the shared two-frame fixture, which is below the four frames
`summarise` needs before it reports halves at all, so it passed by skipping
everything it was meant to check.

**And the third one broke the page.** The halves need a midpoint, and I called
it `cut`:

```python
cut = span["hi"] - 0.1 * (span["hi"] - (span["lo"] or 0.0))   # 13.14 ns
...
for row in motifs:
    tail = connection.execute("... AND time_ns >= ?", (run_id, motif_id, cut))
    cut = len(scores) // 2                                    # 73
```

`cut` already held the final-tenth time boundary, computed once before the
loop. The first motif was summarised while it still did and came out correct.
Every motif after it asked the database for frames past 73 ns of a 14.6 ns run
and got NULL, so nine of eleven motifs silently lost their final-tenth mean,
and the run page raised a 500 on the first one it could not format.

Two things found it independently and the difference between them matters. The
page broke immediately and loudly, which is what a template with no `None`
guard buys you. The suite caught it too, in ten separate tests, which is what
says the fault was real and not a rendering accident. What neither would have
caught is the version of this bug where the overwritten value still formats:
the first motif had a number, and a glance at the page would have shown a
column that looked populated.

The regression test asserts on **every** motif rather than one, because an
error that starts at the second element is invisible to a test that checks the
first.

### What three faults in one change bought

Three mistakes in one small edit is a signal about the change, so the next
thing added was not more code but a check of a different kind: the analysis
JSON a run writes and the SQLite summary the pages read are two independent
paths to the same seven numbers, and they have to agree.

Run against the live database, all eleven motifs across both rows agree on
every one of `start_s`, `final_s`, `final_tenth_mean_s`, the two half-means,
the scatter and the frame count. Zero disagreements, so the ingest is sound end
to end and not only in the parts that have unit tests.

It is a test now, and it was checked the way section 4n's rod test was: putting
the `cut` shadowing back makes it fail, taking it out makes it pass. That
matters more than it usually would here, because this is the check that would
have caught the fault **on the first run** rather than when the page raised a
500 — and, unlike the 500, it would also catch the version of that bug where
the overwritten value still formats into something plausible.

The unit tests in this build mostly check that a function does what it says.
This one checks that two things that should say the same thing do, which is a
different and harder property, and it is the one that a silently wrong number
violates.

It runs on the real data too, not only on a fixture. `flexappeal refstates
report --verify` does the same comparison against whatever is on disk, and the
streaming loop's final pass calls it, so when the four rows land unattended the
run says whether their numbers survived the trip rather than assuming it.
Reported rather than fatal: a disagreement is worth finding in the morning, and
stopping the pass over one would lose the figures and the controls with it.

Against the live database right now it reports **every stored number agrees
with the analysis file it came from**, across eleven motifs and two rows.

There is a test that the verifier can fail, as well as that it passes. Moving
one score in the analysis file by hand makes it complain about exactly that
field, because a checker that cannot fail is worth nothing and this one runs on
four rows with nobody watching.

**And a fourth, found by the same run.** Fixing the shadowing took the suite
from ten failures to four, and the four left were the new column itself:

```jinja
{%- if row.sd_s is not none %}
```

`summarise` omits the scatter for a run too short to halve, which is section
0's rule working exactly as intended. On a dict with no such key Jinja returns
`Undefined` rather than `None`, and `Undefined is not none` is true, so the
guard let the missing key straight through to the format filter. `row.get("sd_s")`
is the version that copes with absent as well as null.

Worth recording because the two halves of that pairing are both correct on
their own. Omitting a key you cannot produce is the right rule for the data.
Treating a missing key as null is the right instinct in a template. Putting
them together needs the template to test for absence and not for null, and
`is not none` reads like it does that and does not.

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
| explorer first paint | < 2 s | **432 ms** (176 to 456, median 396, before the explorer grew) |
| frame scrub | < 100 ms | **0.1 ms** median handler cost, 0.5 ms worst |
| per-system web bundle at 50 ns | ≤ 150 MB | **3.4 MB** kinase, **3.5 MB** receptor |

These were measured once on Tier 1 and again after sections 4n to 4x added a
bilayer toggle, two reference consensus columns and five motif bands to the
ribbon, a badge strip on every slider move, and a thumbnail to each card. The
first-paint cost of all of that is 36 ms.

The bundle figure moved further, from "about 1 MB for a finished row" to 3.5
MB, and both the old number and the new one are inside the budget by a factor
of forty or more. Two things changed: the receptor gained a 64 KB lipid pair,
and the first projection was linear when the directory is not. The topology
PDB is a fixed size whatever the run length and only the coordinates and the
analysis JSON grow, so projecting the whole directory from a 4.3 ns run
overstated one system and understated the other. Section 4k has the split.

First paint is measured from the browser's own `first-contentful-paint` entry
rather than by timing the navigation call, because the target is about when a
reader sees something and not when the last asset arrives. The viewer and the
plots load after first paint by design, so a slow Mol\* does not count against
it.

The bundle is far inside its budget because the viewer coordinates are protein
and ligand only at a 1 ns stride, with the bilayer's 476 phosphorus atoms
alongside for the membrane toggle. The solvated trajectory a membrane row
produces is about a gigabyte, and almost none of it is drawn.

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

---

## 6. The positive controls

Section 9 names six. **Five are decided** and one is pending on data that is
still being generated. Everything here
is `flexappeal refstates controls` output rather than a transcription, and the
command is the thing to re-run rather than this table to trust.

| control | what it asks | result |
|---|---|---|
| **P1** | Reference structures run through the trajectory code path reproduce their own KLIFS DFG labels | **17 of 17** |
| **P2** | The same, for the αC-helix labels | **17 of 17** |
| **P3** | The β2AR cleft threshold separates the catalogue's states, leave-one-out | **100 of 100**, needs 90% |
| **P4** | Each start structure scores as the state it was built from at t = 0 | **4 of 4** |
| **P5** | An inactive-start run stays inactive over its final 20 ns | **pending** |
| **P6** | Every membrane chunk inside section 5.8's windows | **2 of 2 chunks**, APL 67.7 Å² |

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

**P6 is decided**, as of the final-pass rehearsal in section 4af. The streaming
passes run with `--skip-membrane`, so no membrane chunk reached the live
database for most of the build and P6 had nothing to read; running the pass
puts them there. On `b2ar_apo` at 16 ns it reads two chunks, 67.3 and 68.0 Å²
per lipid against an expected 60 to 68, and 38.6 Å of P-P thickness against an
expected 36 to 41. The area per lipid is still falling from where the builder
left it, 68.4 on Tier 1 and 67.2 by the second chunk, which is what section 5.2
predicted an equilibrating bilayer would do.

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
