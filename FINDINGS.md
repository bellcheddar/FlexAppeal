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
5D5A itself as the nearest inactive reference at **0.00 Å** — which is what it
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

## 5. Production runs

Pending. This section will carry the per-row speed, the energy drift, the
membrane quality checks against section 5.8's windows, and the state traces.
