# Gates waiting on Marc

Written per section 0.2's gate protocol. Everything these gates do not block is
still running: Tier 1 is finishing, Tier 2 is chained behind it, and nothing
below needs an answer before those complete.

---

## G4 Publish — waiting

**What it blocks.** Any `git push`. `main` is **74 commits ahead of
`origin/main`** and none has left this machine. That is the whole
reference-states build, from `phase 0: vendor BoltzMaker` onwards: nothing from
it has ever been pushed, because G4 has been closed since it started.

This figure is measured with `git rev-list --count origin/main..HEAD` on each
update, because it has been wrong before: it said 23 twice, written from memory
rather than from `git`, in the one file whose purpose is to carry numbers you
can act on. It was 57 earlier in this build and is 74 now.

**What is done.** Phases 0 to 4 and 6 are complete, Phase 5 is running, Phase 7
is drafted. The spec's own acceptance criteria are met for every finished phase:

| Control | Result |
|---|---|
| P1 DFG label against KLIFS | 17 / 17 (spec wants ≥ 90 %) |
| P2 αC label against KLIFS | 17 / 17 (≥ 85 %) |
| P3 TM3–TM6 state separation, leave-one-out | 100 / 100 at 11.39 Å (≥ 90 %) |
| P4 each start scores as its own state | 4 / 4 |
| P6 membrane integrity | passes on a real bilayer: APL 68.9 Å², thickness 37.6 Å |
| P5 inactive start stays inactive | the one control still genuinely waiting; needs `b2ar_carazolol` or `abl1_imatinib`, which are 3rd and 4th in the queue |
| Manifests re-derived from source | ABL1 20/20, β2AR 116/116, 0 changed |
| Test suite | 880 pass, 0 fail (`pixi run test`) |
| axe-core, 8 pages | zero serious or critical |
| Playwright, 8 pages x 2 viewports | no console errors |
| First paint | 396 ms median against a 2 s budget |
| Energy drift | within noise at 3 sigma; the 2 fs fallback is unused |
| Receptor ingest path | exercised; 5 motifs, 116 references, 100/100 agreeing |
| Bilayer QC rows | exercised on a real membrane trajectory; G6 closed |

The test figure in this file has now been wrong twice, at 833 and at 837, both
times because I wrote it from memory or from arithmetic rather than from a run.
880 is measured, and it moves as tests are added. A bare `pytest` from the repository root also used to collect
the vendored BoltzMaker suite, whose 17 environment failures I briefly read as a
regression of this project; `pytest.ini` now scopes collection so that the plain
command and `pixi run test` agree.

**What unblocks it.** One of:

```
# everything, as 74 commits
git -C /Users/dellboy/Documents/Vibe_Coding/FlexAppeal push origin main
```

or a narrower answer if you would rather keep some of it local. Worth knowing
before you decide:

- The repository is **public**. Nothing in these commits contains a compound
  series, a SMILES string or a private structure; the systems are ABL1 with
  imatinib and β2AR with carazolol, both long-published, and every structure is
  fetched at run time rather than committed.
- `.gitignore` now names `DEEPFLEX_*.md` and `FXA2_BUILD_PROMPT.md`. Those were
  untracked and one `git add -A` from being published, which nearly happened
  earlier in this build.
- `DECISIONS.md` and `FINDINGS.md` describe several of my own mistakes in
  detail, including one case where I misread a comparison and spent four
  diagnostics on the wrong cause. That is deliberate, and it is public if this
  pushes.

**What happens next either way.** Tier 2 keeps running. The analysis, the
explorer pages and the findings all work from the local checkout.

---

## G3 Deploy — waiting, and not urgent

**What it blocks.** rsync to flexappeal.mdeller.com.

**What is done.** A fifth tab, three pages and four JSON endpoints, all reading
one SQLite file. The feature is additive: with no database present the section
renders empty rather than erroring, so deploying the code before the data is
safe.

**What unblocks it.** Approval to run the existing deploy path. Two things to
decide with it:

1. Whether to ship the explorer **now**, with the partial ABL1 run, or wait for
   Tier 2. The pages are honest about it: a run with no `run.json` is labelled
   `running` rather than shown as complete.
2. Whether `data/refstates.db` goes across. It is gitignored, so it travels by
   rsync or not at all. Without it the tab is an empty section.

**What happens next either way.** Nothing blocks. The pages are tested locally
at every route.

**The QC that was outstanding is now done.** Marc approved installing
Playwright, so section 10's checks ran properly rather than being approximated.
Eight pages at both viewports with no console errors, eight through axe-core
with zero serious or critical violations, 25 tests. It found eleven real
problems, one critical, including four contrast failures that all came from the
same mistake: `opacity` on small text, which makes the rendered colour something
no token states. All fixed. FINDINGS section 4b has the detail.

---

## G1 Compute — will not open

Recorded because section 0.2 asks for the projection rather than silence.

**Measured on the Tier 1 run, with the barostat on.** The β2AR system produces
at **60 ns/day**, holding 303.5 to 304.3 K against a 303.15 K target and 1.017
g/mL. That is below the 86 ns/day the heating stages suggested, and the
difference is the `MonteCarloMembraneBarostat`, which costs about 30 per cent.
The earlier figure in this file was an estimate taken before production started
and is replaced rather than kept.

**Corrected to include the stages before production**, which the first version
of this table left out. Section 5.2 puts 10 ns of NPT equilibration in front of
every membrane run and 1 ns in front of every water run, and that is real GPU
time: for the membrane rows it is a fifth of the row again. The gentle 0.1 ns of
heating runs at 1 fs, so it costs four times its length at the production step.

| row | ns/day | heat | equil | production | GPU-hours |
|---|---|---|---|---|---|
| abl1_apo | 213 measured | 0.5 | 1 | 50 | 5.3, with 4.3 ns already done |
| abl1_imatinib | 213 assumed | 0.5 | 1 | 50 | 5.8 |
| b2ar_apo | **59.8 measured** | 0.5 | 10 | 50 | **24.4** |
| b2ar_carazolol | ~61 assumed | 0.5 | 10 | 50 | ~23.9 |
| | | | | | **59.4** |

Against `MAX_COMPUTE_HOURS` of 96, with **36.6 hours spare**. That is still
enough to redo either ABL1 row twice over, or one membrane row with 12 hours to
spare, which is the margin that matters because D009 dropped the replicas and a
failed row now costs the whole row rather than a third of it.

The earlier figure in this file said 50 hours. It was wrong by the
equilibration, and is corrected rather than quietly adjusted, because the whole
point of this section is that the number is measured.

So **G1 does not open** and Tier 3 on RunPod is not needed. `deploy/runpod/`
still gets written and not executed, as section 5.3 requires.

## G2 Secrets — not needed

No credential is required for anything in Tier 2. G2 would only open if Tier 3
were wanted, and on the figures above it is not.
