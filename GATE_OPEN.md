# Gates waiting on Marc

Written per section 0.2's gate protocol. Everything these gates do not block is
still running: Tier 1 is finishing, Tier 2 is chained behind it, and nothing
below needs an answer before those complete.

---

## G4 Publish — waiting

**What it blocks.** Any `git push`. 23 commits sit on `main` locally and none
has left this machine.

**What is done.** Phases 0 to 4 and 6 are complete, Phase 5 is running, Phase 7
is drafted. The spec's own acceptance criteria are met for every finished phase:

| Control | Result |
|---|---|
| P1 DFG label against KLIFS | 17 / 17 (spec wants ≥ 90 %) |
| P2 αC label against KLIFS | 17 / 17 (≥ 85 %) |
| P3 TM3–TM6 state separation, leave-one-out | 100 / 100 at 11.39 Å (≥ 90 %) |
| P4 each start scores as its own state | 4 / 4 |
| Test suite | 789 pass, 0 fail |

**What unblocks it.** One of:

```
# everything, as 23 commits
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

**One QC item is outstanding and needs a small decision with this gate.**
Section 10 asks for Playwright at 1440×900 and 390×844 and axe-core with no
serious violations. Neither can run here: `playwright` is not in the
environment, there is no preinstalled Chromium, and section 10 says never to run
`playwright install`. What ran instead is server-side: every `/static/` asset
referenced by any of these pages is fetched and has to return 200, and the
subset of axe-core's serious rules that live in the markup rather than in the
rendered page. That misses anything only visible at runtime, a JavaScript
exception, contrast, focus order. Unblocking it needs either a Chromium already
on the machine or approval to download one.

---

## G1 Compute — will not open

Recorded because section 0.2 asks for the projection rather than silence.

**Measured on the Tier 1 run, with the barostat on.** The β2AR system produces
at **60 ns/day**, holding 303.5 to 304.3 K against a 303.15 K target and 1.017
g/mL. That is below the 86 ns/day the heating stages suggested, and the
difference is the `MonteCarloMembraneBarostat`, which costs about 30 per cent.
The earlier figure in this file was an estimate taken before production started
and is replaced rather than kept.

| row | atoms | ns/day | 50 ns costs |
|---|---|---|---|
| abl1_apo | 46,721 | 213 measured | 5.6 h, of which 4.3 ns is already done |
| abl1_imatinib | ~46,700 | 213 assumed | 5.6 h |
| b2ar_apo | 194,032 | **60 measured** | **20.0 h** |
| b2ar_carazolol | 188,364 | 60 assumed, slightly fewer atoms | ~19.5 h |
| | | | **≈ 50 h** |

Against `MAX_COMPUTE_HOURS` of 96, with 46 hours spare. One membrane row could
be redone and the build would still come in at 70 hours, which is the margin
that matters because D009 dropped the replicas and a failed row now costs the
whole row rather than a third of it.

So **G1 does not open** and Tier 3 on RunPod is not needed. `deploy/runpod/`
still gets written and not executed, as section 5.3 requires.

## G2 Secrets — not needed

No credential is required for anything in Tier 2. G2 would only open if Tier 3
were wanted, and on the figures above it is not.
