# Tier 3 on RunPod — written, never run

Section 5.3 asks for this to exist and section 11 asks for it not to be
executed, so nothing here has been run and the numbers below are the reason it
has not needed to be.

## Why it is not needed

Tier 2 fits on the Mac Studio. Measured, not projected:

| row | atoms | ns/day | 50 ns |
|---|---|---|---|
| ABL1, water | 46,721 | 213 | 5.6 h |
| β2AR, POPC | 194,032 | 60 | 20.0 h |

Four rows come to about 50 GPU-hours against `MAX_COMPUTE_HOURS` of 96, so
**G1 does not open**. This directory is for the case where it would: the spec's
Tier 3 lengths are 500 ns for ABL1 and 1 µs for β2AR, which is 10 to 20 times
Tier 2 and does not fit on one consumer GPU.

Starting it needs **G1** (budget approval) and **G2** (a RunPod API key, passed
through the environment and never written into a file).

## What the Tier 3 case actually buys

Worth stating before anyone spends money on it. The β2AR active-to-inactive
relaxation this build is looking for was seen by Dror et al. on the **microsecond**
scale. Its absence at 50 or 100 ns is an expected limitation rather than a
negative result, which section 9 is explicit about. So Tier 3 is the only tier
that could turn that particular absence into an observation, and that is the
argument for it; nothing in the Tier 2 results can be fixed by running them
longer on faster hardware.

## Files

- `job.sh` — the whole run on a fresh pod: clone, restore the environment from
  the committed `pixi.lock`, then run the same `tools/queue_md.py` the Mac uses.
- `sync.sh` — pull checkpoints and trajectories back, and push them up again on
  resume.

## Resume

A pod can be reclaimed at any time, and a 1 µs run cannot start again from
zero. `tools/run_md.py` already checkpoints every nanosecond and resumes from
`state.chk`, which is what makes this safe; `job.sh` does nothing clever beyond
restoring the checkpoint directory before it starts and syncing it out
afterwards. The resume path is the same code the Mac runs, so it is exercised
every time a local run is interrupted rather than only in an emergency.
