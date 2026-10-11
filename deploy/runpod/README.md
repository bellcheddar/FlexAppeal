# Tier 3 on RunPod

Written during the build and **not executed**, which is what section 5.3 asks
for. Nothing here has ever run. Two gates stand in front of it:

- **G1 (compute).** Tier 3 is 500 ns per ABL1 row and 1 µs per β2AR row,
  against Tier 2's 50 ns. Twenty times the GPU hours is a budget decision, not
  an assumption, and section 2.2 wants it made on measured ns/day rather than
  on a figure scaled from the Mac's.
- **G2 (secrets).** A RunPod API key and whatever the sync remote needs.
  Section 0.2 is explicit that a credential arrives in the environment and is
  never written to a file, so `.env.example` names the variables and holds no
  values.

## The two scripts

    job.sh    on the pod.  One row, one length. Clones, installs from the
              lock, restores whatever a previous pod got through, runs the
              same queue the Mac runs, pushes the result back.
    sync.sh   on the pod.  One row's checkpoint and trajectory between the pod
              and durable storage, over rclone, so the same script works
              against S3, Backblaze or SFTP without knowing which.

    bash deploy/runpod/job.sh b2ar_apo 1000
    bash deploy/runpod/job.sh abl1_apo 500

There is no launcher here and that is deliberate. Creating the pod is a console
action behind G1 and G2, and a script that creates billable infrastructure is
the wrong thing to have lying about in a public repository next to one that
says it has never been run.

## Three things that are not obvious

**A checkpoint cannot cross the platform that wrote it.** OpenMM refuses one
written on OpenCL if the context is CUDA, with "Checkpoint was created with a
different Platform: OpenCL". So a part-finished Tier 2 row on the Mac cannot be
carried up here and continued, and Tier 3 rows start from their own
equilibration. Found on 11 October, section 4af of `FINDINGS.md`. Between pods
it is not a problem, because both are CUDA, which is what `sync.sh` relies on.

**The fallback ladder is wrong on a rented GPU.** `tools/run_md.py` tries CUDA,
then OpenCL, then CPU. That is right on a machine you own, where a bad GPU costs
time. Here it would drop a 1 µs row to the CPU at a fiftieth of the speed and go
on billing, with one line in a log nobody is reading to say so. `job.sh` passes
`--require-platform CUDA`, which turns the ladder off: a pod that cannot give it
one stops instead.

**A pod's disk does not survive the pod.** Everything that must outlive it goes
through `sync.sh`, and the checkpoint goes first: losing the trajectory costs
the frames since the last sync, losing the checkpoint costs the run. `job.sh`
pushes on both the success and the failure path, because a failed row's
checkpoint is the most valuable thing on the disk.

## Resuming

Run `job.sh` again with the same row and `RUNPOD_SYNC_URL`. It pulls the
checkpoint back, `tools/queue_md.py` skips a row whose `run.json` says it
reached the target, and `tools/run_md.py` resumes the rest. `job.sh` prints
`--preflight` first, which reports what each row would do and takes no GPU, so
the question of whether this resumes or starts again is answered before the
meter runs rather than after.

## The budget, when G1 is opened

Measure, do not assume. Run the Tier 1 smoke on the pod first:

    pixi run python tools/queue_md.py --smoke --require-platform CUDA

It reports ns/day per row and the projected hours for the production length,
which is the number G1 needs. Forty minutes of pod time is the cheapest way to
find out that the budget was wrong.
