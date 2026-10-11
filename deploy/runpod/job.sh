#!/usr/bin/env bash
# One Tier 3 run on a fresh RunPod pod. Never executed; see README.md.
#
#   bash deploy/runpod/job.sh <row-name> <production-ns>
#
# Expects in the environment, never in a file:
#   FLEXAPPEAL_REPO   git URL to clone
#   RUNPOD_SYNC_URL   where checkpoints and trajectories live between pods
#
# Deliberately not a Dockerfile. The environment is already pinned by the
# committed pixi.lock, and a second description of it in a Dockerfile is a
# second thing to keep in step with the first.
set -euo pipefail

ROW="${1:?a row name, for example b2ar_apo}"
NS="${2:?production nanoseconds}"
WORK="${WORK:-/workspace}"
REPO="${FLEXAPPEAL_REPO:?set FLEXAPPEAL_REPO}"

echo "[$(date -u +%H:%M:%S)] pod $(hostname), row $ROW, target ${NS} ns"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

if [[ ! -d "$WORK/FlexAppeal" ]]; then
    git clone --depth 1 "$REPO" "$WORK/FlexAppeal"
fi
cd "$WORK/FlexAppeal"

# From the lock, not from the manifest: a Tier 3 run that resolved different
# package versions from the Tier 2 runs would not be comparable with them, and
# that is the whole point of running it.
if ! command -v pixi >/dev/null; then
    curl -fsSL https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
fi
pixi install --locked

# Whatever a previous pod got through. Restored before the run so run_md.py
# finds state.chk and continues rather than starting again.
if [[ -n "${RUNPOD_SYNC_URL:-}" ]]; then
    bash deploy/runpod/sync.sh pull "$ROW" || echo "nothing to resume from"
fi

mkdir -p "data/md/$ROW"
if [[ -f "data/md/$ROW/state.chk" ]]; then
    echo "[$(date -u +%H:%M:%S)] resuming from a checkpoint"
fi

# What the row will actually do, before anything takes the GPU. Reads only, so
# it costs nothing, and on a pod billed by the hour the answer to "will this
# resume or start again" is worth having before the meter runs rather than
# after.
pixi run python tools/queue_md.py --preflight --production-ns "$NS" || true

# The same queue the Mac runs, so there is one definition of the protocol. The
# membrane rows need --membrane and their prebuilt system, which queue_md.py
# already knows; this only passes the row and the length.
# --require-platform CUDA, because the fallback ladder is wrong here. On a
# machine you own, falling back to the CPU costs time. On a pod billed by the
# hour it drops a microsecond row to a fiftieth of the speed and goes on
# billing, with one line in a log nobody is reading to say so.
#
# `|| status=$?` rather than a bare call followed by `status=$?`. Under the
# `set -e` at the top of this file, a failing command exits the script
# immediately, so the assignment was unreachable and the push below never ran:
# a row that crashed, or a pod reclaimed mid-run, lost the checkpoint it had
# just spent hours earning. Making the call part of a condition is what stops
# `set -e` firing, and the push now happens on both paths.
status=0
pixi run python tools/queue_md.py --only "$ROW" --production-ns "$NS" \
    --require-platform CUDA || status=$?

if [[ -n "${RUNPOD_SYNC_URL:-}" ]]; then
    # Deliberately not guarded by $status. A failed row's checkpoint is the
    # most valuable thing on the disk, because the next pod resumes from it.
    bash deploy/runpod/sync.sh push "$ROW" || echo "push failed; the checkpoint is still on this pod"
fi

echo "[$(date -u +%H:%M:%S)] exit $status"
exit $status
