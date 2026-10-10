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
done_ns=0
if [[ -f "data/md/$ROW/state.chk" ]]; then
    echo "[$(date -u +%H:%M:%S)] resuming from a checkpoint"
fi

# The same queue the Mac runs, so there is one definition of the protocol. The
# membrane rows need --membrane and their prebuilt system, which queue_md.py
# already knows; this only passes the row and the length.
pixi run python tools/queue_md.py --only "$ROW" --production-ns "$NS"
status=$?

if [[ -n "${RUNPOD_SYNC_URL:-}" ]]; then
    bash deploy/runpod/sync.sh push "$ROW"
fi

echo "[$(date -u +%H:%M:%S)] exit $status"
exit $status
