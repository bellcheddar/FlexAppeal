#!/usr/bin/env bash
# Move one row's checkpoint and trajectory between a pod and durable storage.
# Never executed; see README.md.
#
#   bash deploy/runpod/sync.sh pull|push <row-name>
#
# RUNPOD_SYNC_URL is an rclone remote, so the same script works against S3,
# Backblaze or a plain SFTP host without this file knowing which.
set -euo pipefail

ACTION="${1:?pull or push}"
ROW="${2:?a row name}"
REMOTE="${RUNPOD_SYNC_URL:?set RUNPOD_SYNC_URL}"
LOCAL="data/md/$ROW"

# The checkpoint is what makes a reclaimed pod survivable, so it goes first and
# the trajectory second: losing the trajectory costs the frames since the last
# sync, losing the checkpoint costs the run.
case "$ACTION" in
  pull)
      mkdir -p "$LOCAL"
      rclone copy "$REMOTE/$ROW" "$LOCAL" \
             --include "state.chk" --include "system.xml" \
             --include "solvated.pdb" --include "trajectory.dcd" \
             --include "state_data.csv" -P
      ;;
  push)
      rclone copy "$LOCAL" "$REMOTE/$ROW" \
             --include "state.chk" -P
      rclone copy "$LOCAL" "$REMOTE/$ROW" \
             --include "trajectory.dcd" --include "state_data.csv" \
             --include "run.json" -P
      ;;
  *)  echo "pull or push"; exit 2 ;;
esac
