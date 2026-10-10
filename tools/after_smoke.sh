#!/bin/bash
# Wait for the Tier 1 smoke run to finish, then start Tier 2 only if it passed.
# Section 5.3 requires 1 ns of a row before any long run on it, and the point of
# chaining rather than launching both is that a failed smoke must not be
# followed by 24 hours of the same failure.
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal
SMOKE=data/md/b2ar_apo_smoke/run.json

while true; do
  [ -f "$SMOKE" ] && break
  pgrep -f "queue_md.py --smoke" >/dev/null || break
  sleep 30
done

if [ ! -f "$SMOKE" ]; then
  echo "[$(date +%H:%M:%S)] Tier 1 did not produce run.json; Tier 2 not started"
  exit 1
fi

echo "[$(date +%H:%M:%S)] Tier 1 passed:"
cat "$SMOKE"
echo "[$(date +%H:%M:%S)] starting Tier 2, all four rows at 50 ns"
exec caffeinate -dimsu pixi run python tools/queue_md.py
