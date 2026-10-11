#!/usr/bin/env bash
# Emit one line and exit when the queue reaches a terminal state, or when the
# current row stops producing frames.
#
# The existing watchers cover completion and failure, both of which the queue
# writes to its own log. Neither covers a hang: a run_md that stops stepping
# without raising leaves every process alive and every log quiet, and silence
# there is indistinguishable from a slow row. Over a 19-hour row that is the
# failure most worth catching, because it costs the whole row.
#
# Stall is defined on the state data rather than on the trajectory file, whose
# size can sit unchanged while a write is buffered.
set -u
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal

INTERVAL="${INTERVAL:-120}"
# 20 minutes of no new frame. A membrane row writes one every ~2.4 minutes at
# 59 ns/day, so this is eight missed frames and not a slow disk.
PATIENCE="${PATIENCE:-10}"

last=""
quiet=0
while true; do
    if grep -qE "queue done" data/md/production.log 2>/dev/null; then
        echo "QUEUE DONE: $(grep -E 'finished in|FAILED' data/md/production.log | tail -4 | tr '\n' ' ')"
        exit 0
    fi
    if grep -qE "FAILED" data/md/production.log 2>/dev/null; then
        echo "ROW FAILED: $(grep -E 'FAILED' data/md/production.log | tail -2 | tr '\n' ' ')"
        exit 0
    fi

    # The row the queue says it is on, so this follows the queue rather than
    # naming a row and going stale when the queue moves on.
    row=$(grep -E "^\[.*\] starting " data/md/production.log 2>/dev/null \
          | tail -1 | sed -E 's/.*starting ([a-z0-9_]+):.*/\1/')
    [ -z "$row" ] && { sleep "$INTERVAL"; continue; }

    if grep -qE "finished in" data/md/production.log 2>/dev/null; then
        echo "ROW FINISHED: $(grep -E 'finished in' data/md/production.log | tail -1)"
        exit 0
    fi

    current=$(tail -1 "data/md/$row/state_data.csv" 2>/dev/null || echo "")
    if [ -n "$current" ] && [ "$current" = "$last" ]; then
        quiet=$((quiet + 1))
        if [ "$quiet" -ge "$PATIENCE" ]; then
            echo "STALLED: $row has written no frame for $((quiet * INTERVAL / 60)) minutes; last state line: $current"
            exit 0
        fi
    else
        quiet=0
    fi
    last="$current"
    sleep "$INTERVAL"
done
