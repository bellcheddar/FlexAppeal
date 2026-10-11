#!/usr/bin/env bash
# Re-analyse and re-ingest every row with frames, every half hour, for as long
# as production is running. Section 2.3's streaming analysis.
#
#   bash tools/stream_analysis.sh
#
# A loop rather than a daemon or a cron entry: it is tied to the lifetime of the
# thing it is watching, so it stops on its own when the queue does and leaves
# nothing behind to remember to turn off. The ingest is idempotent and mdtraj
# reads whatever frames are on disk, so each pass simply sees more of the run.
#
# **Why the liveness check is not just `pgrep -f tools/run_md.py`.** It was, and
# that guard could never go false. This script's own text contains the pattern,
# so the shell that launched it carried `tools/run_md.py` in its command line
# and matched. That parent lives exactly as long as the loop does, so the loop
# kept itself alive: the condition was true because the loop was running, which
# is a check that measures nothing. The visible symptom would have been the
# worst kind, because the loop would have gone on streaming perfectly well and
# only the final pass after it, the one that computes the bilayer QC and G6,
# would never have run.
#
# So two conditions, both on state rather than on a pattern that matches us:
# the queue's own "queue done" marker, and a process check that excludes this
# script's ancestors by PID.
set -u
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal

INTERVAL="${INTERVAL:-1800}"
QUEUE_LOG="${QUEUE_LOG:-data/md/production.log}"

# This script's own process tree, which matches any pattern naming run_md.py
# because the pattern is written above.
ancestors() {
    local pid=$$
    while [ -n "$pid" ] && [ "$pid" -gt 1 ]; do
        printf '%s\n' "$pid"
        pid=$(ps -o ppid= -p "$pid" 2>/dev/null | tr -d ' ')
    done
}

producing() {
    # The queue says so itself once every row is done or has failed.
    if [ -f "$QUEUE_LOG" ] && grep -q "queue done" "$QUEUE_LOG"; then
        return 1
    fi
    # Otherwise, a real run_md.py that is not one of our own processes. This
    # catches a queue that died without writing its marker.
    local mine
    mine=$(ancestors | tr '\n' '|' | sed 's/|$//')
    pgrep -f "tools/run_md.py" 2>/dev/null \
        | grep -qvE "^(${mine})$"
}

while producing; do
    echo "[$(date +%H:%M:%S)] streaming pass"
    bash tools/analyse_all.sh --skip-membrane 2>&1 | grep -E "analysing|ingesting|motifs|rows|FAILED|no trajectory" || true
    sleep "$INTERVAL"
done

echo "[$(date +%H:%M:%S)] production has stopped; one final pass with the bilayer QC"
bash tools/analyse_all.sh

# The figures last, because they are drawn from the database the pass above has
# just finished filling. Section 11 wants FINDINGS "with figures" and a figure
# drawn from a half-written run is worse than no figure, so this is the one
# place they are regenerated rather than on every streaming pass.
# Section 2.2 wants the measured ns/day recorded rather than assumed, and a
# row's real rate is only known once it has finished, so this refreshes as the
# queue drains rather than being a one-off at the start.
echo "[$(date +%H:%M:%S)] refreshing config/benchmarks.toml"
pixi run python tools/hwprobe.py || echo "hwprobe failed, see above"

echo "[$(date +%H:%M:%S)] redrawing the figures"
pixi run python tools/make_figures.py || echo "figures failed, see above"

# And the controls, which is what says whether P5 and P6 came out. They read
# the same database and cost seconds, so there is no reason to make anyone
# remember to run them.
# Every stored number against the analysis file it came from. The four rows
# land unattended, and a fault in the ingest produces numbers that are wrong
# rather than missing, which nothing downstream would notice. Reported rather
# than fatal: a disagreement is worth finding in the morning, and stopping the
# rest of the pass over it would lose the figures and the controls too.
echo "[$(date +%H:%M:%S)] verifying the stored numbers against the analyses"
pixi run python -m flexappeal.refstates.cli report --verify \
  | grep -E "agrees|disagree|^    [a-z0-9_]+ [a-zA-Z0-9_.]+:" \
  || echo "verification reported a disagreement, see above"

echo "[$(date +%H:%M:%S)] positive controls"
pixi run python -m flexappeal.refstates.cli controls 2>&1   | sed -n '/^P[0-9]:/,$p' || true

echo "[$(date +%H:%M:%S)] final pass done"
