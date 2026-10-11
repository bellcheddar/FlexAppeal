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
#
# **Why the process check needs a minute of patience.** It looked once. The
# queue runs each row as a fresh subprocess, so at every row boundary there is a
# moment with no runner alive and nothing yet saying the queue is done, and one
# look in that window would have ended the streaming for every row after the
# first and fired the final pass on a partial trajectory. The gap is seconds
# against a half-hour interval, so it would have happened rarely and silently.
set -u
# The repository root, overridable. It was a bare absolute path, which is fine
# for the launcher that runs this and wrong for anything else: a test that
# sourced this to exercise one of its helpers was silently moved to the real
# repository and read the real `data/md`, so its own fixtures were never
# consulted and it passed whatever the helper did. Defaulting to the script's
# own location keeps the launcher working from anywhere and lets a test point
# it at a directory it controls.
cd "${FLEXAPPEAL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

INTERVAL="${INTERVAL:-1800}"
QUEUE_LOG="${QUEUE_LOG:-data/md/production.log}"
# A minute of grace before an absent production process counts as the end, as
# twelve looks five seconds apart. Overridable so the test for it does not take
# a minute to decide that nothing is there.
GRACE_CHECKS="${GRACE_CHECKS:-12}"
GRACE_SLEEP="${GRACE_SLEEP:-5}"

# This script's own process tree, which matches any pattern naming run_md.py
# because the pattern is written above.
ancestors() {
    local pid=$$
    while [ -n "$pid" ] && [ "$pid" -gt 1 ]; do
        printf '%s\n' "$pid"
        pid=$(ps -o ppid= -p "$pid" 2>/dev/null | tr -d ' ')
    done
}

# A production process that is not one of our own: either the queue, which
# lives across row boundaries, or a runner. Both patterns, because the queue is
# the process whose lifetime actually matches production and the runner catches
# a queue started by hand for one row. Excluded by PID rather than by pattern
# for the reason above: this script's own text contains both names.
#
# A `--preflight` run would match too. It cannot matter: it is read-only and
# lasts under a second, and by the time the queue is genuinely finished the
# "queue done" marker short-circuits this check before it is reached.
production_alive() {
    local mine
    mine=$(ancestors | tr '\n' '|' | sed 's/|$//')
    pgrep -f "tools/(queue|run)_md\.py" 2>/dev/null \
        | grep -qvE "^(${mine})$"
}

producing() {
    # The queue says so itself once every row is done or has failed.
    if [ -f "$QUEUE_LOG" ] && grep -q "queue done" "$QUEUE_LOG"; then
        return 1
    fi
    # An absence has to persist before it counts as the end of production. The
    # queue starts each row as a fresh subprocess, so there is a moment at
    # every row boundary with no runner alive and no "queue done" written, and
    # concluding from a single look would end the streaming for every row after
    # the first. The symptom would be the bad kind: the final pass, which
    # computes the bilayer QC and G6, would run on a third of a trajectory and
    # nothing would re-run it.
    #
    # The same window opens when the queue is restarted by hand. That happened
    # at 04:19 on 11 October and lasted about a minute; the loop survived it
    # only because it was asleep. Hence a minute of grace, which costs one
    # minute of delay on the one occasion production has really stopped.
    local attempt
    for attempt in $(seq "$GRACE_CHECKS"); do
        if production_alive; then
            return 0
        fi
        sleep "$GRACE_SLEEP"
    done
    return 1
}

# Sourced rather than run, the functions above are available and nothing
# happens. That is how `producing` is tested: a test sources this with a `pgrep`
# stub on PATH and checks that a gap with no production process in it does not
# end the loop. Without the guard, sourcing it would start a half-hour loop
# inside the test suite.
if [ "${BASH_SOURCE[0]}" != "${0}" ]; then
    return 0 2>/dev/null || true
fi

main() {
    while producing; do
        echo "[$(date +%H:%M:%S)] streaming pass"
        bash tools/analyse_all.sh --skip-membrane 2>&1 | grep -E "analysing|ingesting|motifs|rows|FAILED|no trajectory" || true
        sleep "$INTERVAL"
    done

    echo "[$(date +%H:%M:%S)] production has stopped; one final pass with the bilayer QC"
    bash tools/analyse_all.sh

    # Section 2.2 wants the measured ns/day recorded rather than assumed, and a
    # row's real rate is only known once it has finished, so this refreshes as the
    # queue drains rather than being a one-off at the start.
    echo "[$(date +%H:%M:%S)] refreshing config/benchmarks.toml"
    pixi run python tools/hwprobe.py || echo "hwprobe failed, see above"

    # The figures after the analysis, because they are drawn from the database the
    # pass above has just finished filling. Section 11 wants FINDINGS "with
    # figures", and a figure drawn from a half-written run is worse than no figure,
    # so this is the one place they are regenerated rather than on every pass.
    echo "[$(date +%H:%M:%S)] redrawing the figures"
    pixi run python tools/make_figures.py || echo "figures failed, see above"

    # Every stored number against the analysis file it came from. The rows land
    # unattended, and a fault in the ingest produces numbers that are wrong rather
    # than missing, which nothing downstream would notice. Reported rather than
    # fatal: a disagreement is worth finding in the morning, and stopping the rest
    # of the pass over it would lose the controls below as well.
    echo "[$(date +%H:%M:%S)] verifying the stored numbers against the analyses"
    VERIFY_LOG="data/md/verify.log"
    # Both of the verifier's checks are echoed, not just the first. They catch
    # different things: the round-trip proves the stored numbers match the
    # analysis they came from, and the metric check proves there was something
    # to store, which an all-NaN metric passes the round trip without having.
    # A check whose result never reaches this log is close to one that did not
    # run, since this log is the only account of the pass anyone reads.
    if pixi run python -m flexappeal.refstates.cli report --verify \
          > "$VERIFY_LOG" 2>&1; then
        grep -E "agrees with the analysis|has a reading" "$VERIFY_LOG" || true
    else
        echo "VERIFICATION FAILED; the first lines of $VERIFY_LOG:"
        grep -E "disagree|no reading at all|^    [a-z0-9_]+ " "$VERIFY_LOG" \
            | head -20
    fi

    # And the controls, which is what says whether P5 and P6 came out. They read
    # the same database and cost seconds, so there is no reason to make anyone
    # remember to run them.
    echo "[$(date +%H:%M:%S)] positive controls"
    pixi run python -m flexappeal.refstates.cli controls 2>&1 \
      | sed -n '/^P[0-9]:/,$p' || true

    echo "[$(date +%H:%M:%S)] final pass done"
}

main "$@"
