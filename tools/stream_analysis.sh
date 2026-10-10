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
set -u
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal

INTERVAL="${INTERVAL:-1800}"

while pgrep -f "tools/run_md.py" >/dev/null; do
    echo "[$(date +%H:%M:%S)] streaming pass"
    bash tools/analyse_all.sh --skip-membrane 2>&1 | grep -E "analysing|ingesting|motifs|rows|FAILED|no trajectory" || true
    sleep "$INTERVAL"
done

echo "[$(date +%H:%M:%S)] production has stopped; one final pass with the bilayer QC"
bash tools/analyse_all.sh
