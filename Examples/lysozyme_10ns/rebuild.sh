#!/usr/bin/env bash
# Rebuild the worked example from scratch.
#
# Everything the Example tab shows -- the option tables, the terminal captures,
# every analysis panel -- comes from the artefacts this script produces. None of
# it is written by hand, so regenerating is the only way to update the page, and
# a change to the run script or the analysis is reflected the next time this is
# run.
#
# Takes roughly an hour and a half on an Apple M1 Max: about seventy minutes of
# simulation at ~200 ns/day on OpenCL, the rest preparation and analysis. It
# needs the FlexAppeal development environment (pixi install), not the bundle's
# own -- the bundle would solve a second multi-gigabyte environment to do
# exactly the same work.
#
# Usage:  bash Examples/lysozyme_10ns/rebuild.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

BLUE=$'\033[38;2;30;115;190m'; GREEN=$'\033[38;2;0;208;132m'; RESET=$'\033[0m'
step() { printf "%s\n" "${BLUE}→${RESET} $1"; }
ok()   { printf "%s\n" "${GREEN}✓${RESET} $1"; }

command -v pixi >/dev/null 2>&1 || {
  echo "pixi is not installed -- see install.sh"; exit 1; }

step "building the bundle from Examples/lysozyme_10ns/config.json"
pixi run python FlexAppeal.py bundle tests/fixtures/1aki.pdb \
  -c "$HERE/config.json" \
  -o "$HERE/flexappeal_lysozyme_10ns.command"
chmod +x "$HERE/flexappeal_lysozyme_10ns.command"

step "extracting it"
rm -rf "$HERE/run"
pixi run python - "$HERE" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, ".")
from flexappeal import bundle

here = Path(sys.argv[1])
files = bundle.unpack((here / "flexappeal_lysozyme_10ns.command").read_bytes())
out = here / "run"
out.mkdir(parents=True, exist_ok=True)
for name, data in files.items():
    (out / name).write_bytes(data)
print(f"  {len(files)} files")
PY

# script(1) rather than a plain redirect, on purpose: the run must believe it is
# talking to a terminal or rich prints milestone lines instead of the bars, and
# the bars are half of what the Example tab is showing. COLUMNS is pinned so the
# captures are reproducible rather than depending on the window that ran this.
step "running 10 ns -- roughly seventy minutes, progress on screen"
( cd "$HERE/run" && \
  script -q "$HERE/terminal.log" /bin/bash -c \
    "stty cols 120 rows 45 2>/dev/null; COLUMNS=120 '$ROOT/.pixi/envs/default/bin/python' run.py" )

# Captured the same way as the run, and for the same reason: analyse.py draws a
# progress bar over the metric list, and without a pty rich prints milestone
# lines instead. Recording it also means a rebuild reproduces every artefact
# rather than most of them.
step "analysing"
( cd "$HERE/run" && \
  script -q "$HERE/analyse.log" /bin/bash -c \
    "stty cols 120 rows 45 2>/dev/null; COLUMNS=120 '$ROOT/.pixi/envs/default/bin/python' analyse.py" )

step "collecting the results file"
mkdir -p "$HERE/output"
find "$HERE/run" -name '*.fxa' -exec cp {} "$HERE/output/" \;

step "rendering the terminal captures"
pixi run python scripts/terminal_capture.py "$HERE/terminal.log" \
  --out "$HERE/screenshots"

step "compressing the raw captures"
gzip -9 -f -k "$HERE/terminal.log" "$HERE/analyse.log"

ok "done. The Example tab will pick this up on its next request."
printf "   %s\n" "$(du -sh "$HERE/output" | cut -f1) of results, $(ls "$HERE/screenshots" | wc -l | tr -d ' ') captures"
