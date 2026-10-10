#!/usr/bin/env bash
# Rebuild the worked example from scratch.
#
# Everything the Example tab shows -- the option tables, the terminal captures,
# every analysis panel and every binding-site panel -- comes from the artefacts
# this script produces. None of it is written by hand, so regenerating is the
# only way to update the page, and a change to the run script or the analysis is
# reflected the next time this is run.
#
# T4 lysozyme L99A with benzene in its buried cavity (PDB 181L). The ligand is
# the point: hen egg-white lysozyme does not bind benzene, and the apo example
# beside this one cannot exercise a single binding-site panel. This is the
# canonical cavity-binding model system, so the pocket the panels describe is one
# the literature already has an opinion about.
#
# Takes roughly eight hours on an Apple M2 Ultra: about seven of simulation at
# ~320 ns/day on OpenCL, the rest preparation, parameterisation and analysis.
#
# It needs the FlexAppeal development environment, through `pixi run` and NOT by
# calling .pixi/envs/default/bin/python directly. That distinction is the whole
# reason this script exists in this form: the ligand needs AM1-BCC charges, which
# means antechamber and sqm on PATH and AMBERHOME set, and only pixi's activation
# does that. Calling the interpreter directly fails two minutes in with
# "AmberTools is not on PATH", which is how the first run of this died.
#
# Usage:  bash examples/lysozyme_benzene_100ns/rebuild.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

BLUE=$'\033[38;2;30;115;190m'; GREEN=$'\033[38;2;0;208;132m'; RESET=$'\033[0m'
step() { printf "%s\n" "${BLUE}→${RESET} $1"; }
ok()   { printf "%s\n" "${GREEN}✓${RESET} $1"; }

command -v pixi >/dev/null 2>&1 || {
  echo "pixi is not installed -- see install.sh"; exit 1; }

JOB=lysozyme_benzene_100ns

step "building the bundle from $HERE/config.json"
# The chemical definition for benzene is fetched here, at prep time, and travels
# in the bundle: a PDB HETATM record has no bond orders, and the run script
# transfers them from that SDF onto the deposited coordinates.
pixi run python FlexAppeal.py bundle tests/fixtures/181l.pdb \
  -c "$HERE/config.json" \
  -o "$HERE/flexappeal_$JOB.command"
chmod +x "$HERE/flexappeal_$JOB.command"

step "extracting it"
rm -rf "$HERE/run"
pixi run python - "$HERE" "$JOB" <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, ".")
from flexappeal import bundle

here, job = Path(sys.argv[1]), sys.argv[2]
files = bundle.unpack((here / f"flexappeal_{job}.command").read_bytes())
out = here / "run"
for name, data in files.items():
    target = out / name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
print(f"  {len(files)} files")
PY

# script(1) rather than a plain redirect, on purpose: the run must believe it is
# talking to a terminal or rich prints milestone lines instead of the bars, and
# the bars are half of what the Example tab is showing. COLUMNS is pinned so the
# captures are reproducible rather than depending on the window that ran this.
step "running 100 ns -- roughly seven hours, progress on screen"
( cd "$HERE/run" && \
  script -q "$HERE/terminal.log" /bin/bash -c \
    "stty cols 120 rows 45 2>/dev/null; COLUMNS=120 pixi run --manifest-path '$ROOT/pixi.toml' python run.py" )

step "analysing"
( cd "$HERE/run" && \
  script -q "$HERE/analyse.log" /bin/bash -c \
    "stty cols 120 rows 45 2>/dev/null; COLUMNS=120 pixi run --manifest-path '$ROOT/pixi.toml' python analyse.py" )

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
