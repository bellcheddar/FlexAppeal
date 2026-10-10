#!/bin/bash
# Analyse and ingest every row that has a trajectory, then report the database.
#
# Safe to run while production is still going: mdtraj reads whatever frames are
# on disk and the ingest is idempotent, so this can be run repeatedly and each
# pass simply sees more of the run. That is the streaming analysis section 2.3
# asks for, done as a cron-able script rather than as a daemon, because a daemon
# that has to be kept alive across a machine sleeping is more moving parts than
# the job needs.
set -u
cd /Users/dellboy/Documents/Vibe_Coding/FlexAppeal

# --skip-membrane for the streaming passes. The state analysis loads protein and
# ligand only, which is a few megabytes, while the bilayer QC needs the full
# solvated trajectory and that is most of a gigabyte once a 50 ns membrane row
# is under way. Cheap enough to run every half hour one way, not the other.
SKIP_MEMBRANE=""
if [ "${1:-}" = "--skip-membrane" ]; then SKIP_MEMBRANE="--skip-membrane"; fi

# A case statement rather than associative arrays: macOS ships bash 3.2, which
# has none, and `declare -A` there silently creates a plain array so the first
# lookup fails under `set -u`.
row_kind() {
  case "$1" in
    abl1_apo|abl1_imatinib) echo kinase ;;
    b2ar_apo|b2ar_carazolol) echo gpcr ;;
  esac
}
row_start() {
  case "$1" in
    abl1_apo) echo abl1_apo_start.pdb ;;
    abl1_imatinib) echo abl1_imatinib_start.pdb ;;
    b2ar_apo) echo b2ar_apo_start.pdb ;;
    b2ar_carazolol) echo b2ar_carazolol_start.pdb ;;
  esac
}
# The reference entry whose KLIFS pocket map resolves the construct's own motif
# residues. Only the kinase rows need one; the receptor resolves through GPCRdb.
row_refpdb() {
  case "$1" in
    abl1_apo) echo 2GQG ;;
    abl1_imatinib) echo 2HYY ;;
    *) echo "" ;;
  esac
}

ingested=()
for run in abl1_apo abl1_imatinib b2ar_apo b2ar_carazolol; do
  traj="data/md/$run/trajectory.dcd"
  if [ ! -f "$traj" ]; then
    echo "[$(date +%H:%M:%S)] $run: no trajectory yet, skipping"
    continue
  fi
  echo "[$(date +%H:%M:%S)] $run: analysing"
  args=(--name "$run" --start "data/systems/$(row_start "$run")" \
        --kind "$(row_kind "$run")")
  refpdb=$(row_refpdb "$run")
  [ -n "$refpdb" ] && args+=(--reference-pdb "$refpdb")
  if pixi run python tools/analyse_run.py "${args[@]}" > "data/md/$run.analysis.log" 2>&1; then
    ingested+=("$run")
  else
    echo "[$(date +%H:%M:%S)] $run: analysis FAILED, see data/md/$run.analysis.log"
    tail -5 "data/md/$run.analysis.log"
  fi
done

if [ ${#ingested[@]} -eq 0 ]; then
  echo "[$(date +%H:%M:%S)] nothing analysed"
  exit 1
fi

names=()
for run in "${ingested[@]}"; do names+=(--name "$run"); done
echo "[$(date +%H:%M:%S)] ingesting ${ingested[*]}"
pixi run python tools/ingest.py "${names[@]}" $SKIP_MEMBRANE

pixi run python - <<'PY'
import sqlite3
from pathlib import Path
db = Path("data/refstates.db")
if not db.is_file():
    raise SystemExit("no database written")
c = sqlite3.connect(db)
c.row_factory = sqlite3.Row
print("\nruns in the explorer database:")
for row in c.execute("SELECT run_id, start_state, length_ns, status, ns_per_day "
                     "FROM runs ORDER BY run_id"):
    speed = f"{row['ns_per_day']:.0f} ns/day" if row["ns_per_day"] else "speed not recorded"
    length = f"{row['length_ns']:.1f} ns" if row["length_ns"] else "no frames"
    print(f"  {row['run_id']:16s} {str(row['start_state']):10s} {length:>10s}  "
          f"{row['status']:9s} {speed}")
for table in ("motif_frames", "membrane", "reference_structures"):
    n = c.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    print(f"  {table:22s} {n:>8,} rows")
PY
