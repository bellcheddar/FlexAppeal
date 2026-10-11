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
# The repository root, overridable. It was a bare absolute path, which is fine
# for the launcher that runs this and wrong for anything else: a test that
# sourced this to exercise one of its helpers was silently moved to the real
# repository and read the real `data/md`, so its own fixtures were never
# consulted and it passed whatever the helper did. Defaulting to the script's
# own location keeps the launcher working from anywhere and lets a test point
# it at a directory it controls.
cd "${FLEXAPPEAL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# --skip-membrane for the streaming passes. The state analysis loads protein and
# ligand only, which is a few megabytes, while the bilayer QC needs the full
# solvated trajectory and that is most of a gigabyte once a 50 ns membrane row
# is under way. Cheap enough to run every half hour one way, not the other.
#
# It skips the bilayer QC for a row that is **still producing**, which is what
# the expense is actually about, rather than for every row in the pass. A row
# whose run.json says complete will not change again, so its QC is computed as
# soon as it finishes and the row is final from that moment.
#
# Without the distinction a membrane row that finished early in a queue waited
# for the whole queue to drain before anyone could see its bilayer. b2ar_apo
# finished at 17:36 on 11 October with three days of queue behind it, and its
# QC had to be run by hand to avoid exactly that wait. b2ar_carazolol finishes
# on the Monday with a row still to go, so it would have been the same again.
STREAMING=""
if [ "${1:-}" = "--skip-membrane" ]; then STREAMING="yes"; fi

# A row is finished when its own run.json says so. The queue writes that file
# at the end of a row, so its presence with a complete status is the same fact
# the queue itself acts on.
row_finished() {
  [ -f "data/md/$1/run.json" ] \
    && grep -q '"status"[[:space:]]*:[[:space:]]*"complete"' "data/md/$1/run.json"
}

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

# Sourced rather than run, the helpers above are available and nothing happens,
# which is how `row_finished` is tested. Without it, sourcing this would start a
# full analysis pass inside the suite. `stream_analysis.sh` carries the same
# guard for the same reason.
if [ "${BASH_SOURCE[0]}" != "${0}" ]; then
    return 0 2>/dev/null || true
fi

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

# Two ingests rather than one, because the bilayer QC is per row and so is the
# reason to defer it. A pass with nothing in one group simply does not run it.
final=()
producing=()
for run in "${ingested[@]}"; do
  if [ -z "$STREAMING" ] || row_finished "$run"; then
    final+=("$run")
  else
    producing+=("$run")
  fi
done

if [ ${#final[@]} -gt 0 ]; then
  names=()
  for run in "${final[@]}"; do names+=(--name "$run"); done
  echo "[$(date +%H:%M:%S)] ingesting ${final[*]} with the bilayer QC"
  pixi run python tools/ingest.py "${names[@]}"
fi

if [ ${#producing[@]} -gt 0 ]; then
  names=()
  for run in "${producing[@]}"; do names+=(--name "$run"); done
  echo "[$(date +%H:%M:%S)] ingesting ${producing[*]}, still producing so no bilayer QC"
  pixi run python tools/ingest.py "${names[@]}" --skip-membrane
fi

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
