#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation}"
COLLECTOR="$ANALYSIS_ROOT/scripts/update_hourly_report.py"
SOURCE_MONITOR="$ANALYSIS_ROOT/monitoring"
PYTHON_BIN="${PYTHON_BIN:-/home/petadimensionlab/workspace/Geneformer/.venv/bin/python}"
LOCK_FILE="${TMPDIR:-/tmp}/geneformer-nsclc-monitor-refresh.lock"

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  echo "A live-report refresh is already running; skipping this interval."
  exit 0
fi

if ! tmux has-session -t gf_heldout_allgene 2>/dev/null \
  && [[ -f "$SOURCE_MONITOR/latest_status.json" ]] \
  && grep -q '"run_active": false' "$SOURCE_MONITOR/latest_status.json"; then
  echo "The monitored run is already complete; no additional idle snapshot is needed."
  exit 0
fi

if [[ ! -f "$COLLECTOR" ]]; then
  echo "Missing live status collector: $COLLECTOR" >&2
  exit 1
fi

"$PYTHON_BIN" "$COLLECTOR"
cp "$SOURCE_MONITOR/latest_status.json" "$DIR/latest_status.json"
cp "$SOURCE_MONITOR/hourly_history.csv" "$DIR/hourly_history.csv"
cp "$SOURCE_MONITOR/gpu_statistics.png" "$DIR/gpu_statistics.png"
sed -i 's/\r$//' "$DIR/hourly_history.csv"
"$DIR/report_generation_job.sh"
