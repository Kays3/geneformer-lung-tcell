#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation}"
COLLECTOR="$ANALYSIS_ROOT/scripts/update_hourly_report.py"
SOURCE_MONITOR="$ANALYSIS_ROOT/monitoring"
PYTHON_BIN="${PYTHON_BIN:-/home/petadimensionlab/workspace/Geneformer/.venv/bin/python}"
LOCK_FILE="${TMPDIR:-/tmp}/geneformer-nsclc-monitor-refresh.lock"
PUBLISH_TO_GIT="${PUBLISH_TO_GIT:-0}"

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

if [[ "$PUBLISH_TO_GIT" == "1" ]]; then
  CHECKOUT="$(git -C "$DIR" rev-parse --show-toplevel)"
  BRANCH="$(git -C "$CHECKOUT" symbolic-ref --quiet --short HEAD)"
  MONITOR_PATH="current_workflow/monitoring"
  GENERATED_PATHS=(
    "$MONITOR_PATH/GPU_PROGRESS_REPORT.md"
    "$MONITOR_PATH/latest_status.json"
    "$MONITOR_PATH/hourly_history.csv"
    "$MONITOR_PATH/progress_animation.gif"
    "$MONITOR_PATH/progress_animation.svg"
    "$MONITOR_PATH/gpu_statistics.png"
    "$MONITOR_PATH/cell_interaction_diagram.svg"
    "$MONITOR_PATH/snapshot_gallery"
  )

  git -C "$CHECKOUT" add -- "${GENERATED_PATHS[@]}"
  if git -C "$CHECKOUT" diff --cached --quiet -- "${GENERATED_PATHS[@]}"; then
    echo "No generated monitoring changes to publish."
    exit 0
  fi

  STAMP="$(date '+%Y-%m-%d %H:%M %Z')"
  git -C "$CHECKOUT" commit --only -m "monitor: update perturbation status $STAMP" -- "${GENERATED_PATHS[@]}"
  git -C "$CHECKOUT" push origin "HEAD:$BRANCH"
fi
