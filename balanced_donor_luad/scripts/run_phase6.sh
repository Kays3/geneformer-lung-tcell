#!/usr/bin/env bash
# Phase 6 driver (one host): ISP for this host's genes (host map), BOTH operations, all donors; GPU sampled every 60 s.
# Usage: run_phase6.sh <host: thinkstation1|thinkstation2> <repo_worktree> <data_root> <geneformer_root> <venv_python>
set -euo pipefail
HOST="$1"; REPO="$2"; DATA="$3"; GF="$4"; PY="$5"
case "$HOST" in thinkstation1|thinkstation2) ;; *) echo "host must match the host map: $HOST"; exit 2;; esac
BD="$REPO/balanced_donor_luad"; WORK="$DATA/phase6"; mkdir -p "$WORK"
nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$WORK/nvidia_smi_samples_$HOST.csv" 2>&1 &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_start_$HOST.csv"
date -u +%FT%TZ >> "$WORK/started_utc_$HOST.txt"   # one line per run segment (resume appends)
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit_$HOST.txt"
PYTHONPATH="$GF" "$PY" "$BD/scripts/run_isp.py" --manifest "$DATA/goals/donor_manifest.json" \
  --host-map "$BD/controls/pre_gpu_host_map.csv" --host "$HOST" \
  --token-dict "$GF/geneformer/token_dictionary_gc104M.pkl" --bf16-bench "$REPO/sclc_validation/bf16_bench" \
  --scripts "$BD/scripts" --out "$WORK"
date -u +%FT%TZ >> "$WORK/finished_utc_$HOST.txt"
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_end_$HOST.csv"
echo PHASE6_DONE
