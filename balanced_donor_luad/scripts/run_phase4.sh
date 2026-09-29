#!/usr/bin/env bash
# Phase 4 driver (thinkstation1): 5 fold fine-tunes in sequence, GPU sampled every 60 s.
# Usage: run_phase4.sh <repo_worktree> <data_root> <geneformer_root> <venv_python>
set -euo pipefail
REPO="$1"; DATA="$2"; GF="$3"; PY="$4"
BD="$REPO/balanced_donor_luad"; WORK="$DATA/phase4"; mkdir -p "$WORK"
nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 > "$WORK/nvidia_smi_samples.csv" 2>&1 &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_start.csv"
date -u +%FT%TZ > "$WORK/started_utc.txt"
git -C "$REPO" rev-parse HEAD > "$WORK/code_commit.txt"
for k in 0 1 2 3 4; do
  if [ -f "$WORK/fold$k/fold_record.json" ]; then echo "fold $k done, skipping"; continue; fi
  PYTHONPATH="$GF" "$PY" "$BD/scripts/finetune_fold.py" --dataset "$DATA/data/tokenized/balanced_donor_luad_pool300.dataset" \
    --folds "$BD/cohort/folds.json" --fold $k --base-model "$GF/Geneformer-V2-316M" --work "$WORK" \
    --bf16-bench "$REPO/sclc_validation/bf16_bench" > "$WORK/fold$k.log" 2>&1
  echo "fold $k finished $(date -u +%FT%TZ)"
done
date -u +%FT%TZ > "$WORK/finished_utc.txt"
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_end.csv"
echo PHASE4_DONE
