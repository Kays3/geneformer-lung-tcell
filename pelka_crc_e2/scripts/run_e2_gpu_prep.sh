#!/usr/bin/env bash
# E2 GPU steps 1-3 on thinkstation1 (registration s.5): five fold fine-tunes (balanced_donor_luad
# finetune_fold.py, unchanged), the classifier gate (classifier_gate.py, unchanged; a FAIL stops here),
# per-donor goal centroids and ISP inputs (goal_embeddings.py, unchanged), and the no-op gate
# (pre_isp_gates.py --mode noop, unchanged; a FAIL stops here). GPU sampled every 60 s.
# Usage: run_e2_gpu_prep.sh <repo_worktree> <data_root> <geneformer_root> <venv_python>
set -euo pipefail
REPO="$1"; DATA="$2"; GF="$3"; PY="$4"
E2="$REPO/pelka_crc_e2"; BD="$REPO/balanced_donor_luad"; BF="$REPO/sclc_validation/bf16_bench"
DS="$DATA/data/tokenized/pelka_crc_e2_pool300.dataset"
WORK="$DATA/phase4"; GOALS="$DATA/goals"; mkdir -p "$WORK" "$GOALS"
sha256sum -c "$E2/registration/pins_prep.sha256"
nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$DATA/nvidia_smi_samples_prep.csv" 2>&1 &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$DATA/gpu_procs_at_prep_start.csv"
date -u +%FT%TZ >> "$WORK/started_utc.txt"
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit.txt"
for k in 0 1 2 3 4; do
  if [ -f "$WORK/fold$k/fold_record.json" ]; then echo "fold $k done, skipping"; continue; fi
  PYTHONPATH="$GF" "$PY" "$BD/scripts/finetune_fold.py" --dataset "$DS" --folds "$E2/cohort/folds.json" --fold $k \
    --base-model "$GF/Geneformer-V2-316M" --work "$WORK" --bf16-bench "$BF" > "$WORK/fold$k.log" 2>&1
  echo "fold $k finished $(date -u +%FT%TZ)"
done
date -u +%FT%TZ > "$WORK/finished_utc.txt"
if ! "$PY" "$BD/scripts/classifier_gate.py" --work "$WORK" --out "$WORK/classifier_gate.json"; then
  echo "STOP: classifier gate FAILED (registration s.6.1); no goals, no ISP"; exit 10
fi
if [ ! -f "$GOALS/donor_manifest.json" ]; then
  date -u +%FT%TZ > "$GOALS/started_utc.txt"
  PYTHONPATH="$GF" "$PY" "$BD/scripts/goal_embeddings.py" --dataset "$DS" --folds "$E2/cohort/folds.json" \
    --phase4 "$WORK" --bf16-bench "$BF" --scripts "$BD/scripts" --out "$GOALS" > "$GOALS/goals.log" 2>&1
  date -u +%FT%TZ > "$GOALS/finished_utc.txt"
fi
FIRST_DONOR=$("$PY" -c "import json; print(sorted(json.load(open('$GOALS/donor_manifest.json')))[0])")
NOOP_GENE=ENSG00000166710   # B2M, the gene the balanced-donor no-op gate used
PYTHONPATH="$GF" "$PY" "$BD/scripts/pre_isp_gates.py" --mode noop --manifest "$GOALS/donor_manifest.json" \
  --donor "$FIRST_DONOR" --gene "$NOOP_GENE" --bf16-bench "$BF" --scripts "$BD/scripts" --out "$DATA/noop_gate" \
  > "$DATA/noop_gate.log" 2>&1
"$PY" -c "import json,sys; r=json.load(open('$DATA/noop_gate/noop_gate.json')); print(r); sys.exit(0 if r['PASS'] else 11)" \
  || { echo "STOP: no-op gate FAILED"; exit 11; }
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$DATA/gpu_procs_at_prep_end.csv"
echo E2_PREP_DONE
