#!/usr/bin/env bash
# Continuation of the 316M second-stage chain (2026-09-18), stages 2+3
# only: state-embedding rebuild -> full panel arm. Split out from
# run_316m_stage_chain.sh's stage 1 because stage 1 (fine-tune) already
# completed successfully in a prior run -- training itself succeeded;
# only the post-training eval-metrics computation crashed on a real
# Geneformer library bug (see eval_tie_fix.py), recovered without
# retraining via recover_316m_test_metrics.py.
#
# Per god's ruling: the ONE manual boundary is a human/agent reading the
# recovered macro_f1 and launching this script if it clears the 0.60
# sanity threshold. Once launched, this is detached and unattended for
# both remaining stages, same discipline as run_316m_stage_chain.sh's
# original three-stage design -- no interactive stage boundary inside
# this script.
set -uo pipefail

HOME_DIR="$(cd ~ && pwd)"
REPO_ROOT="$HOME_DIR/workspace/geneformer-lung-tcell"
BF16_BENCH="$REPO_ROOT/sclc_validation/bf16_bench"
TARGETED_PANEL_DIR="$REPO_ROOT/sclc_validation/perturbation_workflow/targeted_panel"
CHAIN_ROOT="$BF16_BENCH/runs/316m_bf16_chain"
MARKERS="$CHAIN_ROOT/markers"
mkdir -p "$MARKERS"

FINETUNE_WORK_DIR="$BF16_BENCH/runs/316m_bf16_finetune"
CENTROIDS_DIR="$FINETUNE_WORK_DIR/state_embeddings"
CENTROIDS_FILE="$CENTROIDS_DIR/real_316m_centroids.pkl"

VENV_PY="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin/python3"
GENEFORMER_ROOT="$HOME_DIR/workspace/geneformer-uv-starter/Geneformer"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

MODEL_DIR="$(cat "$FINETUNE_WORK_DIR/runs/MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt")"
log "Resuming with recovered checkpoint: $MODEL_DIR"

# ---------- Stage 2: state-embedding rebuild ----------
log "STAGE 2/3: state-embedding rebuild starting"
STAGE2_START=$(date +%s)
PATH="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin:$PATH" \
  GENEFORMER_ROOT="$GENEFORMER_ROOT" \
  MODEL_DIRECTORY="$MODEL_DIR" \
  OUTPUT_DIR="$CENTROIDS_DIR" \
  OUTPUT_PREFIX="real_316m_centroids" \
  bash "$BF16_BENCH/power_sample.sh" 316m_bf16_state_emb "$CENTROIDS_DIR/power" \
  -- "$VENV_PY" "$BF16_BENCH/calc_state_embeddings_316m.py" \
  > "$CHAIN_ROOT/stage2_state_emb.log" 2>&1
STAGE2_EXIT=$?
STAGE2_WALL=$(( $(date +%s) - STAGE2_START ))

if [ "$STAGE2_EXIT" -ne 0 ]; then
  log "STAGE 2 FAILED (exit $STAGE2_EXIT) -- aborting chain"
  python3 -c "
import json
json.dump({'stage': 'state_emb', 'status': 'failed', 'exit_code': $STAGE2_EXIT, 'wall_s': $STAGE2_WALL}, open('$MARKERS/stage2_state_emb.json', 'w'), indent=2)
"
  exit 3
fi
log "STAGE 2 done: wall=${STAGE2_WALL}s centroids=$CENTROIDS_FILE"
python3 -c "
import json
json.dump({'stage': 'state_emb', 'status': 'done', 'wall_s': $STAGE2_WALL, 'model_directory': '$MODEL_DIR', 'centroids_file': '$CENTROIDS_FILE'}, open('$MARKERS/stage2_state_emb.json', 'w'), indent=2)
"

# ---------- Stage 3: full panel arm ----------
log "STAGE 3/3: full 316M-bf16 panel arm starting (150 units, overexpress-only, max-ncells 300)"
STAGE3_START=$(date +%s)
cd "$TARGETED_PANEL_DIR"
PATH="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin:$PATH" \
  GENEFORMER_ROOT="$GENEFORMER_ROOT" \
  HTAN_FINETUNE_ROOT="$FINETUNE_WORK_DIR" \
  STATE_EMB_FILE_OVERRIDE="$CENTROIDS_FILE" \
  bash "$BF16_BENCH/power_sample.sh" 316m_bf16_panel "$BF16_BENCH/runs/316m_bf16/panel_power" \
  -- "$VENV_PY" run_targeted_panel.py --dtype bf16 --run-tag 316m_bf16 --max-ncells 300 --perturb-types overexpress \
  > "$CHAIN_ROOT/stage3_panel.log" 2>&1
STAGE3_EXIT=$?
STAGE3_WALL=$(( $(date +%s) - STAGE3_START ))

if [ "$STAGE3_EXIT" -ne 0 ]; then
  log "STAGE 3 FAILED (exit $STAGE3_EXIT)"
  python3 -c "
import json
json.dump({'stage': 'panel', 'status': 'failed', 'exit_code': $STAGE3_EXIT, 'wall_s': $STAGE3_WALL}, open('$MARKERS/stage3_panel.json', 'w'), indent=2)
"
  exit 4
fi
log "STAGE 3 done: wall=${STAGE3_WALL}s"
python3 -c "
import json
json.dump({'stage': 'panel', 'status': 'done', 'wall_s': $STAGE3_WALL}, open('$MARKERS/stage3_panel.json', 'w'), indent=2)
"

python3 -c "
import json
json.dump({'chain': '316m_bf16_stages_2_3', 'status': 'done', 'stages': ['state_emb', 'panel']}, open('$MARKERS/CHAIN_DONE.json', 'w'), indent=2)
"
log "CHAIN COMPLETE. Canary is the manual next step -- report to god, then select top-10 |shift| genes from the panel results."
