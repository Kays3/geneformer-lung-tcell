#!/usr/bin/env bash
# One chained, detached script for the 316M second stage (2026-09-18),
# per god's ruling #3: fine-tune -> sanity-abort check -> state-embedding
# rebuild -> full panel arm, no interactive stage boundaries (the recorded
# lesson from three prior compaction-driven stalls on the 104M stage's
# individually-launched arms). The precision canary is the one intentional
# manual boundary after this (it needs the panel's own top-|shift| genes
# for gene selection) and is run separately once this chain reports done.
#
# Each stage is power-sampled independently and writes its own JSON
# completion marker under $CHAIN_ROOT/markers/ so progress survives even
# if the launching session goes idle between polls.
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
MACRO_F1_MIN="0.60"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

# ---------- Stage 1: fine-tune ----------
log "STAGE 1/3: fine-tune starting"
STAGE1_START=$(date +%s)
PATH="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin:$PATH" \
  bash "$BF16_BENCH/power_sample.sh" 316m_bf16_finetune "$FINETUNE_WORK_DIR/power" \
  -- "$VENV_PY" "$BF16_BENCH/run_finetune_316m.py" \
  > "$CHAIN_ROOT/stage1_finetune.log" 2>&1
STAGE1_EXIT=$?
STAGE1_WALL=$(( $(date +%s) - STAGE1_START ))

if [ "$STAGE1_EXIT" -ne 0 ]; then
  log "STAGE 1 FAILED (exit $STAGE1_EXIT) -- aborting chain"
  python3 -c "
import json
json.dump({'stage': 'finetune', 'status': 'failed', 'exit_code': $STAGE1_EXIT, 'wall_s': $STAGE1_WALL}, open('$MARKERS/stage1_finetune.json', 'w'), indent=2)
"
  exit 1
fi

TEST_METRICS="$FINETUNE_WORK_DIR/tables/test_metrics.json"
MACRO_F1=$(python3 -c "
import json
d = json.load(open('$TEST_METRICS'))
print(d['test_metrics']['macro_f1'])
")
ACC=$(python3 -c "
import json
d = json.load(open('$TEST_METRICS'))
print(d['test_metrics']['acc'])
")
log "STAGE 1 done: wall=${STAGE1_WALL}s macro_f1=${MACRO_F1} acc=${ACC}"

SANITY_PASS=$(python3 -c "print(1 if float('$MACRO_F1') >= $MACRO_F1_MIN else 0)")
python3 -c "
import json
json.dump({
    'stage': 'finetune', 'status': 'done', 'wall_s': $STAGE1_WALL,
    'macro_f1': $MACRO_F1, 'acc': $ACC,
    'sanity_threshold': $MACRO_F1_MIN, 'sanity_pass': bool($SANITY_PASS),
}, open('$MARKERS/stage1_finetune.json', 'w'), indent=2)
"

if [ "$SANITY_PASS" -eq 0 ]; then
  log "SANITY ABORT: macro_f1=${MACRO_F1} < ${MACRO_F1_MIN} -- broken-training guard tripped, stopping chain (state-emb + panel arm NOT run)"
  exit 2
fi
log "Sanity check passed (macro_f1=${MACRO_F1} >= ${MACRO_F1_MIN}); continuing chain"

# ---------- Stage 2: state-embedding rebuild ----------
log "STAGE 2/3: state-embedding rebuild starting"
MODEL_DIR="$(cat "$FINETUNE_WORK_DIR/runs/MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt")"
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
json.dump({'chain': '316m_bf16', 'status': 'done', 'stages': ['finetune', 'state_emb', 'panel']}, open('$MARKERS/CHAIN_DONE.json', 'w'), indent=2)
"
log "CHAIN COMPLETE. Canary is the manual next step -- report to god, then select top-10 |shift| genes from the panel results."
