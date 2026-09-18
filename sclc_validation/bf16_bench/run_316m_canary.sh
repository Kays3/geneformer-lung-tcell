#!/usr/bin/env bash
# 316M precision canary (2026-09-18), per god's ruling: 10 genes (top-10 by
# MAX |Shift_to_goal_end| across source-goal pairs in the 316M-bf16 panel),
# x 3 sources x 2 dtypes = 60 units. One detached run covering BOTH dtype
# passes -- no interactive stage boundary between them. The per-architecture
# 316M centroids (real_316m_centroids.pkl) are reused for both dtypes, same
# choice already validated by the 104M stage's frozen-gate pass. Explicitly
# a post-hoc, risk-bounding canary -- not a G5 gate, no top-20 verdict at
# this n.
set -uo pipefail

HOME_DIR="$(cd ~ && pwd)"
REPO_ROOT="$HOME_DIR/workspace/geneformer-lung-tcell"
BF16_BENCH="$REPO_ROOT/sclc_validation/bf16_bench"
TARGETED_PANEL_DIR="$REPO_ROOT/sclc_validation/perturbation_workflow/targeted_panel"
CHAIN_ROOT="$BF16_BENCH/runs/316m_bf16_chain"
MARKERS="$CHAIN_ROOT/markers"
mkdir -p "$MARKERS"

FINETUNE_WORK_DIR="$BF16_BENCH/runs/316m_bf16_finetune"
CENTROIDS_FILE="$FINETUNE_WORK_DIR/state_embeddings/real_316m_centroids.pkl"
CANARY_GENES="$BF16_BENCH/runs/316m_canary_genes.json"

VENV_PY="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin/python3"
GENEFORMER_ROOT="$HOME_DIR/workspace/geneformer-uv-starter/Geneformer"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

run_canary_dtype() {
  local dtype="$1"
  local run_tag="316m_canary_${dtype}"
  log "CANARY ${dtype}: starting (10 genes x 3 sources, overexpress, max-ncells 300)"
  local start
  start=$(date +%s)
  cd "$TARGETED_PANEL_DIR"
  PATH="$HOME_DIR/workspace/geneformer-uv-starter/sclc_analysis/.venv/bin:$PATH" \
    GENEFORMER_ROOT="$GENEFORMER_ROOT" \
    HTAN_FINETUNE_ROOT="$FINETUNE_WORK_DIR" \
    STATE_EMB_FILE_OVERRIDE="$CENTROIDS_FILE" \
    TARGET_GENES_FILE_OVERRIDE="$CANARY_GENES" \
    bash "$BF16_BENCH/power_sample.sh" "$run_tag" "$BF16_BENCH/runs/$run_tag/power" \
    -- "$VENV_PY" run_targeted_panel.py --dtype "$dtype" --run-tag "$run_tag" --max-ncells 300 --perturb-types overexpress \
    > "$CHAIN_ROOT/canary_${dtype}.log" 2>&1
  local exit_code=$?
  local wall=$(( $(date +%s) - start ))
  if [ "$exit_code" -ne 0 ]; then
    log "CANARY ${dtype} FAILED (exit $exit_code)"
    python3 -c "
import json
json.dump({'stage': 'canary_${dtype}', 'status': 'failed', 'exit_code': $exit_code, 'wall_s': $wall}, open('$MARKERS/canary_${dtype}.json', 'w'), indent=2)
"
    return 1
  fi
  log "CANARY ${dtype} done: wall=${wall}s"
  python3 -c "
import json
json.dump({'stage': 'canary_${dtype}', 'status': 'done', 'wall_s': $wall}, open('$MARKERS/canary_${dtype}.json', 'w'), indent=2)
"
  return 0
}

run_canary_dtype fp32 || exit 1
run_canary_dtype bf16 || exit 2

python3 -c "
import json
json.dump({'chain': '316m_canary', 'status': 'done', 'stages': ['canary_fp32', 'canary_bf16']}, open('$MARKERS/CANARY_DONE.json', 'w'), indent=2)
"
log "CANARY COMPLETE (both dtypes). Compute rho/sign-agreement/max-delta from the two stats roots next."
