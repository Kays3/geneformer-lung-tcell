#!/usr/bin/env bash
# S100 LUAD run, 2026-09-28 (Amendment 8, 082da381). Launch only after Stanley's sign-off 2.
set -euo pipefail
R=/home/kaisar/workspace/geneformer-lung-tcell
B=$R/sclc_validation/bf16_bench/runs
H=/home/kaisar/workspace
P=$B/s100_luad_preflight_pergene_20260928
export GENEFORMER_ROOT=$H/geneformer-uv-starter/Geneformer
export SCLC_PERTURBATION_ROOT=$H/KD/sclc_luad_normal_htan_heldout_allgene_perturbation
export HTAN_FINETUNE_ROOT=$B/316m_bf16_finetune
export STATE_EMB_FILE_OVERRIDE=$B/316m_bf16_finetune/state_embeddings/real_316m_centroids.pkl
PY=$H/geneformer-uv-starter/sclc_analysis/.venv/bin/python3
cd $R/sclc_validation/perturbation_workflow/targeted_panel
[ "$(sha256sum run_targeted_panel.py | cut -c1-64)" = 974535b7b9ad6b49e33e4dfaf2a84e7070f0cc63e64bd8b271233cf77b60f754 ] || { echo "runner hash mismatch"; exit 2; }
[ "$(sha256sum $P/s100_luad_run_panel_pergene_20260928.json | cut -c1-64)" = 8a0b36683542e51fa81a655505a348f9286669df9adac22859d5f6406f34afba ] || { echo "run panel hash mismatch"; exit 2; }
date -u +%FT%TZ > $B/s100_luad_20260928.started_utc
TARGET_GENES_FILE_OVERRIDE=$R/sclc_validation/perturbation_workflow/s100_isp/s100_gene_panel_20260922.json \
  $PY run_targeted_panel.py --dtype bf16 --sources luad --run-tag s100_luad_noop_20260928 --perturb-types noop
TARGET_GENES_FILE_OVERRIDE=$P/s100_luad_run_panel_pergene_20260928.json \
  $PY run_targeted_panel.py --dtype bf16 --sources luad --run-tag s100_luad_20260928 --perturb-types delete overexpress
date -u +%FT%TZ > $B/s100_luad_20260928.finished_utc
