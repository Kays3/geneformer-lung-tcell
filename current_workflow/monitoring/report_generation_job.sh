#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DEFAULT_STATS_DIR="/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/stats"

if [[ -z "${PERTURBATION_STATS_DIR:-}" && -d "$DEFAULT_STATS_DIR" ]]; then
  export PERTURBATION_STATS_DIR="$DEFAULT_STATS_DIR"
fi


exec "$PYTHON_BIN" "$DIR/generate_progress_report.py"
