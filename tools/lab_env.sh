#!/usr/bin/env bash
# Resolve every machine-specific path the analyses need, in one place.
#
# Source it, do not execute it:
#
#   source tools/lab_env.sh
#   python sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py
#
# Resolution order, first hit wins per variable:
#   1. anything already exported in the environment  (explicit override)
#   2. ~/.config/geneformer-lung-tcell/paths.env     (per-machine, untracked)
#   3. the /srv/lab defaults below                   (post-migration layout)
#
# The scripts in this repository all read these variables and fall back to their
# own defaults, so sourcing this is optional but makes a machine's layout
# explicit in one auditable place instead of scattered across call sites. Paths
# are never hardcoded in tracked analysis code -- see WORKFLOW.md.

_LAB_ROOT="${LAB_ROOT:-/srv/lab}"
_USER_ENV="${LAB_ENV_FILE:-$HOME/.config/geneformer-lung-tcell/paths.env}"

# Per-machine overrides, if present.
# shellcheck disable=SC1090
[[ -r "$_USER_ENV" ]] && source "$_USER_ENV"

# Defaults only fill variables the caller has not already set.
: "${SCLC_PERTURBATION_ROOT:=$_LAB_ROOT/KD/sclc_luad_normal_htan_heldout_allgene_perturbation}"
: "${TARGETED_PANEL_RUN_DIR:=$_LAB_ROOT/KD/sclc_luad_normal_htan_targeted_panel_perturbation}"
: "${HTAN_H5AD:=$_LAB_ROOT/KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad}"
: "${GSE263196_RAW_DIR:=$_LAB_ROOT/spatial_raw/GSE263196_RAW}"
: "${GENEFORMER_TOKEN_DICT:=$_LAB_ROOT/geneformer-uv-starter/geneformer-workspace/Geneformer/geneformer/token_dictionary_gc104M.pkl}"
: "${PYTHON_BIN:=$_LAB_ROOT/geneformer-uv-starter/.venv/bin/python}"

export SCLC_PERTURBATION_ROOT TARGETED_PANEL_RUN_DIR HTAN_H5AD \
       GSE263196_RAW_DIR GENEFORMER_TOKEN_DICT PYTHON_BIN

lab_env_check() {
    # Report which resolved paths actually exist. Missing entries are printed
    # rather than exiting, because a machine legitimately holds only the assets
    # for the arm it ran.
    local name value missing=0
    printf '\n\033[1mResolved lab paths\033[0m\n'
    for name in SCLC_PERTURBATION_ROOT TARGETED_PANEL_RUN_DIR HTAN_H5AD \
                GSE263196_RAW_DIR GENEFORMER_TOKEN_DICT PYTHON_BIN; do
        value="${!name}"
        if [[ -e "$value" ]]; then
            printf '  \033[32mok     \033[0m %-24s %s\n' "$name" "$value"
        else
            printf '  \033[31mMISSING\033[0m %-24s %s\n' "$name" "$value"
            missing=$((missing + 1))
        fi
    done
    if [[ $missing -gt 0 ]]; then
        printf '\n  %d path(s) missing.\n' "$missing"
        printf '  If the migration has not run yet: sudo bash tools/migrate_to_srv_lab.sh --inventory\n'
        printf '  If this machine legitimately lacks them, set them in %s\n' "$_USER_ENV"
    fi
    return 0
}

# Allow `bash tools/lab_env.sh` as a quick check without sourcing.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    lab_env_check
fi
