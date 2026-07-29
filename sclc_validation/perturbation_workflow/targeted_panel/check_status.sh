#!/usr/bin/env bash
set -euo pipefail

run_dir="${TARGETED_PANEL_RUN_DIR:-$HOME/workspace/KD/sclc_luad_normal_htan_targeted_panel_perturbation}"
log_file="$run_dir/logs/targeted_panel.log"
config_file="$run_dir/tables/run_config.json"
runner="$run_dir/scripts/run_targeted_panel.py"

if [[ ! -d "$run_dir" ]]; then
    printf 'Run directory not found: %s\n' "$run_dir" >&2
    exit 2
fi

target_genes=50
if command -v jq >/dev/null 2>&1 && [[ -f "$config_file" ]]; then
    target_genes="$(jq -r '.n_target_genes // 50' "$config_file")"
fi

printf 'LAST STAGE: %s\n' "$(grep -E '===' "$log_file" 2>/dev/null | tail -n 1 || true)"
printf 'LAST COMPLETION: %s\n' \
    "$(grep -E 'complete in|Targeted panel perturbation complete' "$log_file" 2>/dev/null | tail -n 1 || true)"

completed_total=0
for perturb_type in delete overexpress; do
    for source in normal sclc luad; do
        output_dir="$run_dir/raw/$perturb_type/$source"
        completed="$(find "$output_dir" -maxdepth 1 -type f -name 'targeted_*.complete.json' 2>/dev/null | wc -l)"
        completed_total=$((completed_total + completed))
        printf '%s/%s: %s/%s genes\n' "$perturb_type" "$source" "$completed" "$target_genes"
    done
done

stats_count="$(find "$run_dir/stats" -type f -name '*.csv' 2>/dev/null | wc -l)"
expected_total=$((target_genes * 6))
printf 'STATS: %s/12 tables\n' "$stats_count"

recent_failures="$(
    tail -n 250 "$log_file" 2>/dev/null \
        | grep -E 'Traceback|CUDA out of memory|Cannot re-initialize CUDA|Killed|AssertionError' \
        || true
)"
if [[ -n "$recent_failures" ]]; then
    printf 'RECENT FAILURES:\n%s\n' "$recent_failures"
fi

if pgrep -f -- "$runner" >/dev/null 2>&1; then
    printf 'RUN: running\n'
elif [[ "$completed_total" -eq "$expected_total" ]] \
    && [[ "$stats_count" -eq 12 ]] \
    && grep -Fq 'Targeted panel perturbation complete.' "$log_file"; then
    printf 'RUN: complete\n'
elif [[ -n "$recent_failures" ]]; then
    printf 'RUN: failed\n'
    exit 1
else
    printf 'RUN: incomplete (%s/%s gene-runs)\n' "$completed_total" "$expected_total"
    exit 1
fi
