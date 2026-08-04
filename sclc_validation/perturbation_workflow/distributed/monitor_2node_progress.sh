#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd)
WORKSPACE_ROOT=$(cd -- "$REPO_ROOT/.." && pwd)
ANALYSIS_ROOT="$WORKSPACE_ROOT/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"
REMOTE_USER=${REMOTE_USER:-thinkstation1}
REMOTE_IP=${REMOTE_IP:-192.168.100.1}
REMOTE_SSH="$REMOTE_USER@$REMOTE_IP"
REMOTE_ROOT=${REMOTE_ROOT:-/home/thinkstation1/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation}
LOG=${PROGRESS_LOG:-$REPO_ROOT/sclc_validation/primary_test_perturbation/progress.log}
INTERVAL=${INTERVAL:-300}

mkdir -p "$(dirname -- "$LOG")"
while :; do
    timestamp=$(date --iso-8601=seconds)
    local_counts=$(for arm in delete overexpress; do for source in normal sclc luad; do n=$(find "$ANALYSIS_ROOT/raw/$arm/$source" -maxdepth 1 -name 'heldout_*complete.json' 2>/dev/null | wc -l); printf '%s/%s=%s ' "$arm" "$source" "$n"; done; done)
    remote_counts=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE_SSH" "for arm in delete overexpress; do for source in normal sclc luad; do n=\$(find '$REMOTE_ROOT/raw/'\$arm/'/'\$source -maxdepth 1 -name 'heldout_*complete.json' 2>/dev/null | wc -l); printf '%s/%s=%s ' \"\$arm\" \"\$source\" \"\$n\"; done; done" 2>/dev/null || printf 'unavailable')
    local_gpu=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,power.draw --format=csv,noheader 2>/dev/null | tr '\n' ';' || printf 'unavailable')
    remote_gpu=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE_SSH" 'nvidia-smi --query-gpu=utilization.gpu,memory.used,power.draw --format=csv,noheader 2>/dev/null | tr "\n" ";"' 2>/dev/null || printf 'unavailable')
    local_pid=$(if [[ -f "$ANALYSIS_ROOT/logs/2node/local_delete.pid" ]]; then ps -p "$(<"$ANALYSIS_ROOT/logs/2node/local_delete.pid")" -o pid=,stat= 2>/dev/null | tr -s ' '; else printf 'none'; fi)
    remote_pid=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE_SSH" "if [[ -f '$REMOTE_ROOT/logs/2node/remote_overexpress.pid' ]]; then ps -p \$(cat '$REMOTE_ROOT/logs/2node/remote_overexpress.pid') -o pid=,stat= 2>/dev/null | tr -s ' '; else printf 'none'; fi" 2>/dev/null || printf 'unavailable')
    printf '%s | monitor: local_pid=%s remote_pid=%s | local_shards=%s|remote_shards=%s|local_gpu=%s|remote_gpu=%s\n' "$timestamp" "$local_pid" "$remote_pid" "$local_counts" "$remote_counts" "$local_gpu" "$remote_gpu" >> "$LOG"
    sleep "$INTERVAL"
done
