#!/usr/bin/env bash
set -euo pipefail

# Run the existing resumable Geneformer SCLC all-gene screen on two independent
# GPUs. Geneformer's InSilicoPerturber is not a torchrun/DDP workload; splitting
# perturbation types avoids CUDA/NCCL contention while using the direct link for
# the one-time asset transfer and result collection.

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd)
WORKSPACE_ROOT=$(cd -- "$REPO_ROOT/.." && pwd)
ANALYSIS_ROOT="$WORKSPACE_ROOT/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"
FINETUNE_ROOT="$WORKSPACE_ROOT/KD/sclc_luad_normal_htan_finetune"
RUNNER="$ANALYSIS_ROOT/scripts/run_heldout_allgene.py"

REMOTE_USER=${REMOTE_USER:-thinkstation1}
REMOTE_IP=${REMOTE_IP:-192.168.100.1}
REMOTE_HOME=${REMOTE_HOME:-/home/thinkstation1}
REMOTE_WORKSPACE=${REMOTE_WORKSPACE:-$REMOTE_HOME/workspace}
REMOTE_PYTHON=${REMOTE_PYTHON:-$REMOTE_WORKSPACE/geneformer-uv-starter/.venv/bin/python}
LOCAL_PYTHON=${LOCAL_PYTHON:-$WORKSPACE_ROOT/geneformer-uv-starter/.venv/bin/python}
IFACE=${IFACE:-enp1s0f0np0}
REMOTE_SSH="$REMOTE_USER@$REMOTE_IP"
REMOTE_ANALYSIS_ROOT="$REMOTE_WORKSPACE/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"
REMOTE_FINETUNE_ROOT="$REMOTE_WORKSPACE/KD/sclc_luad_normal_htan_finetune"
LOG_ROOT="$ANALYSIS_ROOT/logs/2node"
PROGRESS_LOG=${PROGRESS_LOG:-$REPO_ROOT/sclc_validation/primary_test_perturbation/progress.log}

log_event() {
    mkdir -p "$(dirname -- "$PROGRESS_LOG")"
    printf '%s | %s\n' "$(date --iso-8601=seconds)" "$*" | tee -a "$PROGRESS_LOG"
}

usage() {
    printf '%s\n' \
        "Usage: $0 {prepare|start|sync|stats|status}" \
        "  prepare  verify link, transfer assets, and run remote smoke test" \
        "  start    launch deletion locally and overexpression remotely" \
        "  sync     copy remote overexpression raw outputs back locally" \
        "  stats    calculate local stats after both arms are synchronized" \
        "  status   show local and remote process/output status"
}

ssh_remote() {
    ssh -o BatchMode=yes -o ConnectTimeout=10 -o ServerAliveInterval=30 "$REMOTE_SSH" "$@"
}

sync_assets() {
    ssh_remote "mkdir -p '$REMOTE_WORKSPACE/KD'"
    rsync -aH --partial --info=progress2 \
        "$FINETUNE_ROOT/" "$REMOTE_SSH:$REMOTE_FINETUNE_ROOT/"
    rsync -aH --partial --info=progress2 \
        --exclude='raw/' --exclude='stats/' --exclude='logs/' \
        "$ANALYSIS_ROOT/" "$REMOTE_SSH:$REMOTE_ANALYSIS_ROOT/"
    rsync -aH --partial --info=progress2 \
        "$RUNNER" "$REMOTE_SSH:$REMOTE_ANALYSIS_ROOT/scripts/run_heldout_allgene.py"
    printf '%s\n' "$REMOTE_FINETUNE_ROOT/runs/260728_geneformer_cellClassifier_sclc_luad_normal_htan/ksplit1" \
        | ssh_remote "cat > '$REMOTE_FINETUNE_ROOT/runs/MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt'"
}

prepare() {
    log_event "prepare: validating direct link $IFACE -> $REMOTE_IP"
    ip route get "$REMOTE_IP" | grep -q "dev $IFACE"
    ping -c 3 -W 1 "$REMOTE_IP"
    ssh_remote "test -x '$REMOTE_PYTHON'"
    log_event "prepare: synchronizing SCLC assets to $REMOTE_SSH"
    sync_assets
    log_event "prepare: running remote overexpression smoke test"
    ssh_remote "NCCL_SOCKET_IFNAME=$IFACE '$REMOTE_PYTHON' '$REMOTE_ANALYSIS_ROOT/scripts/run_heldout_allgene.py' smoke-test --perturb-types overexpress --forward-batch-size 16 --nproc 4"
    log_event "prepare: running local deletion smoke test"
    "$LOCAL_PYTHON" "$RUNNER" smoke-test --perturb-types delete --forward-batch-size 16 --nproc 4
    log_event "prepare: both smoke tests passed"
}

start() {
    mkdir -p "$LOG_ROOT"
    log_event "start: synchronizing assets and launching primary test-set perturbations"
    sync_assets
    nohup env CUDA_VISIBLE_DEVICES=0 NCCL_SOCKET_IFNAME="$IFACE" \
        "$LOCAL_PYTHON" "$RUNNER" perturb --perturb-types delete \
        --sources normal sclc luad --forward-batch-size 16 --nproc 4 \
        >"$LOG_ROOT/local_delete.log" 2>&1 &
    printf '%s\n' "$!" >"$LOG_ROOT/local_delete.pid"
    ssh_remote "mkdir -p '$REMOTE_ANALYSIS_ROOT/logs/2node'; nohup env CUDA_VISIBLE_DEVICES=0 NCCL_SOCKET_IFNAME=$IFACE '$REMOTE_PYTHON' '$REMOTE_ANALYSIS_ROOT/scripts/run_heldout_allgene.py' perturb --perturb-types overexpress --sources normal sclc luad --forward-batch-size 16 --nproc 4 >'$REMOTE_ANALYSIS_ROOT/logs/2node/remote_overexpress.log' 2>&1 & echo \$! >'$REMOTE_ANALYSIS_ROOT/logs/2node/remote_overexpress.pid'"
    log_event "start: local deletion pid $(<"$LOG_ROOT/local_delete.pid"); remote overexpression launched"
    printf '%s\n' "Started local deletion and remote overexpression. Logs: $LOG_ROOT and $REMOTE_ANALYSIS_ROOT/logs/2node"
}

sync_results() {
    log_event "sync: collecting remote overexpression raw outputs"
    mkdir -p "$ANALYSIS_ROOT/raw/overexpress"
    rsync -aH --partial --info=progress2 \
        "$REMOTE_SSH:$REMOTE_ANALYSIS_ROOT/raw/overexpress/" \
        "$ANALYSIS_ROOT/raw/overexpress/"
    log_event "sync: remote overexpression outputs synchronized"
}

stats() {
    log_event "stats: building Geneformer primary arm statistics"
    "$LOCAL_PYTHON" "$RUNNER" stats --perturb-types delete overexpress --sources normal sclc luad
    log_event "stats: Geneformer primary arm statistics complete"
}

status() {
    printf '%s\n' 'Local:'
    if [[ -f "$LOG_ROOT/local_delete.pid" ]]; then
        ps -p "$(<"$LOG_ROOT/local_delete.pid")" -o pid=,stat=,etime=,cmd= || true
    fi
    for type in delete overexpress; do
        for source in normal sclc luad; do
            count=$(find "$ANALYSIS_ROOT/raw/$type/$source" -maxdepth 1 -name 'heldout_*complete.json' 2>/dev/null | wc -l)
            printf '  %s/%s complete markers: %s\n' "$type" "$source" "$count"
        done
    done
    printf '%s\n' 'Remote:'
    ssh_remote "if [[ -f '$REMOTE_ANALYSIS_ROOT/logs/2node/remote_overexpress.pid' ]]; then ps -p \$(cat '$REMOTE_ANALYSIS_ROOT/logs/2node/remote_overexpress.pid') -o pid=,stat=,etime=,cmd= || true; fi; for source in normal sclc luad; do printf '  overexpress/%s complete markers: '; find '$REMOTE_ANALYSIS_ROOT/raw/overexpress/'\$source -maxdepth 1 -name 'heldout_*complete.json' 2>/dev/null | wc -l; done"
}

case "${1:-}" in
    prepare) prepare ;;
    start) start ;;
    sync) sync_results ;;
    stats) stats ;;
    status) status ;;
    *) usage >&2; exit 2 ;;
esac
