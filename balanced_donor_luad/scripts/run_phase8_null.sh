#!/usr/bin/env bash
# Amendment 4 driver: runs balanced_donor_luad/scripts/run_isp.py, UNCHANGED, against the frozen
# null_genes.json draw order, under a hard GPU-arm budget ceiling.
#
# The ceiling is enforced by reading the SAME "seconds" field run_isp.py writes into every
# *.complete.json marker (ISP-STD-1 E.3: count GPU arms, not wall clock). It never inspects
# nvidia-smi or wall-clock elapsed time to decide whether to stop. nvidia-smi is still sampled
# for the audit trail, exactly as Phase 6 did, but is not part of the stop decision.
#
# Usage: run_phase8_null.sh <host> <repo_worktree> <data_root> <geneformer_root> <venv_python> <ceiling_seconds>
set -euo pipefail
HOST="$1"; REPO="$2"; DATA="$3"; GF="$4"; PY="$5"; CEILING="${6:-36000}"
BD="$REPO/balanced_donor_luad"
NULLDIR="$DATA/phase8_null"; WORK="$NULLDIR/out"
mkdir -p "$WORK"

GENES=$(python3 -c "import json; print(','.join(json.load(open('$NULLDIR/null_genes.json'))['draw_order']))")

nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$WORK/nvidia_smi_samples_$HOST.csv" 2>&1 &
SAMPLER=$!
date -u +%FT%TZ >> "$WORK/started_utc_$HOST.txt"
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit_$HOST.txt"

# Launch run_isp.py detached in its own process group (setsid), so the kill below targets exactly
# this subtree and never the ssh/bash shell that launched it (the pgrep -f self-match trap).
setsid env PYTHONPATH="$GF" "$PY" "$BD/scripts/run_isp.py" \
  --manifest "$DATA/goals/donor_manifest.json" --host-map "$BD/controls/pre_gpu_host_map.csv" \
  --host "$HOST" --genes "$GENES" \
  --token-dict "$GF/geneformer/token_dictionary_gc104M.pkl" --bf16-bench "$REPO/sclc_validation/bf16_bench" \
  --scripts "$BD/scripts" --out "$WORK" > "$WORK/run_log_stdout_$HOST.log" 2>&1 &
ISP_PID=$!
ISP_PGID=$(ps -o pgid= -p "$ISP_PID" | tr -d ' ')
echo "isp_pid=$ISP_PID isp_pgid=$ISP_PGID" | tee "$WORK/pgid_$HOST.txt"

STOP_REASON="completed"
while kill -0 "$ISP_PID" 2>/dev/null; do
  sleep 30
  SPENT=$(python3 -c "
import glob, json
s = 0.0
for f in glob.glob('$WORK/*/*/*.complete.json'):
    try:
        s += json.load(open(f)).get('seconds', 0.0)
    except Exception:
        pass
print(s)
")
  # Verify the PID we are about to kill really is the run_isp.py we launched, not a coincidental reuse.
  REAL=$(ps -o pid,pgid,cmd -p "$ISP_PID" 2>/dev/null | tail -n +2 | grep -c "run_isp.py" || true)
  echo "$(date -u +%FT%TZ) spent_s=$SPENT ceiling_s=$CEILING pid_verified=$REAL" >> "$WORK/ceiling_monitor_$HOST.log"
  if [ "$REAL" = "0" ]; then
    break  # process already gone or PID reused by something else; the while-loop's kill -0 will resolve this
  fi
  if python3 -c "exit(0 if $SPENT >= $CEILING else 1)"; then
    echo "$(date -u +%FT%TZ) CEILING REACHED spent_s=$SPENT >= $CEILING -- stopping pgid $ISP_PGID" >> "$WORK/ceiling_monitor_$HOST.log"
    kill -TERM -"$ISP_PGID" 2>/dev/null || true
    sleep 5
    kill -KILL -"$ISP_PGID" 2>/dev/null || true
    STOP_REASON="ceiling"
    break
  fi
done
wait "$ISP_PID" 2>/dev/null || true

kill "$SAMPLER" 2>/dev/null || true
date -u +%FT%TZ >> "$WORK/finished_utc_$HOST.txt"
echo "$STOP_REASON" > "$WORK/stop_reason_$HOST.txt"
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_end_$HOST.csv"
echo "PHASE8_NULL_DONE reason=$STOP_REASON"
