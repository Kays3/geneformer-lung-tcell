#!/usr/bin/env bash
# E2 GPU step 4 (registration s.5.4): run_isp.py (balanced_donor_luad, UNCHANGED) over the frozen run order
# (controls/isp_run_order.json: 100 null genes, then eligible Panel B genes, then controls), both operations,
# all 19 donors, thinkstation1 only. Hard ceiling on GPU-arm seconds summed from the markers' own "seconds"
# field (same mechanism as run_phase8_null.sh). After the run: no-op spot checks (noop_spotcheck.py,
# unchanged) on the frozen spot-check genes, first donor.
# Usage: run_e2_isp.sh <repo_worktree> <data_root> <geneformer_root> <venv_python> <ceiling_seconds>
set -euo pipefail
REPO="$1"; DATA="$2"; GF="$3"; PY="$4"; CEILING="$5"
HOST=thinkstation1
E2="$REPO/pelka_crc_e2"; BD="$REPO/balanced_donor_luad"; BF="$REPO/sclc_validation/bf16_bench"
ORDER="$E2/controls/isp_run_order.json"; MANIFEST="$DATA/goals/donor_manifest.json"
WORK="$DATA/isp/out"; mkdir -p "$WORK"
[ "$(hostname)" = "thinkstationpgx-a6e2" ] || { echo "not thinkstation1"; exit 2; }
sha256sum -c "$E2/registration/pins_isp.sha256"
sha256sum "$MANIFEST" "$DATA"/phase4/fold*/fold_record.json > "$WORK/post_prep_inputs.sha256"
GENES=$("$PY" -c "import json; print(','.join(json.load(open('$ORDER'))['genes']))")
nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$WORK/nvidia_smi_samples.csv" 2>&1 &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT
nvidia-smi --query-compute-apps=pid,process_name --format=csv > "$WORK/gpu_procs_at_start.csv"
date -u +%FT%TZ >> "$WORK/started_utc.txt"
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit.txt"
echo "$CEILING" > "$WORK/ceiling_seconds.txt"
RUN_ISP=("$PY" "$BD/scripts/run_isp.py" --manifest "$MANIFEST" --host-map "$E2/controls/pre_gpu_host_map.csv"
         --token-dict "$GF/geneformer/token_dictionary_gc104M.pkl" --bf16-bench "$BF" --scripts "$BD/scripts")
setsid env PYTHONPATH="$GF" "${RUN_ISP[@]}" --host "$HOST" --genes "$GENES" --out "$WORK" \
  > "$WORK/run_isp_stdout.log" 2>&1 &
ISP_PID=$!
ISP_PGID=""
for _ in $(seq 1 50); do
  CAND=$(ps -o pgid= -p "$ISP_PID" 2>/dev/null | tr -d ' ')
  if [ -n "$CAND" ] && [ "$CAND" = "$ISP_PID" ]; then ISP_PGID="$CAND"; break; fi
  sleep 0.1
done
[ -n "$ISP_PGID" ] || { kill -KILL "$ISP_PID" 2>/dev/null || true; echo "pgid never matched pid"; exit 5; }
echo "isp_pid=$ISP_PID isp_pgid=$ISP_PGID" > "$WORK/pgid.txt"
STOP_REASON="completed"
while kill -0 "$ISP_PID" 2>/dev/null; do
  sleep 60
  SPENT=$("$PY" -c "
import glob, json
s = 0.0
for f in glob.glob('$WORK/*/*/*.complete.json'):
    try: s += json.load(open(f)).get('seconds', 0.0)
    except Exception: pass
print(s)")
  echo "$(date -u +%FT%TZ) spent_s=$SPENT ceiling_s=$CEILING" >> "$WORK/ceiling_monitor.log"
  if "$PY" -c "import sys; sys.exit(0 if $SPENT >= $CEILING else 1)"; then
    echo "$(date -u +%FT%TZ) CEILING REACHED" >> "$WORK/ceiling_monitor.log"
    kill -TERM -"$ISP_PGID" 2>/dev/null || true; sleep 5; kill -KILL -"$ISP_PGID" 2>/dev/null || true
    STOP_REASON="ceiling"; break
  fi
done
wait "$ISP_PID" 2>/dev/null || true
date -u +%FT%TZ >> "$WORK/finished_utc.txt"
echo "$STOP_REASON" > "$WORK/stop_reason.txt"
NOOP_DIR="$DATA/isp/noop_spotchecks"; mkdir -p "$NOOP_DIR"; : > "$NOOP_DIR/results.jsonl"
FIRST_DONOR=$("$PY" -c "import json; print(sorted(json.load(open('$MANIFEST')))[0])")
BASE_ARGS=$("$PY" -c "import json,sys; print(json.dumps(sys.argv[1:]))" "${RUN_ISP[@]}")
for g in $("$PY" -c "import json; print(' '.join(json.load(open('$ORDER'))['noop_spotcheck_genes']))"); do
  if [ -f "$WORK/delete/$g/$FIRST_DONOR.complete.json" ] && [ -f "$WORK/overexpress/$g/$FIRST_DONOR.complete.json" ]; then
    PYTHONPATH="$GF" "$PY" "$BD/scripts/noop_spotcheck.py" --gene "$g" --donor "$FIRST_DONOR" --host "$HOST" \
      --run-isp-args "$BASE_ARGS" --work "$NOOP_DIR/${g}_work" --out "$NOOP_DIR/${g}_result.json" \
      >> "$NOOP_DIR/run.log" 2>&1 && cat "$NOOP_DIR/${g}_result.json" >> "$NOOP_DIR/results.jsonl" \
      || echo "{\"gene\":\"$g\",\"pass\":false,\"error\":\"noop_spotcheck.py failed\"}" >> "$NOOP_DIR/results.jsonl"
  fi
done
echo "E2_ISP_DONE reason=$STOP_REASON"
