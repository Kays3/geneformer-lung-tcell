#!/usr/bin/env bash
# Amendment 6 driver: runs balanced_donor_luad/scripts/run_isp.py, UNCHANGED, against ONLY the 100
# NEW estimable genes (positions 985-1791 of the seed-20260929 order --
# null_genes_200_estimable.json's draw_order[100:200]). A5's markers for positions 1-984 are REUSED,
# never re-run (Amendment 6 s.6.2). Implements Stanley's run-package requirements (gate on commit
# 037683c's registration; this is the separate run-package gate):
#
#   R-B(i)   bridge check: 3 named (gene,donor) pairs already completed by A5, re-run fresh into a
#            scratch dir, and required to be bitwise-identical to A5's actual output pickles --
#            BEFORE any new gene is run. If they differ, A5 and A6 markers are not poolable and the
#            extension refuses to launch.
#   R-B(ii)  environment identity: goals/xfer_sha256.txt (the Amendments 1-3 pin covering every
#            per-donor goals.pkl, isp_input dataset, and phase4 fold-model file actually consumed by
#            run_isp.py) is re-verified in full via `sha256sum -c`.
#   R-B(iii) pip freeze + torch/CUDA/driver versions are recorded in run_config for comparison
#            against A5's nvidia-smi samples.
#   R-D      a SEPARATE output directory (phase8_null_ext/out, not phase8_null/out) so the ceiling
#            monitor's seconds-sum can never double-count A5's ~6.7 GPU-h against this extension's
#            10.20 GPU-h ceiling. null_analysis.py reads BOTH trees for the combined N=200 analysis.
#   R-E      no-op due positions are the OVERALL order positions 120,140,160,180,200 (new-list index
#            20,40,60,80,100 within these 100 new genes) -- distinct from A5's already-checked
#            overall positions 20,40,60,80,100.
#
# Stanley's run-package re-gate (2026-09-29T17-41-03Z-a6-runpkg-gate) added:
#   B2 goals/xfer_sha256.txt is itself pinned in EXPECT (a `sha256sum -c` pass proves nothing if the
#      manifest file itself can silently change).
#   B3 refuses unless A5 has actually finished: phase8_null/out/stop_reason_<host>.txt must exist AND
#      no run_isp.py process may be running (checked via `ps`, never `pgrep -f`) -- this is the very
#      first check, before anything else, so an unfinished A5 can never be bridged against and two GPU
#      runs can never overlap.
#   B4 noop_spotcheck.py (executed by this driver) and null_analysis_combined.py (the combined-analysis
#      code this run package exists to feed) are both pinned in EXPECT too.
#
# Usage: run_phase8_null_extension.sh <host> <repo_worktree> <data_root> <geneformer_root> <venv_python> <ceiling_seconds>
set -euo pipefail
HOST="$1"; REPO="$2"; DATA="$3"; GF="$4"; PY="$5"; CEILING="${6:-36720}"
BD="$REPO/balanced_donor_luad"

if [ "$HOST" != "thinkstation1" ]; then
  echo "REFUSING: Amendment 6 registers this extension ts1-only. Got HOST=$HOST." >&2
  exit 2
fi

A5_NULL_GENES="$BD/phase8_null/null_genes_100_estimable.json"
A5_FROZEN="$BD/phase8_null/frozen_100_estimable.json"
EXT_NULL_GENES_200="$BD/phase8_null/null_genes_200_estimable.json"
EXT_FROZEN_200="$BD/phase8_null/frozen_200_estimable.json"
DONOR_MANIFEST="$DATA/goals/donor_manifest.json"
XFER_MANIFEST="$DATA/goals/xfer_sha256.txt"

A5_DIR="$DATA/phase8_null"; A5_WORK="$A5_DIR/out"
EXTDIR="$DATA/phase8_null_ext"; WORK="$EXTDIR/out"
mkdir -p "$WORK"

# --- B3: refuse unless A5 has actually finished -- stop_reason file exists AND no run_isp.py process
#     is running (ps, never pgrep -f). This is the very first check: an unfinished A5 must never be
#     bridged against, and two GPU runs on this host must never overlap. ---
A5_STOP_REASON="$A5_WORK/stop_reason_$HOST.txt"
if [ ! -f "$A5_STOP_REASON" ]; then
  echo "REFUSING (B3): $A5_STOP_REASON does not exist -- A5 has not finished yet." >&2
  exit 1
fi
A5_STILL_RUNNING=$(ps -eo pid,cmd | grep -c "[r]un_isp\.py" || true)
if [ "$A5_STILL_RUNNING" != "0" ]; then
  echo "REFUSING (B3): a run_isp.py process is still running on $HOST -- A5 (or something else) is" \
       "not finished. Refusing to bridge or launch until it exits." >&2
  exit 1
fi
echo "B3 verified: A5 finished (stop_reason=$(cat "$A5_STOP_REASON")), no run_isp.py process running"

# --- gene list check: 200 total, first 100 identical to A5's committed list ---
python3 -c "
import json
full = json.load(open('$EXT_NULL_GENES_200'))
a5 = json.load(open('$A5_NULL_GENES'))
assert full['n'] == 200 and len(full['draw_order']) == 200, full['n']
assert full['draw_order'][:100] == a5['draw_order'], 'first 100 of the 200 must equal A5 exactly'
print('gene list verified: N=200, first 100 == A5 exactly')
"
NEW_GENES_CSV=$(python3 -c "import json; print(','.join(json.load(open('$EXT_NULL_GENES_200'))['draw_order'][100:200]))")

# --- single-host claim, its own lock (separate run tree from A5) ---
CLAIM="$EXTDIR/claimed_by_host.txt"
if [ -f "$CLAIM" ] && [ "$(cat "$CLAIM")" != "$HOST" ]; then
  echo "REFUSING: $EXTDIR already claimed by host $(cat "$CLAIM"), this is $HOST." >&2
  exit 3
fi
echo "$HOST" > "$CLAIM"

# --- hash-check every pinned script/data input: A5's full EXPECT list plus the new gene-list files ---
declare -A EXPECT=(
  ["$BD/scripts/run_isp.py"]="f440ec9c98bd56ee2ea6425abad7f1705e43779762990265e27cde50acbf2910"
  ["$BD/scripts/analyse.py"]="ab7c552f38b29432445d0a85bef21553dc693fb5ebe3c75c1f26773c64628b36"
  ["$REPO/sclc_validation/bf16_bench/dtype_cast.py"]="ac7acc383786a23c0dfa97f23a21b4a3057cb8b538125b64c81cd5112ae3ec89"
  ["$BD/scripts/model_cache.py"]="3cbaaf8e333f39d26372703934578cd0d7f5348b70d86cc464d785259b31a0e8"
  ["$BD/scripts/inproc_map.py"]="7d419a36ead0bbf6ecda469f27470c4b571e0592f714ed19b3ffe31c0998f517"
  ["$GF/geneformer/token_dictionary_gc104M.pkl"]="67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f"
  ["$DONOR_MANIFEST"]="3cdd83a648c67cbd71d522ee59f1173e33a717b1f21cd93d3453ee32af794623"
  ["$A5_NULL_GENES"]="246b17bf279ff1a2028bffc98681416dd3162b0497b61191c57946a7c5902642"
  ["$A5_FROZEN"]="5bd018b039831d4d04815a5dbee2b3758b51aff031fe05affea2edb9cee08976"
  ["$EXT_NULL_GENES_200"]="9c3ad015e56c90277401e83f892eea2cf9127071ac24eaaa0b76789bf44c0f02"
  ["$EXT_FROZEN_200"]="da7bf511a2016b153132ad2530a6b8564b5ec2d50e6098c3f7263250161bd919"
  ["$XFER_MANIFEST"]="91f22298e7525ca957e72b837f9e6d67a0c29487c62c30e94dca9c1058c0b46a"
  ["$BD/scripts/noop_spotcheck.py"]="e9bbdfcc7e8797e437d1687e9c74c5ce49de8c57ecdfce9e04daefcffc13d5be"
  # Stanley's run-package re-gate-2 (B5/B6/m1): null_analysis.py's apply_noop_gate fix changed its hash;
  # it was NOT previously in this EXPECT list even though null_analysis_combined.py imports it -- added
  # per Stanley's "its EXPECT must pin both analysis files' new hashes".
  ["$BD/scripts/null_analysis.py"]="bf12d145ead61d11ae71786f76c64a7260967a916e39a4faeeca88290f1b5c9c"
  ["$BD/scripts/null_analysis_combined.py"]="6a372f2f7fbd95fb217849f6f10f08856314512090e5a0f527e56116fd9be949"
)
RUN_CONFIG="$EXTDIR/run_config_$HOST.json"
python3 - "$RUN_CONFIG" "$HOST" "$CEILING" "${!EXPECT[@]}" <<'PYEOF' "${EXPECT[@]}"
import hashlib, json, sys
out_path, host, ceiling = sys.argv[1], sys.argv[2], sys.argv[3]
n = (len(sys.argv) - 4) // 2
paths = sys.argv[4:4+n]
expected = sys.argv[4+n:4+2*n]
mismatches = []
recorded = {}
for path, exp in zip(paths, expected):
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    recorded[path] = {"expected": exp, "actual": h, "match": h == exp}
    if h != exp:
        mismatches.append(path)
json.dump({"host": host, "ceiling_seconds": float(ceiling), "hashes": recorded, "mismatches": mismatches,
           "note_rd": "the ceiling monitor below globs $WORK (phase8_null_ext/out) only -- A5's "
                      "phase8_null/out tree is never summed, by construction of the separate directory.",
           "note_re": "no-op due positions for this run are OVERALL order 120,140,160,180,200 "
                      "(new-list index 20,40,60,80,100), distinct from A5's already-checked 20,40,60,80,100.",
           "note_ceiling_scope": "Stanley's minor note: the bridge check (3 calls) and the no-op spot "
                      "checks (~5 calls) run outside this ceiling_seconds figure, same disclosed pattern "
                      "as A5's C3/note_c3 -- their GPU seconds are small relative to the ceiling and are "
                      "not counted toward it."},
           open(out_path, "w"), indent=1)
if mismatches:
    print("HASH MISMATCH, refusing to launch:", mismatches, file=sys.stderr)
    sys.exit(4)
print("all pinned script/data inputs verified:", list(recorded.keys()))
PYEOF

# --- R-B(ii): re-verify the FULL Amendments 1-3 environment pin (goals.pkl, isp_input datasets, fold models) ---
HOME_REL_MANIFEST="workspace/$(basename "$DATA")/goals/xfer_sha256.txt"
( cd "$HOME" && sha256sum -c "$HOME_REL_MANIFEST" --quiet > "$EXTDIR/xfer_verify_$HOST.log" 2>&1 )
if [ -s "$EXTDIR/xfer_verify_$HOST.log" ]; then
  echo "R-B(ii) FAILED: goals/xfer_sha256.txt manifest mismatch, see $EXTDIR/xfer_verify_$HOST.log" >&2
  exit 6
fi
echo "R-B(ii) verified: full Amendments 1-3 environment pin (198 files: goals + fold models) unchanged"

# --- R-B(iii): record the software environment for comparison against A5's ---
"$PY" -m pip freeze > "$EXTDIR/pip_freeze_$HOST.txt"
"$PY" -c "import torch; print('torch', torch.__version__); print('cuda', torch.version.cuda)" > "$EXTDIR/torch_cuda_$HOST.txt" 2>&1 || true
nvidia-smi --query-gpu=driver_version --format=csv,noheader > "$EXTDIR/driver_version_$HOST.txt"

# --- R-B(i): bridge check -- 3 named A5 (gene,donor) pairs, both ops each, bitwise vs A5's actual output ---
# Named per Stanley's suggestion: first estimable gene, one mid-list gene, one with the highest n_token_cells
# (computed read-only from frozen_100_estimable.json: position 1, position 50, position 6 respectively).
BRIDGE_DIR="$EXTDIR/bridge_check_$HOST"; mkdir -p "$BRIDGE_DIR"
declare -a BRIDGE_GENES=("ENSG00000115009" "ENSG00000104979" "ENSG00000162244")
declare -a BRIDGE_DONORS=("Leader_Merad_2021_570" "Kim_Lee_2020_P0006" "Kim_Lee_2020_P0006")
: > "$EXTDIR/bridge_check_$HOST.log"
for i in 0 1 2; do
  g="${BRIDGE_GENES[$i]}"; d="${BRIDGE_DONORS[$i]}"
  OUTDIR="$BRIDGE_DIR/pair_$i"
  mkdir -p "$OUTDIR"
  PYTHONPATH="$GF" "$PY" "$BD/scripts/run_isp.py" \
    --manifest "$DONOR_MANIFEST" --host-map "$BD/controls/pre_gpu_host_map.csv" --host "$HOST" \
    --genes "$g" --donors "$d" \
    --token-dict "$GF/geneformer/token_dictionary_gc104M.pkl" --bf16-bench "$REPO/sclc_validation/bf16_bench" \
    --scripts "$BD/scripts" --out "$OUTDIR" > "$BRIDGE_DIR/pair_${i}.log" 2>&1
  for op in delete overexpress; do
    A5_PICKLE=$(ls "$A5_WORK/$op/$g/in_silico_${op}_${d}"*_raw.pickle 2>/dev/null | head -1)
    NEW_PICKLE=$(ls "$OUTDIR/$op/$g/in_silico_${op}_${d}"*_raw.pickle 2>/dev/null | head -1)
    if [ -z "$A5_PICKLE" ] || [ -z "$NEW_PICKLE" ]; then
      echo "BRIDGE CHECK FAILED: missing pickle for $op $g $d (A5='$A5_PICKLE' NEW='$NEW_PICKLE')" >&2
      exit 7
    fi
    if ! cmp -s "$A5_PICKLE" "$NEW_PICKLE"; then
      echo "BRIDGE CHECK FAILED: $op $g $d differs bitwise between A5 and a fresh re-run -- environments" \
           "are NOT poolable, refusing to launch the extension." >&2
      exit 7
    fi
    echo "bridge check OK: $op $g $d bitwise-identical ($A5_PICKLE)" >> "$EXTDIR/bridge_check_$HOST.log"
  done
done
echo "R-B(i) bridge check PASSED: 6/6 (gene,donor,op) calls bitwise-identical to A5"

# --- main run: only the 100 NEW genes, into the separate WORK tree ---
nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$WORK/nvidia_smi_samples_$HOST.csv" 2>&1 &
SAMPLER=$!
date -u +%FT%TZ >> "$WORK/started_utc_$HOST.txt"
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit_$HOST.txt"

BASE_RUN_ISP_ARGS=$(python3 -c "
import json
print(json.dumps(['$PY', '$BD/scripts/run_isp.py',
  '--manifest', '$DONOR_MANIFEST',
  '--host-map', '$BD/controls/pre_gpu_host_map.csv',
  '--token-dict', '$GF/geneformer/token_dictionary_gc104M.pkl',
  '--bf16-bench', '$REPO/sclc_validation/bf16_bench', '--scripts', '$BD/scripts']))
")

setsid env PYTHONPATH="$GF" "$PY" "$BD/scripts/run_isp.py" \
  --manifest "$DONOR_MANIFEST" --host-map "$BD/controls/pre_gpu_host_map.csv" \
  --host "$HOST" --genes "$NEW_GENES_CSV" \
  --token-dict "$GF/geneformer/token_dictionary_gc104M.pkl" --bf16-bench "$REPO/sclc_validation/bf16_bench" \
  --scripts "$BD/scripts" --out "$WORK" > "$WORK/run_log_stdout_$HOST.log" 2>&1 &
ISP_PID=$!

ISP_PGID=""
for _ in $(seq 1 50); do
  CAND=$(ps -o pgid= -p "$ISP_PID" 2>/dev/null | tr -d ' ')
  if [ -n "$CAND" ] && [ "$CAND" = "$ISP_PID" ]; then ISP_PGID="$CAND"; break; fi
  sleep 0.1
done
if [ -z "$ISP_PGID" ]; then
  echo "PGID never matched PID after 5s (got '$CAND') -- refusing to arm the ceiling kill against the wrong" \
       "process group. Killing the child directly and aborting." >&2
  kill -KILL "$ISP_PID" 2>/dev/null || true
  exit 5
fi
echo "isp_pid=$ISP_PID isp_pgid=$ISP_PGID (verified pgid==pid)" | tee "$WORK/pgid_$HOST.txt"

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
  REAL=$(ps -o pid,pgid,cmd -p "$ISP_PID" 2>/dev/null | tail -n +2 | grep -c "run_isp.py" || true)
  echo "$(date -u +%FT%TZ) spent_s=$SPENT ceiling_s=$CEILING pid_verified=$REAL" >> "$WORK/ceiling_monitor_$HOST.log"
  if [ "$REAL" = "0" ]; then break; fi
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

# --- R-E: no-op spot check on the NEW genes only. Due = overall positions 120,140,160,180,200
#     (new-list index 20,40,60,80,100 within these 100 genes) -- distinct from A5's 20,40,60,80,100. ---
NOOP_EVERY=20
NOOP_DIR="$EXTDIR/noop_spotchecks_$HOST"
mkdir -p "$NOOP_DIR"
FIRST_DONOR=$(python3 -c "import json; print(sorted(json.load(open('$DONOR_MANIFEST')))[0])")
python3 -c "
genes = '$NEW_GENES_CSV'.split(',')
print(','.join(genes[i] for i in range(len(genes)) if (i+1) % $NOOP_EVERY == 0))
" > "$NOOP_DIR/spotcheck_genes.txt"
NOOP_RESULTS="$NOOP_DIR/results.jsonl"
: > "$NOOP_RESULTS"
IFS=',' read -ra NOOP_GENES < "$NOOP_DIR/spotcheck_genes.txt"
for g in "${NOOP_GENES[@]}"; do
  [ -z "$g" ] && continue
  MARK_D="$WORK/delete/$g/${FIRST_DONOR}.complete.json"
  MARK_O="$WORK/overexpress/$g/${FIRST_DONOR}.complete.json"
  if [ -f "$MARK_D" ] && [ -f "$MARK_O" ]; then
    RESULT_FILE="$NOOP_DIR/${g}_result.json"
    PYTHONPATH="$GF" "$PY" "$BD/scripts/noop_spotcheck.py" --gene "$g" --donor "$FIRST_DONOR" --host "$HOST" \
      --run-isp-args "$BASE_RUN_ISP_ARGS" --work "$NOOP_DIR/${g}_work" --out "$RESULT_FILE" \
      >> "$NOOP_DIR/run.log" 2>&1 && cat "$RESULT_FILE" >> "$NOOP_RESULTS" || echo "{\"gene\":\"$g\",\"pass\":false,\"error\":\"noop_spotcheck.py failed, see run.log\"}" >> "$NOOP_RESULTS"
  fi
done

echo "PHASE8_NULL_EXT_DONE reason=$STOP_REASON noop_spotchecks=$NOOP_RESULTS"
