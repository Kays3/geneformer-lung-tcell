#!/usr/bin/env bash
# Amendment 5 driver: runs balanced_donor_luad/scripts/run_isp.py, UNCHANGED, against the frozen
# null_genes_100_estimable.json draw order (the 100 ESTIMABLE genes, Amendment 5 B1 -- NOT
# Amendment 4's original 200-gene null_genes.json, which has only 12 estimable), under a hard GPU-arm
# budget ceiling. Fixes applied after Stanley's 2026-09-29 gates: B1 (this file now points at the
# right gene list), P6 (PGID race), P7/B4 (hash-verified inputs incl. the new gene list, the frozen
# freeze, and donor_manifest.json + run_config.json), P8/B3 (no-op spot check every 20th gene, run
# with the pinned venv python, not bare python3), P9 (single-host claim lock).
#
# The ceiling is enforced by reading the SAME "seconds" field run_isp.py writes into every
# *.complete.json marker (ISP-STD-1 E.3: count GPU arms, not wall clock).
#
# Usage: run_phase8_null.sh <host> <repo_worktree> <data_root> <geneformer_root> <venv_python> <ceiling_seconds>
set -euo pipefail
HOST="$1"; REPO="$2"; DATA="$3"; GF="$4"; PY="$5"; CEILING="${6:-57600}"
BD="$REPO/balanced_donor_luad"
COMMITTED_NULL_GENES="$BD/phase8_null/null_genes_100_estimable.json"   # B1: the 100 ESTIMABLE genes, git-tracked
FROZEN="$BD/phase8_null/frozen_100_estimable.json"                     # B2: gate-identity source, git-tracked
DONOR_MANIFEST="$DATA/goals/donor_manifest.json"
NULLDIR="$DATA/phase8_null"; WORK="$NULLDIR/out"
mkdir -p "$WORK"

N_CHECK=$(python3 -c "import json; d=json.load(open('$COMMITTED_NULL_GENES')); assert d['n']==100 and len(d['draw_order'])==100, d['n']; print(d['n'])")
echo "gene list verified: N=$N_CHECK (must be 100)"

# --- P9: single-host claim, refuse if another host already claimed this run ---
CLAIM="$NULLDIR/claimed_by_host.txt"
if [ -f "$CLAIM" ] && [ "$(cat "$CLAIM")" != "$HOST" ]; then
  echo "REFUSING: $NULLDIR already claimed by host $(cat "$CLAIM"), this is $HOST. Amendment 5 registers a" \
       "single host (or an explicit gene-disjoint split, not this script as written). Aborting." >&2
  exit 3
fi
echo "$HOST" > "$CLAIM"

# --- P7: hash-check every pinned input before any GPU call; refuse on mismatch ---
sha256_of() { shasum -a 256 "$1" 2>/dev/null | awk '{print $1}' || sha256sum "$1" | awk '{print $1}'; }
declare -A EXPECT=(
  ["$BD/scripts/run_isp.py"]="f440ec9c98bd56ee2ea6425abad7f1705e43779762990265e27cde50acbf2910"
  ["$BD/scripts/analyse.py"]="ab7c552f38b29432445d0a85bef21553dc693fb5ebe3c75c1f26773c64628b36"
  ["$REPO/sclc_validation/bf16_bench/dtype_cast.py"]="ac7acc383786a23c0dfa97f23a21b4a3057cb8b538125b64c81cd5112ae3ec89"
  ["$BD/scripts/model_cache.py"]="3cbaaf8e333f39d26372703934578cd0d7f5348b70d86cc464d785259b31a0e8"
  ["$BD/scripts/inproc_map.py"]="7d419a36ead0bbf6ecda469f27470c4b571e0592f714ed19b3ffe31c0998f517"
  ["$GF/geneformer/token_dictionary_gc104M.pkl"]="67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f"
  ["$COMMITTED_NULL_GENES"]="246b17bf279ff1a2028bffc98681416dd3162b0497b61191c57946a7c5902642"
  ["$FROZEN"]="5bd018b039831d4d04815a5dbee2b3758b51aff031fe05affea2edb9cee08976"
  ["$DONOR_MANIFEST"]="3cdd83a648c67cbd71d522ee59f1173e33a717b1f21cd93d3453ee32af794623"
)
# B4, noted rather than added here: the ovx index for these 100 genes cannot be pinned pre-launch --
# it is built from Phase 8's own GPU output (build_ovx_index.py reads *.complete.json markers and
# pickle lengths that do not exist until after this run), so it is inherently a post-GPU artifact.
# Its integrity gate is build_ovx_index.py's own count checks (0 failures required, unchanged from
# Amendment 3h), and its hash is recorded in run_config.json AFTER it is built, before null_analysis.py
# reads it -- not in this pre-launch EXPECT. Fold-model and goal-centroid identity are pinned in
# Amendments 1-3 (phase4_results/, goals/) and unchanged by this Amendment.
RUN_CONFIG="$NULLDIR/run_config_$HOST.json"
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
           "note_c3": "the ~10 no-op spot-check calls (P8) run in a separate output directory "
                      "(noop_spotchecks_<host>/) AFTER the main run stops, and their GPU seconds are NOT "
                      "counted toward ceiling_seconds above -- disclosed per Stanley's C3, not fixed, since "
                      "the no-op overhead is small (~10 calls) relative to the ceiling."},
           open(out_path, "w"), indent=1)
if mismatches:
    print("HASH MISMATCH, refusing to launch:", mismatches, file=sys.stderr)
    sys.exit(4)
print("all pinned inputs verified:", list(recorded.keys()))
PYEOF

GENES=$(python3 -c "import json; print(','.join(json.load(open('$COMMITTED_NULL_GENES'))['draw_order']))")

nvidia-smi --query-gpu=timestamp,name,driver_version,utilization.gpu,memory.used,power.draw \
  --format=csv -l 60 >> "$WORK/nvidia_smi_samples_$HOST.csv" 2>&1 &
SAMPLER=$!
date -u +%FT%TZ >> "$WORK/started_utc_$HOST.txt"
git -C "$REPO" rev-parse HEAD >> "$WORK/code_commit_$HOST.txt"

BASE_RUN_ISP_ARGS=$(python3 -c "
import json
print(json.dumps(['$PY', '$BD/scripts/run_isp.py',
  '--manifest', '$DATA/goals/donor_manifest.json',
  '--host-map', '$BD/controls/pre_gpu_host_map.csv',
  '--token-dict', '$GF/geneformer/token_dictionary_gc104M.pkl',
  '--bf16-bench', '$REPO/sclc_validation/bf16_bench', '--scripts', '$BD/scripts']))
")

# --- P6: launch in its own session; wait until the child's own pgid matches its pid before trusting it ---
setsid env PYTHONPATH="$GF" "$PY" "$BD/scripts/run_isp.py" \
  --manifest "$DATA/goals/donor_manifest.json" --host-map "$BD/controls/pre_gpu_host_map.csv" \
  --host "$HOST" --genes "$GENES" \
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

NOOP_EVERY=20
GENE_IDX=0
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

# --- P8: no-op spot check, every NOOP_EVERY-th completed gene (first donor), after the main run stops ---
NOOP_DIR="$NULLDIR/noop_spotchecks_$HOST"
mkdir -p "$NOOP_DIR"
FIRST_DONOR=$(python3 -c "import json; print(sorted(json.load(open('$DATA/goals/donor_manifest.json')))[0])")
python3 -c "
import json
genes = '$GENES'.split(',')
print(','.join(genes[i] for i in range(len(genes)) if (i+1) % $NOOP_EVERY == 0))
" > "$NOOP_DIR/spotcheck_genes.txt"
NOOP_RESULTS="$NOOP_DIR/results.jsonl"
: > "$NOOP_RESULTS"
IFS=',' read -ra NOOP_GENES < "$NOOP_DIR/spotcheck_genes.txt"
for g in "${NOOP_GENES[@]}"; do
  [ -z "$g" ] && continue
  # only spot-check genes that actually completed both markers for the first donor
  MARK_D="$WORK/delete/$g/${FIRST_DONOR//\//_}.complete.json"
  MARK_O="$WORK/overexpress/$g/${FIRST_DONOR//\//_}.complete.json"
  if [ -f "$MARK_D" ] && [ -f "$MARK_O" ]; then
    RESULT_FILE="$NOOP_DIR/${g}_result.json"
    PYTHONPATH="$GF" "$PY" "$BD/scripts/noop_spotcheck.py" --gene "$g" --donor "$FIRST_DONOR" --host "$HOST" \
      --run-isp-args "$BASE_RUN_ISP_ARGS" --work "$NOOP_DIR/${g}_work" --out "$RESULT_FILE" \
      >> "$NOOP_DIR/run.log" 2>&1 && cat "$RESULT_FILE" >> "$NOOP_RESULTS" || echo "{\"gene\":\"$g\",\"pass\":false,\"error\":\"noop_spotcheck.py failed, see run.log\"}" >> "$NOOP_RESULTS"
  fi
done

echo "PHASE8_NULL_DONE reason=$STOP_REASON noop_spotchecks=$NOOP_DIR/results.jsonl"
