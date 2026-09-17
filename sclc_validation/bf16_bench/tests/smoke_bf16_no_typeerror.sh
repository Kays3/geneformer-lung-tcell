#!/usr/bin/env bash
# Regression test (needs ts1: GPU + the pinned Geneformer env + the bf16
# patch already applied). NOT part of the pytest suite -- real GPU time; run
# only after checkpoint 2 (diff review) clears GPU time for this task.
#
# Proves the 4-site .float() patch actually fixes the documented failure:
# runs one T4 unit with --dtype bf16 against a PATCHED checkout and confirms
# it completes (no `TypeError: Got unsupported ScalarType BFloat16`) and
# produces a valid completion marker. Run smoke_fp32_bitwise_identical.sh
# first so the checkout is already patched.
#
# Usage (on ts1, after the patch is applied):
#   GENEFORMER_ROOT=/home/kaisar/workspace/geneformer-uv-starter/Geneformer \
#   ./smoke_bf16_no_typeerror.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$HERE")"
REPO_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"
T4_RUNNER="$REPO_ROOT/sclc_validation/immune_axis_test/run_t4_overexpression.py"
GENEFORMER_ROOT="${GENEFORMER_ROOT:?set GENEFORMER_ROOT to the patched Geneformer checkout}"
ITEM="${SMOKE_ITEM:-program_exhaustion}"
SOURCE="${SMOKE_SOURCE:-normal}"

if ! "$BENCH_DIR/apply_bf16_patch.sh" "$GENEFORMER_ROOT" | grep -qE "Applied|already applied"; then
    echo "error: patch is not applied to $GENEFORMER_ROOT -- run smoke_fp32_bitwise_identical.sh first." >&2
    exit 1
fi

run_dir="$BENCH_DIR/runs/_smoke/bf16"
rm -rf "$run_dir"

echo "=== bf16 smoke: patched checkout, --dtype bf16 ==="
set +e
GENEFORMER_ROOT="$GENEFORMER_ROOT" python3 "$T4_RUNNER" \
    --run-dir "$run_dir" --phase program --item "$ITEM" --source "$SOURCE" --dtype bf16 \
    2>&1 | tee "$run_dir.log"
rc=$?
set -e

if grep -q "TypeError: Got unsupported ScalarType BFloat16" "$run_dir.log" 2>/dev/null; then
    echo "FAIL: the documented numpy/bf16 TypeError still occurs -- patch did not take effect." >&2
    exit 1
fi
if [[ $rc -ne 0 ]]; then
    echo "FAIL: run exited $rc for a reason other than the known TypeError -- see $run_dir.log" >&2
    exit 1
fi

marker=$(python3 -c "import glob; print(glob.glob('$run_dir/markers/**/*.complete.json', recursive=True)[0])")
python3 -c "
import json
d = json.load(open('$marker'))
assert d['dtype'] == 'bf16', d
print('PASS: bf16 unit completed cleanly.')
print('  marker:', '$marker')
print('  elapsed_seconds:', d['elapsed_seconds'])
print('  peak_gpu_mem_gib:', d.get('peak_gpu_mem_gib'))
"
