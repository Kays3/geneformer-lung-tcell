#!/usr/bin/env bash
# Regression test (needs ts1: GPU + the pinned Geneformer env). NOT part of
# the pytest suite -- requires real GPU time, so it is not run until the
# implementation diff is reviewed and GPU time is authorized (see the task's
# checkpoint 2: "I review your implementation diff BEFORE GPU time is
# spent").
#
# Proves the bf16 patch + --dtype fp32 wrapper are a true no-op for fp32:
# runs the same single T4 unit twice -- once against an unpatched Geneformer
# checkout (today's behaviour), once against the patched checkout with
# --dtype fp32 explicitly -- and diffs the two runs' raw_sha256 completion
# markers. A mismatch means the patch or the dtype wrapper changed fp32
# output, which would be a real regression, not a no-op.
#
# Usage (on ts1):
#   GENEFORMER_ROOT=/home/kaisar/workspace/geneformer-uv-starter/Geneformer \
#   ./smoke_fp32_bitwise_identical.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(dirname "$HERE")"
REPO_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"
T4_RUNNER="$REPO_ROOT/sclc_validation/immune_axis_test/run_t4_overexpression.py"
GENEFORMER_ROOT="${GENEFORMER_ROOT:?set GENEFORMER_ROOT to the pinned Geneformer checkout}"
ITEM="${SMOKE_ITEM:-program_exhaustion}"   # any single manifest item id
SOURCE="${SMOKE_SOURCE:-normal}"

run_dir_a="$BENCH_DIR/runs/_smoke/fp32_unpatched"
run_dir_b="$BENCH_DIR/runs/_smoke/fp32_patched_wrapped"
rm -rf "$run_dir_a" "$run_dir_b"

echo "=== A: unpatched checkout, no wrapper (today's behaviour) ==="
git -C "$GENEFORMER_ROOT" status --porcelain -- geneformer/emb_extractor.py geneformer/perturber_utils.py | grep -q . && {
    echo "error: $GENEFORMER_ROOT already has changes to the patch-target files; start clean." >&2; exit 1; }
GENEFORMER_ROOT="$GENEFORMER_ROOT" python3 "$T4_RUNNER" \
    --run-dir "$run_dir_a" --phase program --item "$ITEM" --source "$SOURCE" --dtype fp32

echo "=== apply patch ==="
EXPECTED_GENEFORMER_COMMIT="${EXPECTED_GENEFORMER_COMMIT:-f45a6c7}" "$BENCH_DIR/apply_bf16_patch.sh" "$GENEFORMER_ROOT"

echo "=== B: patched checkout, --dtype fp32 wrapper installed ==="
GENEFORMER_ROOT="$GENEFORMER_ROOT" python3 "$T4_RUNNER" \
    --run-dir "$run_dir_b" --phase program --item "$ITEM" --source "$SOURCE" --dtype fp32

sha_a=$(python3 -c "import json,glob; print(json.load(open(glob.glob('$run_dir_a/markers/**/*.complete.json', recursive=True)[0]))['raw_sha256'])")
sha_b=$(python3 -c "import json,glob; print(json.load(open(glob.glob('$run_dir_b/markers/**/*.complete.json', recursive=True)[0]))['raw_sha256'])")

echo "unpatched raw_sha256: $sha_a"
echo "patched+fp32 raw_sha256: $sha_b"
if [[ "$sha_a" == "$sha_b" ]]; then
    echo "PASS: bitwise identical -- patch + --dtype fp32 is a true no-op."
else
    echo "FAIL: fp32 output changed after patching -- investigate before trusting bf16 results." >&2
    exit 1
fi
