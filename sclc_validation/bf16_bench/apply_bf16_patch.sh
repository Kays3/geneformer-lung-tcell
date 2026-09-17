#!/usr/bin/env bash
# Apply the 4-site .float() bf16-export fix to a Geneformer checkout.
#
# The fix upcasts bfloat16 tensors to float32 immediately before the four
# `.cpu().numpy()` calls that would otherwise raise
# `TypeError: Got unsupported ScalarType BFloat16` (numpy has no bf16 dtype).
# It is a no-op for fp32 tensors -- see tests/test_bf16_patch.py, which
# proves fp32 markers are bitwise-identical with and without this patch.
#
# Usage:
#   GENEFORMER_ROOT=/path/to/Geneformer ./apply_bf16_patch.sh
#   ./apply_bf16_patch.sh /path/to/Geneformer
#
# Verifies the checkout is at the expected pinned commit before patching
# (override with EXPECTED_GENEFORMER_COMMIT if the pin changes). Idempotent:
# a second run detects the patch is already applied and exits 0 without
# re-patching.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH_FILE="${HERE}/geneformer-bf16-export.patch"
GENEFORMER_ROOT="${1:-${GENEFORMER_ROOT:-}}"
EXPECTED_GENEFORMER_COMMIT="${EXPECTED_GENEFORMER_COMMIT:-f45a6c7}"

if [[ -z "${GENEFORMER_ROOT}" ]]; then
    echo "usage: $0 <path-to-Geneformer-checkout>  (or set GENEFORMER_ROOT)" >&2
    exit 2
fi
if [[ ! -d "${GENEFORMER_ROOT}/.git" ]]; then
    echo "error: ${GENEFORMER_ROOT} is not a git checkout" >&2
    exit 2
fi

actual_commit="$(git -C "${GENEFORMER_ROOT}" rev-parse --short HEAD)"
if [[ "${actual_commit}" != ${EXPECTED_GENEFORMER_COMMIT}* ]]; then
    echo "error: ${GENEFORMER_ROOT} is at ${actual_commit}, expected ${EXPECTED_GENEFORMER_COMMIT}." >&2
    echo "Set EXPECTED_GENEFORMER_COMMIT to override if the pin is intentionally different." >&2
    exit 1
fi

if patch -p1 -d "${GENEFORMER_ROOT}" --dry-run --forward -s < "${PATCH_FILE}" >/dev/null 2>&1; then
    patch -p1 -d "${GENEFORMER_ROOT}" < "${PATCH_FILE}"
    echo "Applied bf16 export patch to ${GENEFORMER_ROOT} (commit ${actual_commit})."
elif patch -p1 -d "${GENEFORMER_ROOT}" --dry-run --forward -R -s < "${PATCH_FILE}" >/dev/null 2>&1; then
    echo "bf16 export patch already applied to ${GENEFORMER_ROOT} (commit ${actual_commit}); nothing to do."
else
    echo "error: patch does not apply cleanly (forward or reverse) to ${GENEFORMER_ROOT}." >&2
    echo "The checkout may have diverged from the expected pin in the two touched files." >&2
    exit 1
fi
