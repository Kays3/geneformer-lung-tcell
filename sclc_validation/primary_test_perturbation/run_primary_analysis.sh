#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${PYTHON_BIN:-python}

"$PYTHON_BIN" "$HERE/scripts/build_primary_report.py"
"$PYTHON_BIN" "$HERE/scripts/build_notebook.py"

printf '%s\n' "Primary report: $HERE/reports/primary_test_perturbation_report.html"
printf '%s\n' "Notebook:       $HERE/primary_test_perturbation.ipynb"
