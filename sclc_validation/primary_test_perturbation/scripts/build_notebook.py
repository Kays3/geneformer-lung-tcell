#!/usr/bin/env python3
"""Generate the reproducible primary-results review notebook."""

from pathlib import Path
import json

HERE = Path(__file__).resolve().parents[1]
def markdown(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(True), "attachments": {}}


def code(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source.splitlines(True)}


nb = {"cells": [
    markdown("# Primary donor-held-out SCLC perturbation analysis\n\nThis notebook reviews only test-cell perturbations with training-only Geneformer V2 reference centroids. It must not be used to pool training, evaluation, or whole-cohort results into the primary estimate."),
    code("from pathlib import Path\nimport pandas as pd\nHERE = Path.cwd()\nTABLES = HERE / 'tables'\nsummary = pd.read_csv(TABLES / 'primary_arm_summary.csv')\naudit = pd.read_csv(TABLES / 'coverage_audit.csv')\nhits = pd.read_csv(TABLES / 'primary_concordant_hits.csv')\nsummary"),
    markdown("## Completeness checks\n\nAll six directional comparisons are required for each arm. Missing files are pending computation, not zero-effect findings."),
    code("assert set(summary.arm) <= {'delete', 'overexpress'}\nassert len(audit) == 12\nprint(audit[['arm', 'comparison', 'status', 'rows']].to_string(index=False))\nprint('Complete tables:', (audit.status == 'ok').sum(), '/', len(audit))"),
    markdown("## Primary concordance\n\nA hit requires FDR < 0.05 in both arms, opposite-signed shifts, and at least 25 deletion detections. These filters do not replace donor-level consistency or biological validation."),
    code("hits.sort_values(['comparison', 'abs_delete_shift'], ascending=[True, False]).head(50)"),
    markdown("## Interpretation boundary\n\nGeneformer deletion and overexpression are rank-encoding interventions. They measure model representation sensitivity and should be described as candidate prioritization, not as validated gene knockout effects."),
], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
output = HERE / 'primary_test_perturbation.ipynb'
output.write_text(json.dumps(nb, indent=1) + "\n")
print(f'Wrote {output}')
