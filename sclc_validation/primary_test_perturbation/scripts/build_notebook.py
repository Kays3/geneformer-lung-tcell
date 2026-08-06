#!/usr/bin/env python3
"""Generate the reproducible primary-results review notebook."""

from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
FIGURES = HERE / 'figures'
FIGURES.mkdir(exist_ok=True)
def markdown(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(True), "attachments": {}}


def code(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source.splitlines(True)}


summary = pd.read_csv(HERE / 'tables' / 'primary_arm_summary.csv')
audit = pd.read_csv(HERE / 'tables' / 'coverage_audit.csv')
hits = pd.read_csv(HERE / 'tables' / 'primary_concordant_hits.csv')

# Save compact, reproducible figures alongside the notebook.
plot = summary.pivot(index='comparison', columns='arm', values='qualified_primary')
ax = plot[['delete', 'overexpress']].plot(kind='bar', figsize=(10, 5), color=['#3b6ea8', '#d97732'])
ax.set_title('Primary qualified genes by disease comparison')
ax.set_xlabel('Comparison')
ax.set_ylabel('Qualified primary genes')
ax.tick_params(axis='x', rotation=35)
ax.legend(title='Perturbation')
ax.grid(axis='y', alpha=0.25)
plt.tight_layout()
plt.savefig(FIGURES / 'qualified_genes_by_comparison.png', dpi=180)
plt.close()

heat = summary.pivot(index='comparison', columns='arm', values='qualified_primary')
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(heat[['delete', 'overexpress']].values, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(2), ['delete', 'overexpress'])
ax.set_yticks(range(len(heat.index)), heat.index)
ax.set_title('Qualified primary genes')
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        ax.text(j, i, f'{heat.iloc[i, j]:,}', ha='center', va='center', fontsize=9)
fig.colorbar(im, ax=ax, label='Genes')
fig.tight_layout()
fig.savefig(FIGURES / 'qualified_genes_heatmap.png', dpi=180)
plt.close(fig)

nb = {"cells": [
    markdown("# Primary donor-held-out SCLC perturbation analysis\n\nThis notebook reviews only test-cell perturbations with training-only Geneformer V2 reference centroids. It must not be used to pool training, evaluation, or whole-cohort results into the primary estimate."),
    code("from pathlib import Path\nfrom IPython.display import display, Image\nimport pandas as pd\nHERE = Path.cwd()\nTABLES = HERE / 'tables'\nFIGURES = HERE / 'figures'\nsummary = pd.read_csv(TABLES / 'primary_arm_summary.csv')\naudit = pd.read_csv(TABLES / 'coverage_audit.csv')\nhits = pd.read_csv(TABLES / 'primary_concordant_hits.csv')\ndisplay(summary)") ,
    markdown("## Completeness checks\n\nAll six directional comparisons are required for each arm. Missing files are pending computation, not zero-effect findings."),
    code("assert set(summary.arm) <= {'delete', 'overexpress'}\nassert len(audit) == 12\nassert (audit.status == 'ok').all()\nprint(audit[['arm', 'comparison', 'status', 'rows']].to_string(index=False))\nprint('Complete tables:', (audit.status == 'ok').sum(), '/', len(audit))\ndisplay(audit)"),
    markdown("## Figures\n\nThe figures below are generated from the final `primary_arm_summary.csv` table."),
    code("display(Image(filename=str(FIGURES / 'qualified_genes_by_comparison.png')))\ndisplay(Image(filename=str(FIGURES / 'qualified_genes_heatmap.png')))"),
    markdown("## Primary concordance\n\nA hit requires FDR < 0.05 in both arms, opposite-signed shifts, and at least 25 deletion detections. These filters do not replace donor-level consistency or biological validation."),
    code("top_hits = hits.sort_values(['comparison', 'abs_delete_shift'], ascending=[True, False]).head(50)\ndisplay(top_hits)") ,
    markdown("## Interpretation boundary\n\nGeneformer deletion and overexpression are rank-encoding interventions. They measure model representation sensitivity and should be described as candidate prioritization, not as validated gene knockout effects."),
], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
output = HERE / 'primary_test_perturbation.ipynb'
output.write_text(json.dumps(nb, indent=1) + "\n")
print(f'Wrote {output}')
