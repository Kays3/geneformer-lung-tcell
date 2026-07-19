#!/usr/bin/env python3
"""Create and execute the reproducible perturbation-analysis notebook."""

from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient


HERE = Path(__file__).resolve().parent
notebook = nbf.v4.new_notebook()
notebook["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3"},
}
notebook["cells"] = [
    nbf.v4.new_markdown_cell(
        "# Held-out all-gene perturbation statistics and pathway analysis\n\n"
        "**Objective.** Audit the six completed directional Geneformer goal-state result tables, "
        "separate FDR significance from positive goal movement, rank coverage-qualified genes, "
        "and review enrichment with explicit biological caveats.\n\n"
        "**Primary ranking definition.** `Goal_end_FDR < 0.05`, `Shift_to_goal_end > 0`, and "
        "`N_Detections >= 25`."
    ),
    nbf.v4.new_code_cell(
        "from pathlib import Path\n"
        "import pandas as pd\n"
        "from IPython.display import Image, display\n\n"
        "HERE = Path.cwd()\n"
        "TABLES = HERE / 'source_tables'\n"
        "summary = pd.read_csv(TABLES / 'comparison_summary.csv')\n"
        "top_genes = pd.read_csv(TABLES / 'top_goal_shift_genes.csv')\n"
        "pathways = pd.read_csv(TABLES / 'pathway_enrichment.csv')\n"
        "themes = pd.read_csv(TABLES / 'pathway_theme_summary.csv')\n"
        "summary"
    ),
    nbf.v4.new_markdown_cell(
        "## Data-quality checks\n\n"
        "The source audit below checks row counts and missing values in the derived tables. "
        "The build script separately verifies one gene row per input comparison."
    ),
    nbf.v4.new_code_cell(
        "quality = pd.DataFrame({\n"
        "    'table': ['comparison_summary', 'top_goal_shift_genes', 'pathway_enrichment', 'pathway_theme_summary'],\n"
        "    'rows': [len(summary), len(top_genes), len(pathways), len(themes)],\n"
        "    'missing_cells': [summary.isna().sum().sum(), top_genes.isna().sum().sum(), pathways.isna().sum().sum(), themes.isna().sum().sum()],\n"
        "})\n"
        "assert len(summary) == 6\n"
        "assert summary['genes_tested'].between(11_000, 15_000).all()\n"
        "assert (summary['qualified_toward_goal'] <= summary['significant_toward_goal']).all()\n"
        "quality"
    ),
    nbf.v4.new_markdown_cell(
        "## Statistical results\n\n"
        "FDR significance alone includes both positive and negative centroid shifts. The qualified "
        "set adds direction and coverage requirements so the biological ranking cannot silently mix "
        "movement toward and away from the requested goal."
    ),
    nbf.v4.new_code_cell(
        "summary_view = summary.assign(\n"
        "    qualified_percent=lambda x: (100 * x.qualified_fraction).round(2),\n"
        ")[['comparison_label', 'genes_tested', 'fdr_significant', 'significant_toward_goal',\n"
        "   'significant_away_from_goal', 'qualified_toward_goal', 'qualified_percent',\n"
        "   'max_positive_shift', 'median_qualified_shift']]\n"
        "summary_view"
    ),
    nbf.v4.new_code_cell(
        "display(Image(filename=str(HERE / 'figures' / 'goal_shift_top_genes.png')))"
    ),
    nbf.v4.new_markdown_cell(
        "## Top ranked genes by direction\n\n"
        "The table retains effect size, FDR, detection count, rank, and stable Ensembl identifier. "
        "A positive shift means movement toward a centroid after deletion—not increased expression "
        "in the goal state and not evidence of therapeutic benefit."
    ),
    nbf.v4.new_code_cell(
        "top_genes.query('rank <= 10')[['comparison_label', 'rank', 'Gene_name', 'Shift_to_goal_end', 'Goal_end_FDR', 'N_Detections']]"
    ),
    nbf.v4.new_markdown_cell(
        "## Pathway analysis\n\n"
        "g:Profiler enrichment used the tested genes in each comparison as the custom background. "
        "Reactome, KEGG, and GO Biological Process terms were FDR corrected by g:Profiler. The "
        "theme table consolidates redundant significant terms with predefined text patterns; it is "
        "descriptive rather than a new inferential layer."
    ),
    nbf.v4.new_code_cell(
        "themes[['comparison_label', 'theme', 'best_adjusted_p_value', 'matching_term_count',\n"
        "        'representative_source', 'representative_term', 'intersection_size']]"
    ),
    nbf.v4.new_code_cell(
        "display(Image(filename=str(HERE / 'figures' / 'pathway_enrichment.png')))"
    ),
    nbf.v4.new_markdown_cell(
        "## Biological interpretation and decision\n\n"
        "- **Recurring programs:** translation/ribosome and immune/cytokine terms dominate; oxidative "
        "phosphorylation is strongest in LUAD → LUSC and also appears in NORMAL → LUSC.\n"
        "- **Cell-identity warning:** SFTPB/SFTPC/NAPSA, MUC1/PIGR/keratins, and FBLN1/DCN/ACKR1 "
        "are difficult to reconcile as purely intrinsic T-cell programs. Ambient RNA, doublets, and "
        "tissue composition are plausible drivers.\n"
        "- **Decision:** share as a prioritization screen with caveats; do not nominate causal targets "
        "until donor-stratified stability and decontamination sensitivity are demonstrated.\n\n"
        "### Next analyses\n\n"
        "1. Compute per-donor effects and leave-one-donor-out ranks.\n"
        "2. Repeat after ambient-RNA correction and doublet/high-epithelial-burden exclusion.\n"
        "3. Compare perturbation direction with source-versus-goal differential expression.\n"
        "4. Validate prioritized T-cell-specific genes in an independent donor-balanced cohort."
    ),
]

output = HERE / "perturbation_statistics.ipynb"
client = NotebookClient(notebook, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(HERE)}})
client.execute()
nbf.write(notebook, output)
print(f"Wrote and executed {output}")
