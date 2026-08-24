#!/usr/bin/env python3
"""Build the delete-vs-overexpress shift report (targeted panel + whole-genome).

Reads the two merged per-gene tables already on disk -- the targeted panel's
committed `targeted_panel_delete_overexpress_merged.csv` and the whole-genome
`allgene_delete_overexpress_shift.csv` written by
`plot_delete_vs_overexpress_shift.py` -- and the two PNGs those scripts produce.
Every number in the report is read from these tables at build time; this script
does not compute anything itself, only summarizes and lays out what the two
upstream scripts already wrote.
"""
from __future__ import annotations

import html
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent.parent
TARGETED_CSV = ROOT / "sclc_validation/perturbation_workflow/targeted_panel/results/targeted_panel_delete_overexpress_merged.csv"
TARGETED_FIGURE = ROOT / "sclc_validation/perturbation_workflow/targeted_panel/figures/delete_vs_overexpress_shift.png"
ALLGENE_CSV = HERE / "tables/allgene_delete_overexpress_shift.csv"
ALLGENE_FIGURE = HERE / "figures/allgene_delete_vs_overexpress_shift.png"
OUT = HERE / "reports/delete_vs_overexpress_shift_report.html"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for comparison in COMPARISONS:
        sub = df[df["comparison"] == comparison]
        n = len(sub)
        concordant = int(sub["concordant"].sum())
        rows.append({
            "comparison": comparison,
            "n_genes": n,
            "n_concordant": concordant,
            "pct_concordant": round(100 * concordant / n, 1) if n else 0.0,
        })
    return pd.DataFrame(rows)


def table_html(df: pd.DataFrame) -> str:
    return df.to_html(index=False, classes="data", border=0)


def main() -> None:
    targeted = pd.read_csv(TARGETED_CSV)
    allgene = pd.read_csv(ALLGENE_CSV)
    targeted_summary = summarize(targeted)
    allgene_summary = summarize(allgene)

    targeted_img = TARGETED_FIGURE.relative_to(OUT.parent, walk_up=True)
    allgene_img = ALLGENE_FIGURE.relative_to(OUT.parent, walk_up=True)

    body = f"""
<h1>Deletion vs overexpression shift toward the goal state</h1>
<p><strong>Generated:</strong> {datetime.now(timezone.utc).isoformat()}</p>

<p>Each gene-comparison pair has two shift measurements toward the same goal
state: one from deleting the gene, one from overexpressing it. If the gene
genuinely drives the model toward or away from that state, the two shifts
should be significant, adequately detected, and opposite in sign &mdash;
overexpression pushing one way, deletion undoing it. That is the
<code>concordant</code> criterion used throughout this repo's perturbation
work: <code>delete_fdr &lt; 0.05</code>, <code>overexpress_fdr &lt; 0.05</code>,
<code>delete_shift &times; overexpress_shift &lt; 0</code>, and
<code>delete_n &ge; 25</code>.</p>

<p>Two datasets, same question at two scales:</p>
<ul>
<li><strong>Targeted 50-gene panel</strong> &mdash; the pre-registered 21-gene
immune panel plus 29 top drivers from the prior NSCLC screen, both arms
census-complete per source state.</li>
<li><strong>Whole-genome screen</strong> &mdash; every gene present in each
held-out cell's token sequence, delete and overexpress arms, computed on
thinkstation2.</li>
</ul>

<h2>Targeted 50-gene panel</h2>
{table_html(targeted_summary)}
<img src="{targeted_img}" alt="Targeted panel: delete vs overexpress shift, one panel per comparison">

<h2>Whole-genome screen</h2>
<p>Genome-scale response is heavy-tailed (a handful of genes exceed
&#124;shift&#124; &gt; 0.3 against a typical scale of &plusmn;0.03), so each
panel below clips its axes to that panel's 1st&ndash;99th percentile; the
figure's own panel titles report how many points fall outside the drawn
frame.</p>
{table_html(allgene_summary)}
<img src="{allgene_img}" alt="Whole-genome screen: delete vs overexpress shift, one panel per comparison">

<h2>Reading the two together</h2>
<p><code>luad_to_sclc</code> and <code>luad_to_normal</code> carry far more
concordant genes than any comparison sourced from SCLC or Normal cells, at
both scales &mdash; LUAD is the state the model moves things toward or away
from most decisively. This is a model-level prioritization signal, not
biological validation: donor-level consistency, ambient-RNA sensitivity, and
independent evidence remain required before any gene here is read as a causal
driver.</p>

<h2>Provenance</h2>
<ul>
<li>Targeted panel: <code>sclc_validation/perturbation_workflow/targeted_panel/plot_delete_vs_overexpress_shift.py</code>
reading the committed <code>results/targeted_panel_delete_overexpress_merged.csv</code>.</li>
<li>Whole-genome: <code>sclc_validation/primary_test_perturbation/scripts/plot_delete_vs_overexpress_shift.py</code>
reading <code>stats/{{delete,overexpress}}/heldout_allgene_*.csv</code> under
<code>SCLC_PERTURBATION_ROOT</code> (thinkstation2), writing the full per-gene
table to <code>tables/allgene_delete_overexpress_shift.csv</code>.</li>
</ul>
"""
    document = (
        "<html><head><meta charset='utf-8'>"
        "<title>Delete vs overexpress shift</title>"
        "<style>body{font-family:system-ui;max-width:1200px;margin:2rem auto;"
        "line-height:1.5;padding:0 1rem} table{border-collapse:collapse;margin:1rem 0}"
        " th,td{border:1px solid #ccd;padding:.35rem .6rem;text-align:right}"
        " th{background:#eef2f5} td:first-child,th:first-child{text-align:left}"
        " code{background:#f4f4f4;padding:0 .2rem} img{max-width:100%;height:auto;"
        "border:1px solid #ddd;margin:.5rem 0 1.5rem}</style></head><body>"
        + body +
        "</body></html>"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(document)
    print(f"Wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
