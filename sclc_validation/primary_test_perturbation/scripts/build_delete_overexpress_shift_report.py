#!/usr/bin/env python3
"""Build the shift report: goal-vs-alt specificity, and delete-vs-overexpress
concordance, each at targeted-panel and whole-genome scale.

Reads the merged per-gene tables already on disk -- the targeted panel's
committed `targeted_panel_delete_overexpress_merged.csv`, its per-arm
per-comparison `results/{arm}/targeted_{arm}_{comparison}.csv` tables, and the
whole-genome `allgene_delete_overexpress_shift.csv` / `allgene_goal_vs_alt_shift.csv`
written by the two `plot_*.py` scripts in this repo -- plus the PNGs those
scripts produce. Every number in the report is read from these tables at
build time; this script does not compute anything beyond summary statistics
(counts, Pearson r) from what the upstream scripts already wrote.
"""
from __future__ import annotations

import base64
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent.parent
TARGETED_DIR = ROOT / "sclc_validation/perturbation_workflow/targeted_panel"
TARGETED_CSV = TARGETED_DIR / "results/targeted_panel_delete_overexpress_merged.csv"
TARGETED_GOAL_ALT_FIGURES = TARGETED_DIR / "figures/goal_vs_alt_shift"
TARGETED_DVO_FIGURES = TARGETED_DIR / "figures/delete_vs_overexpress_shift"
ALLGENE_CSV = HERE / "tables/allgene_delete_overexpress_shift.csv"
ALLGENE_GOAL_ALT_CSV = HERE / "tables/allgene_goal_vs_alt_shift.csv"
ALLGENE_GOAL_ALT_FIGURES = HERE / "figures/goal_vs_alt_shift"
ALLGENE_DVO_FIGURES = HERE / "figures/delete_vs_overexpress_shift"
OUT = HERE / "reports/delete_vs_overexpress_shift_report.html"


def embed(path: Path) -> str:
    """Base64 data URI. A local HTML file's relative <img> links are silently
    dropped by Safari's/WebKit's file:// sandbox for anything outside the
    report's own directory tree, so this is the only cross-browser way to make
    a self-contained report that just works when double-clicked."""
    return f"data:image/png;base64,{base64.b64encode(path.read_bytes()).decode('ascii')}"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]
ARMS = ("delete", "overexpress")


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


def targeted_goal_alt_summary() -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        for comparison in COMPARISONS:
            df = pd.read_csv(TARGETED_DIR / "results" / arm / f"targeted_{arm}_{comparison}.csv")
            alt_col = next(c for c in df.columns if c.startswith("Shift_to_alt_end_"))
            sig = df[df["Sig"] == 1]
            rows.append({
                "arm": arm, "comparison": comparison, "n_genes": len(df),
                "n_significant": len(sig),
                "pearson_r_goal_alt": round(sig["Shift_to_goal_end"].corr(sig[alt_col]), 3) if len(sig) > 2 else float("nan"),
            })
    return pd.DataFrame(rows)


def allgene_goal_alt_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        for comparison in COMPARISONS:
            sub = df[(df["arm"] == arm) & (df["comparison"] == comparison)]
            sig = sub[sub["sig"] == 1]
            rows.append({
                "arm": arm, "comparison": comparison, "n_genes": len(sub),
                "n_significant": len(sig),
                "pearson_r_goal_alt": round(sig["shift_to_goal"].corr(sig["shift_to_alt"]), 3) if len(sig) > 2 else float("nan"),
            })
    return pd.DataFrame(rows)


def goal_alt_image_grid(figures_dir: Path) -> str:
    rows = []
    for comparison in COMPARISONS:
        cells = []
        for arm in ARMS:
            path = figures_dir / arm / f"{comparison}.png"
            cells.append(f'<td><img src="{embed(path)}" alt="{arm} {comparison} goal vs alt shift"></td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    header = "".join(f"<th>{arm}</th>" for arm in ARMS)
    return f'<table class="figure-grid"><thead><tr>{header}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def delete_overexpress_image_grid(figures_dir: Path, columns: int = 2) -> str:
    cells = [
        f'<td><img src="{embed(figures_dir / f"{c}.png")}" alt="{c} delete vs overexpress shift"></td>'
        for c in COMPARISONS
    ]
    rows = ["<tr>" + "".join(cells[i:i + columns]) + "</tr>" for i in range(0, len(cells), columns)]
    return f'<table class="figure-grid"><tbody>{"".join(rows)}</tbody></table>'


def main() -> None:
    targeted = pd.read_csv(TARGETED_CSV)
    allgene = pd.read_csv(ALLGENE_CSV)
    targeted_summary = summarize(targeted)
    allgene_summary = summarize(allgene)
    targeted_goal_alt = targeted_goal_alt_summary()
    allgene_goal_alt = allgene_goal_alt_summary(pd.read_csv(ALLGENE_GOAL_ALT_CSV))

    body = f"""
<h1>T-cell perturbation shift report: goal-vs-alt specificity and delete-vs-overexpress concordance</h1>
<p><strong>Generated:</strong> {datetime.now(timezone.utc).isoformat()}</p>

<p><strong>Scope: T cells only.</strong> "Normal", "SCLC", and "LUAD" below name the
tumor context a T cell was sampled from &mdash; T cells recovered from normal lung,
from SCLC tumors, or from LUAD tumors &mdash; not tumor or epithelial cells
themselves. Every shift, classifier, and gene hit in this report describes the
T-cell compartment.</p>

<p>Two datasets, same underlying perturbations, at two scales:</p>
<ul>
<li><strong>Targeted 50-gene panel</strong> &mdash; the pre-registered 21-gene
immune panel plus 29 top drivers from the prior NSCLC screen, both arms
census-complete per source state.</li>
<li><strong>Whole-genome screen</strong> &mdash; every gene present in each
held-out cell's token sequence, delete and overexpress arms, computed on
thinkstation2.</li>
</ul>

<h1>Part 1 &mdash; shift to goal vs shift to alt state</h1>
<p>Every perturbed cell has two off-target-free readouts from the <em>same</em>
edit: how far it moved toward the comparison's goal state, and how far it
moved toward the third, unnamed "alt" state. A gene whose effect is specific
to the goal sits near the x-axis (large goal shift, &asymp;0 alt shift); a gene
near the dashed y&nbsp;=&nbsp;x line is moving the cell toward both other
states about equally &mdash; not specific to the goal at all. Marker area
scales with N detections (per-gene cell count); color marks whether the goal
shift is FDR-significant. The most-displaced significant genes in each panel
are labeled.</p>

<h2>Targeted 50-gene panel</h2>
{table_html(targeted_goal_alt)}
{goal_alt_image_grid(TARGETED_GOAL_ALT_FIGURES)}

<h2>Whole-genome screen</h2>
<p>Axes in each panel are clipped to that panel's 1st&ndash;99th percentile
(heavy-tailed response); the panel title reports how many points fall outside
the drawn frame.</p>
{table_html(allgene_goal_alt)}
{goal_alt_image_grid(ALLGENE_GOAL_ALT_FIGURES)}

<h1>Part 2 &mdash; deletion vs overexpression concordance</h1>
<p>Each gene-comparison pair also has two shift measurements toward the same
goal state: one from deleting the gene, one from overexpressing it. If the
gene genuinely drives the model toward or away from that state, the two
shifts should be significant, adequately detected, and opposite in sign
&mdash; overexpression pushing one way, deletion undoing it. That is the
<code>concordant</code> criterion used throughout this repo's perturbation
work: <code>delete_fdr &lt; 0.05</code>, <code>overexpress_fdr &lt; 0.05</code>,
<code>delete_shift &times; overexpress_shift &lt; 0</code>, and
<code>delete_n &ge; 25</code>.</p>

<h2>Targeted 50-gene panel</h2>
{table_html(targeted_summary)}
{delete_overexpress_image_grid(TARGETED_DVO_FIGURES)}

<h2>Whole-genome screen</h2>
{table_html(allgene_summary)}
{delete_overexpress_image_grid(ALLGENE_DVO_FIGURES)}

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
<li>Targeted panel, goal vs alt: <code>sclc_validation/perturbation_workflow/targeted_panel/plot_goal_vs_alt_shift.py</code>
reading the committed <code>results/{{delete,overexpress}}/targeted_*.csv</code> tables.</li>
<li>Whole-genome, goal vs alt: <code>sclc_validation/primary_test_perturbation/scripts/plot_goal_vs_alt_shift.py</code>
reading <code>stats/{{delete,overexpress}}/heldout_allgene_*.csv</code> under
<code>SCLC_PERTURBATION_ROOT</code> (thinkstation2), writing
<code>tables/allgene_goal_vs_alt_shift.csv</code>.</li>
<li>Targeted panel, delete vs overexpress: <code>sclc_validation/perturbation_workflow/targeted_panel/plot_delete_vs_overexpress_shift.py</code>
reading the committed <code>results/targeted_panel_delete_overexpress_merged.csv</code>.</li>
<li>Whole-genome, delete vs overexpress: <code>sclc_validation/primary_test_perturbation/scripts/plot_delete_vs_overexpress_shift.py</code>
reading the same <code>stats/</code> tables, writing
<code>tables/allgene_delete_overexpress_shift.csv</code>.</li>
</ul>
"""
    document = (
        "<html><head><meta charset='utf-8'>"
        "<title>T-cell delete vs overexpress shift</title>"
        "<style>body{font-family:system-ui;max-width:1200px;margin:2rem auto;"
        "line-height:1.5;padding:0 1rem} table{border-collapse:collapse;margin:1rem 0}"
        " th,td{border:1px solid #ccd;padding:.35rem .6rem;text-align:right}"
        " th{background:#eef2f5} td:first-child,th:first-child{text-align:left}"
        " code{background:#f4f4f4;padding:0 .2rem} img{max-width:100%;height:auto;"
        "border:1px solid #ddd;margin:.5rem 0 1.5rem}"
        " table.figure-grid{width:100%} table.figure-grid td,table.figure-grid th{border:none;"
        "padding:.25rem;text-align:center;text-transform:capitalize}"
        " table.figure-grid img{margin:0}</style></head><body>"
        + body +
        "</body></html>"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(document)
    print(f"Wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
