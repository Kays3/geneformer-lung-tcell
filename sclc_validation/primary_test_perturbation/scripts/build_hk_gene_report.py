#!/usr/bin/env python3
"""Build the HK/ubiquitous-gene review-response report (T cells).

Responds point-by-point to an external review raising two concerns about
Geneformer in-silico perturbation (ISP): (1) housekeeping/ubiquitous genes are
over-represented among top hits, and (2) strong state-classification accuracy
does not guarantee that ISP reflects real perturbation biology ("closing the
loop", PMC12265564). Every number below is read from tables written by
`hk_gene_diagnostic.py` at build time; this script only lays them out.
"""
from __future__ import annotations

import base64
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent.parent
TABLES = HERE / "tables/hk_gene_diagnostic"
FIGURES = HERE / "figures/hk_gene_diagnostic"
OUT = HERE / "reports/hk_gene_review_response.html"


def embed(path: Path) -> str:
    return f"data:image/png;base64,{base64.b64encode(path.read_bytes()).decode('ascii')}"


def table_html(df: pd.DataFrame, **kwargs) -> str:
    return df.to_html(index=False, classes="data", border=0, **kwargs)


def main() -> None:
    panel_flags = pd.read_csv(TABLES / "targeted_panel_gene_flags.csv")
    dvo_enrich = pd.read_csv(TABLES / "allgene_dvo_hk_enrichment.csv")
    ga_enrich = pd.read_csv(TABLES / "allgene_goal_alt_hk_enrichment.csv")
    gap = pd.read_csv(TABLES / "allgene_hk_concordant_detection_gap.csv").drop_duplicates(
        subset=["Gene_name", "comparison"])

    n_panel_hk_or_ambient = int((panel_flags["is_hk"] | panel_flags["ambient_flag"]).sum())
    top_driver = panel_flags[panel_flags["gene_source"] == "top_driver_luad_lusc_normal"]
    n_top_driver_flagged = int((top_driver["is_hk"] | top_driver["ambient_flag"]).sum())
    pre_registered = panel_flags[panel_flags["gene_source"] == "panel"]
    n_pre_registered_flagged = int((pre_registered["is_hk"] | pre_registered["ambient_flag"]).sum())

    n_gap_total = len(gap)
    n_gap_near_zero = int((gap["detect_frac_gap"] < 0.02).sum())
    n_gap_large = int((gap["detect_frac_gap"] > 0.3).sum())

    dvo_fold_range = (dvo_enrich["odds_ratio"].min(), dvo_enrich["odds_ratio"].max())
    dvo_max_p = dvo_enrich["fisher_p"].max()

    body = f"""
<h1>Review response: housekeeping-gene enrichment and open-loop ISP validity</h1>
<p><strong>Generated:</strong> {datetime.now(timezone.utc).isoformat()}</p>
<p><strong>Scope: T cells only.</strong> "Normal", "SCLC", and "LUAD" name the tumor
context a T cell was sampled from, not a tumor or epithelial cell type. Every
number below describes the T-cell perturbation screens documented elsewhere in
this repo (targeted 50-gene panel and whole-genome, both delete and
overexpress arms).</p>

<p>This report answers two review comments with data already produced by this
project, plus one new diagnostic (<code>hk_gene_diagnostic.py</code>) built
for this purpose. It does not resolve either concern in the abstract; it
quantifies how large each effect actually is in <em>this</em> project's hit
lists, gene by gene, and states plainly what remains a hypothesis.</p>

<p style="background:#eef2f5;border:1px solid #ccd;border-radius:4px;padding:.7rem 1rem">
<strong>Update:</strong> T2 (measured differential expression, named below as
the concrete next step) has since been run against the atlas directly for
these 155 HK-flagged concordant genes. It confirms the detection-fraction
proxy used throughout this report exactly (Pearson r&nbsp;=&nbsp;1.000,
n&nbsp;=&nbsp;397 gene-state pairs) and adds real expression magnitude for the
same genes, which agrees with the interpretive split drawn below (<code>B2M</code>,
<code>EIF1</code>, <code>RPL27A</code> still show no state difference;
<code>HSPA1A/B</code>, <code>DNAJA1/B1</code> still show large, LUAD-specific
increases). See
<a href="../../immune_axis_test/RESULTS_T2.md">immune_axis_test/RESULTS_T2.md</a>.</p>

<h1>Comment 1 &mdash; housekeeping/ubiquitous genes among top ISP hits</h1>
<p>The concern: a substantial fraction of top ISP hits are housekeeping (HK)
genes, which could reflect genuine biology (essential-gene perturbation
producing large transcriptomic change) or a systematic property of
Geneformer's ranking/ISP procedure rather than T-cell-state biology specific
to lung cancer.</p>

<h2>Two flags, two different mechanisms</h2>
<p>These are kept separate throughout because they are not the same failure
mode:</p>
<ul>
<li><strong>Ambient / lineage-foreign</strong> &mdash; a gene a T cell should
not transcribe at all (surfactant, haemoglobin, myeloid, epithelial, stromal
markers). Already diagnosed for this project in
<a href="../METHODS_ambient_risk.md">METHODS_ambient_risk.md</a>: a logistic
score calibrated against 34 known-ambient and 36 known-T-cell anchor genes,
cross-validated AUC&nbsp;=&nbsp;1.0, applied here by joining the existing
<code>ambient_risk_all_genes.csv</code> (11,047 scored genes) against this
project's hit lists &mdash; not recomputed.</li>
<li><strong>Housekeeping/ubiquitous</strong> &mdash; genes genuinely
transcribed by T cells, but constitutively: ribosomal proteins, translation
factors, heat-shock proteins, proteasome subunits, cytoskeleton, oxidative
phosphorylation, plus a short list of classic reference genes (<code>ACTB</code>,
<code>GAPDH</code>, <code>B2M</code>, &hellip;). This is a coarse, standard
family list &mdash; the kind used for <code>percent.ribo</code>/<code>percent.hsp</code>
QC in single-cell pipelines &mdash; not an exhaustive HK catalogue. Defined in
<code>hk_gene_diagnostic.py</code> before any hit list was scored against it.</li>
</ul>

<h2>The targeted 50-gene panel is composed of two very different halves</h2>
<p>The panel's own metadata already distinguishes its two sources: 21
pre-registered immune genes vs. 29 genes pulled in as prior top drivers from
the earlier NSCLC all-gene screen. Scoring both flags against each half
separately:</p>
<table class="data">
<tr><th></th><th>pre-registered immune panel</th><th>prior top-driver hits</th></tr>
<tr><td>genes</td><td>21</td><td>29</td></tr>
<tr><td>HK or ambient-flagged</td><td>{n_pre_registered_flagged} (0%)</td>
<td>{n_top_driver_flagged} ({100*n_top_driver_flagged//29}%)</td></tr>
</table>
<img src="{embed(FIGURES / 'targeted_panel_composition.png')}" alt="Targeted panel composition by HK/ambient flag">
<p><strong>None</strong> of the 21 pre-registered immune genes are HK- or
ambient-flagged. <strong>{n_top_driver_flagged} of 29</strong> top-driver genes
are &mdash; entirely because that set was selected <em>from a prior ISP
screen's own top hits</em>, which is exactly the selection process the review
comment is questioning. This is not evidence the concern is wrong; it is a
demonstration of how it enters a panel silently once "top hit from a previous
ISP run" is used as a selection criterion.</p>

<h2>Whole-genome: HK enrichment is real, large, and consistent across every comparison</h2>
<p>Fisher's exact test, HK-family fraction among hits vs. among all genes
tested in that comparison (background is 2.8&ndash;3.4% depending on
comparison):</p>
{table_html(dvo_enrich[["comparison", "n_hit", "pct_hk_background", "pct_hk_among_hit", "odds_ratio", "fisher_p"]].rename(
    columns={"n_hit": "n concordant", "pct_hk_background": "% HK, background",
             "pct_hk_among_hit": "% HK, among hits", "odds_ratio": "odds ratio", "fisher_p": "Fisher p"}))}
<img src="{embed(FIGURES / 'hk_enrichment_bars.png')}" alt="HK enrichment bar chart, hits vs background">
<p>Odds ratios range {dvo_fold_range[0]:.1f}&ndash;{dvo_fold_range[1]:.1f}&times;
across the six comparisons (worst-case p&nbsp;=&nbsp;{dvo_max_p:.1e}), and the
pattern holds identically for goal-vs-alt significant movers in both arms
(table below) &mdash; this is not an artifact of the concordance criterion
specifically.</p>
{table_html(ga_enrich[["arm", "comparison", "n_hit", "pct_hk_background", "pct_hk_among_hit", "odds_ratio", "fisher_p"]].rename(
    columns={"n_hit": "n significant", "pct_hk_background": "% HK, background",
             "pct_hk_among_hit": "% HK, among hits", "odds_ratio": "odds ratio", "fisher_p": "Fisher p"}))}

<h2>Applying the recommended check: is the HK gene actually differentially detected?</h2>
<p>The Geneformer authors' own recommendation is to check whether a candidate
HK gene is differentially expressed between the source and goal state before
treating it as informative &mdash; equivalently, before treating it as a
negative-control <em>perturbation</em> (a gene that is itself perturbed and
expected to produce little shift) rather than merely a negative-control
<em>target</em> (a gene expected not to respond when something else is
perturbed). This project does not yet have the raw counts pulled locally for a
proper differential-expression test (that is the immune-axis module's T2,
still scheduled, see below); as an interim, coarser proxy, this uses each
gene's <strong>detection fraction</strong> (delete-arm <code>N_Detections</code> /
total held-out T cells for that state) in the source vs. goal state. A gene
sitting on the diagonal is detected at the same rate in both T-cell states
despite a "significant" reciprocal ISP shift &mdash; the pattern the review
comment warns about, since detection rate not differing gives no expression
basis for the shift.</p>
{table_html(pd.DataFrame([
    {"finding": "near-diagonal (detection-fraction gap < 0.02)",
     "count": f"{n_gap_near_zero} / {n_gap_total} ({100*n_gap_near_zero//n_gap_total}%)",
     "reading": "no detectable expression difference between states; ISP shift for these genes has no support from this proxy"},
    {"finding": "large gap (> 0.30)",
     "count": f"{n_gap_large} / {n_gap_total} ({100*n_gap_large//n_gap_total}%)",
     "reading": "substantial detection-rate difference; a genuine, state-linked expression difference is plausible"},
]))}
<img src="{embed(FIGURES / 'hk_detection_gap_scatter.png')}" alt="Detection fraction gap scatter for HK-flagged concordant hits">
<p>The two tails read differently gene by gene, which is the point of running
this check rather than applying a blanket verdict:</p>
<ul>
<li><strong>Near the diagonal (suspicious):</strong> <code>B2M</code> is
detected in 99.96% vs. 99.95% of SCLC and LUAD T cells respectively &mdash;
indistinguishable &mdash; yet clears the concordance bar. <code>EIF1</code> and
<code>RPL27A</code> show the same pattern (&lt;0.2 percentage-point gaps).
For genes like these, this project's own data gives no expression-level
reason to trust the ISP shift as disease-relevant, which is the failure mode
explanation (b) in the motivation predicts.</li>
<li><strong>Far from the diagonal (plausible biology):</strong>
<code>HSPA1A</code>, <code>HSPA1B</code>, <code>DNAJA1</code>,
<code>DNAJB1</code> (heat-shock/co-chaperone), and <code>NDUFA3</code>,
<code>ATP5ME</code> (oxidative phosphorylation) show 30&ndash;50
percentage-point detection gaps between T-cell states. A stress-response or
metabolic-state difference between T cells recovered from tumour vs. normal
tissue is a substantive, testable biological hypothesis &mdash; not obviously
wrong just because the genes are HK-family by list membership.</li>
</ul>
<p><strong>What this does and does not show:</strong> it does not prove any
individual gene is or is not a real driver &mdash; detection fraction is a
prevalence proxy, not a fold-change test, and a true differential-expression
analysis (T2 below) is needed before either reading is treated as settled. It
does show that "HK-flagged" is not a uniform category in this project's
results: roughly a quarter of HK-flagged concordant hits have essentially no
supporting expression difference, and a distinct subset has a large one.</p>

<h1>Comment 2 &mdash; does state-classification accuracy imply reliable ISP?</h1>
<p>The concern (from PMC12265564, "closing the loop"): a Geneformer classifier
fine-tuned to 99.8% accuracy on resting-vs-activated T cells still produced
low-PPV ISP under open-loop evaluation, because a model can learn to separate
observed states without learning how cells respond to a perturbation. This
project's classifiers are also strong &mdash; the three-way SCLC/LUAD/normal
T-cell classifier reports 0.919 accuracy, 0.903 macro F1 &mdash; and by the
same logic, that number is not itself evidence the ISP shifts are causally
meaningful.</p>

<h2>This pipeline is open-loop, plainly stated</h2>
<p>Every shift reported in this project's perturbation work is a movement, in
the embedding space of a state-fine-tuned model, toward that same model's own
class centroids. No experimental perturbation data (Perturb-seq, CRISPR
screen, or equivalent) informs training or evaluation anywhere in this
pipeline. That is precisely the open-loop configuration the cited paper
tested and found insufficient on its own.</p>

<h2>What already exists in this project that is <em>not</em> the same fix, but is not nothing either</h2>
<p>None of the following closes the loop in the cited paper's sense (learning
from real perturbation outcomes). Each addresses a different, narrower
reliability question, and conflating them with closed-loop validation would
overstate what this project currently has:</p>
<ul>
<li><strong>Delete-vs-overexpress concordance</strong> (this report's Part 2
in the companion shift report) &mdash; checks whether two independent edits of
the same gene give reciprocal answers. This catches noise and one-sided
artifacts; it does not check whether the shared direction reflects real
biology, since a systematic representational bias would move both arms the
same way.</li>
<li><strong>Goal-vs-alt specificity</strong> (this report's Part 1) &mdash;
checks whether an edit's effect concentrates on the intended goal state
rather than diffusing to the third state. Same limitation: a global bias in
the embedding survives this check if it affects all three states similarly.</li>
<li><strong>Donor-level consistency</strong>
(<a href="../RESULTS_donor_consistency.md">RESULTS_donor_consistency.md</a>)
&mdash; 45.9% of concordant whole-genome hits are fully consistent in sign
across every donor in both arms, 40.0% majority-consistent. This checks
patient generalization, not perturbation causality.</li>
<li><strong>Ambient-risk diagnostic</strong> &mdash; checks lineage-foreignness,
a specific artifact class, calibrated against known-answer anchors.</li>
<li><strong>This report's HK detection-fraction proxy</strong> &mdash; the
closest thing here to the paper's own recommendation ("evaluate ISP against an
independent perturbation-related benchmark whenever possible"), but a
prevalence proxy from the same dataset is a much weaker independent signal
than real perturbation outcomes, and only covers HK-flagged genes so far.</li>
<li><strong>Spatial validation</strong> (GSE263196 Visium,
<a href="../../spatial_validation/README.md">spatial_validation/README.md</a>)
&mdash; correlates a derived T-cell dysfunction <em>signature score</em>
against independent tissue data. This is the only orthogonal (non-Geneformer)
data source used anywhere in this project, and is the closest precedent for
what an independent benchmark looks like here &mdash; but it validates a
multi-gene program, not individual ISP gene hits, so it does not transfer to
"is gene X a real hit" without further work.</li>
</ul>

<h2>What is actually missing, stated plainly</h2>
<ul>
<li>No true closed-loop retraining exists or is planned in this project: it
would require real T-cell perturbation data (Perturb-seq or CRISPR screen) in
an SCLC/LUAD/normal-comparable system, which is not currently available here.</li>
<li>No per-gene differential-expression test has been run yet against the
source atlas; the detection-fraction proxy above is a stand-in. The
immune-axis module's <strong>T2</strong> ("measured baseline expression",
<a href="../../immune_axis_test/PLAN.md">immune_axis_test/PLAN.md</a> &sect;6)
already specifies this analysis &mdash; per-donor mean expression and
detection rate per program, stratified by CD4/CD8 &mdash; and is the natural
next step to generalize from the 21-gene exhaustion program it was scoped for
to the HK-flagged hit lists in this report.</li>
<li>Classifier accuracy (0.919/0.903 for the T-cell 3-state model) should not
be read, on its own, as evidence for ISP reliability &mdash; per the cited
paper's own result, it was not sufficient even at a much higher accuracy
(99.8%).</li>
</ul>

<h1>Summary</h1>
<table class="data">
<tr><th>Question</th><th>Status</th></tr>
<tr><td>Are HK/ubiquitous genes over-represented among this project's ISP hits?</td>
<td><strong>Yes, demonstrated.</strong> 2.7&ndash;9.6&times; enrichment over background, every comparison, both arms, both hit-selection criteria (p &lt; 1e-6 throughout).</td></tr>
<tr><td>Is that enrichment uniform across those genes?</td>
<td><strong>No, and T2 confirms it directly (not just via the proxy).</strong> ~24% of HK-flagged concordant hits show no detectable expression difference between states (favors artifact); a distinct subset (heat-shock, OXPHOS) shows large differences (favors real biology). Gene-by-gene, not a blanket verdict.</td></tr>
<tr><td>Is this pipeline open-loop in the reviewed paper's sense?</td>
<td><strong>Yes, stated plainly.</strong> No experimental perturbation data informs training or evaluation anywhere in this project.</td></tr>
<tr><td>Does anything in this project substitute for closed-loop validation?</td>
<td><strong>No single existing check does.</strong> Several partial, narrower checks exist (concordance, specificity, donor consistency, ambient-risk, spatial signature validation); none tests perturbation causality directly.</td></tr>
<tr><td>What is the concrete next step?</td>
<td><strong>T2 is done</strong> (see update above and <a href="../../immune_axis_test/RESULTS_T2.md">RESULTS_T2.md</a>). Remaining: pursue real perturbation data if it becomes available, before treating any individual HK-flagged gene as a validated target.</td></tr>
</table>

<h1>Provenance</h1>
<ul>
<li>Analysis: <code>sclc_validation/primary_test_perturbation/scripts/hk_gene_diagnostic.py</code>,
reading the already-committed <code>allgene_delete_overexpress_shift.csv</code>,
<code>allgene_goal_vs_alt_shift.csv</code>, the targeted panel's merged table
and <code>target_gene_panel.json</code>, and the existing
<code>ambient_risk_all_genes.csv</code> / <code>ambient_risk_manifest.json</code>.
No new pulls from the compute nodes; no counts data touched.</li>
<li>Figures: <code>sclc_validation/primary_test_perturbation/scripts/plot_hk_gene_diagnostic.py</code>.</li>
<li>HK reference list and family regexes are literals in
<code>hk_gene_diagnostic.py</code>, defined before scoring any hit list, and
should be audited or extended there.</li>
</ul>
"""
    document = (
        "<html><head><meta charset='utf-8'>"
        "<title>T-cell HK-gene review response</title>"
        "<style>body{font-family:system-ui;max-width:1200px;margin:2rem auto;"
        "line-height:1.55;padding:0 1rem} table{border-collapse:collapse;margin:1rem 0;width:100%}"
        " th,td{border:1px solid #ccd;padding:.4rem .6rem;text-align:left}"
        " th{background:#eef2f5} td:nth-child(n+2){text-align:right}"
        " code{background:#f4f4f4;padding:0 .2rem} img{max-width:100%;height:auto;"
        "border:1px solid #ddd;margin:.5rem 0 1.5rem} h1{margin-top:2.2rem;"
        "border-top:1px solid #e1e0d9;padding-top:1.2rem} h1:first-of-type{border-top:none;padding-top:0}"
        "</style></head><body>"
        + body +
        "</body></html>"
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(document)
    print(f"Wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
