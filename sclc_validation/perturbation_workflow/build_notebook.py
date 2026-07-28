#!/usr/bin/env python3
"""Build the executable SCLC/LUAD/normal perturbation workflow notebook."""

from pathlib import Path

import nbformat as nbf


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "notebooks" / "SCLC_LUAD_normal_perturbation_workflow.ipynb"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

nb = nbf.v4.new_notebook()
nb["metadata"]["kernelspec"] = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3",
}
nb["metadata"]["language_info"] = {"name": "python", "version": "3"}
nb["cells"] = [
    nbf.v4.new_markdown_cell(
        """# SCLC / LUAD / normal T-cell fine-tuning and perturbation

## TL;DR

Fine-tuning and held-out classification are **complete**: 0.919 accuracy,
0.903 macro F1 on the donor-held-out HTAN test set, notably higher than the
prior LUAD/LUSC/normal model (0.783 / 0.758). The donor-aware, all-gene
deletion + overexpression in-silico perturbation screen is **designed and
launched, currently paused** pending a runtime/scope decision (the full
752-shard sweep -- both perturbation types across all three sources --
projects to roughly 7-8 days on the available GPU; batch-size tuning only
bought a modest gain since the workload is compute-bound, not
batch-limited)."""
    ),
    nbf.v4.new_markdown_cell(
        """## Context and methods

This is the SCLC-inclusive counterpart to the repository's LUAD/LUSC/normal
T-cell classifier and all-gene perturbation screen, built on the HTAN
T-cell cohort identified by `sclc_validation/audit/`. Full methodology,
including the donor-disjoint-across-disease guard for three patients who
each contribute both LUAD tumor and normal-tissue T cells, is in
`../METHODS.md`. This notebook reads the machine-generated result tables in
`../results/`."""
    ),
    nbf.v4.new_code_cell(
        """from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt

WORKFLOW_DIR = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd() / "sclc_validation/perturbation_workflow"
RESULTS = WORKFLOW_DIR / "results"

test_metrics = json.loads((RESULTS / "test_metrics.json").read_text())
classification_report = json.loads((RESULTS / "classification_metrics_report.json").read_text())
confusion = pd.read_csv(RESULTS / "test_confusion_matrix.csv", index_col=0)
donor_split = pd.read_csv(RESULTS / "htan_donor_split_assignment.csv")
cohort_summary = pd.read_csv(RESULTS / "htan_cohort_characterization_summary.csv")
cohort_summary"""
    ),
    nbf.v4.new_markdown_cell(
        """## Cohort and split

HTAN/CELLxGENE T-cell object only (46,140 T cells, 42 donors) -- no other
dataset is blended in. The normal class is thin by design (4 donors total,
2 train / 1 eval / 1 test); this is a named limitation, not silently treated
as equivalent in power to the SCLC/LUAD arms."""
    ),
    nbf.v4.new_code_cell(
        """assert test_metrics["donor_leakage_check"] == "PASS"
cross_disease_donors = (
    donor_split.groupby("individual")["disease"].nunique()
)
assert (cross_disease_donors > 1).sum() == 3, "Expected 3 cross-disease (paired tumor/normal) donors"
print("Donor-leakage check (including cross-disease guard): PASS")
print("Cross-disease donors:", cross_disease_donors[cross_disease_donors > 1].index.tolist())"""
    ),
    nbf.v4.new_markdown_cell("## Classification results"),
    nbf.v4.new_code_cell(
        """print(f'Held-out accuracy: {classification_report["overall_accuracy"]:.3f}')
print(f'Held-out macro F1:  {classification_report["overall_macro_f1"]:.3f}')
print()
print("Prior LUAD/LUSC/normal model (calibration reference, different classes/cohort):")
prior = classification_report["prior_workflow_context"]
print(f'  accuracy {prior["held_out_accuracy"]:.3f}, macro F1 {prior["held_out_macro_f1"]:.3f}')
pd.DataFrame(classification_report["per_class"]).sort_values("disease")"""
    ),
    nbf.v4.new_code_cell(
        """fig, ax = plt.subplots(figsize=(5.5, 4.8))
im = ax.imshow(confusion.values, cmap="Blues")
ax.set_xticks(range(len(confusion.columns)), confusion.columns, rotation=30, ha="right")
ax.set_yticks(range(len(confusion.index)), confusion.index)
for i in range(confusion.shape[0]):
    for j in range(confusion.shape[1]):
        ax.text(j, i, f"{confusion.values[i, j]:,}", ha="center", va="center",
                color="white" if confusion.values[i, j] > confusion.values.max() / 2 else "black")
ax.set_xlabel("Predicted")
ax.set_ylabel("True")
ax.set_title("Held-out confusion matrix")
plt.tight_layout()
plt.show()"""
    ),
    nbf.v4.new_markdown_cell(
        """SCLC's main confusion is with LUAD (~20% of SCLC test cells called LUAD,
almost none with normal) -- the biologically closer pair of the two tumor
types. Normal's near-perfect recall is a single-donor result and should be
read as one data point, not a population estimate."""
    ),
    nbf.v4.new_markdown_cell(
        """## In-silico perturbation -- design and status

Every gene token in every held-out test cell is perturbed individually,
once per type (**delete**: removed from the rank-value encoding;
**overexpress**: moved to the front, simulating maximal relative
expression), and scored against training-cells-only disease reference
embeddings. Three source screens (SCLC, LUAD, normal) recover all six
directional comparisons per type. A gene whose deletion and overexpression
move a cell in opposite directions is the strongest evidence tier.

An earlier plan to restrict this screen to the top 3000 variable genes was
reverted: Geneformer's `InSilicoPerturber` has no built-in "sweep
individually, restricted to a list" mode -- passing a gene list perturbs the
whole list as one combined event, not one gene at a time. Running all genes
avoids a reduced-context-embedding deviation from the prior workflow's
methodology, at the cost of a substantially larger compute budget.

**Status at notebook build time**: the sweep is designed, launched, and
currently paused after a throughput-tuning experiment (forward_batch_size
128 vs. 16 only reduced per-shard time from ~15 min to ~12.9 min -- the
workload is GPU-compute-bound, not batch-limited). The full 752-shard sweep
(2 perturbation types x 3 sources) projects to roughly 7-8 days; results will
be added to `../RESULTS.md` once a runtime/scope decision is made and the
sweep completes."""
    ),
    nbf.v4.new_markdown_cell(
        """## Takeaways so far

1. The SCLC/LUAD/normal classifier clears a materially higher held-out bar
   (0.919 accuracy / 0.903 macro F1) than the prior LUAD/LUSC/normal model,
   though the two are not a head-to-head comparison (different classes,
   different cohort).
2. SCLC-vs-LUAD is the classifier's main source of confusion; SCLC-vs-normal
   is essentially clean.
3. The normal class's single test donor means its metrics should be treated
   as a single-patient data point, not validated population performance.
4. The perturbation screen's full scope (all genes, both delete and
   overexpress) is real compute, not a quick follow-up -- a scope or
   runtime decision is needed before it can complete."""
    ),
]

nbf.write(nb, OUTPUT)
print(OUTPUT)
