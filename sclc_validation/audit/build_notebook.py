#!/usr/bin/env python3
"""Build the executable SCLC data-feasibility audit notebook."""

from pathlib import Path

import nbformat as nbf


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "notebooks" / "SCLC_data_feasibility_audit.ipynb"
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
        """# SCLC data feasibility audit

## TL;DR

**GO, with conditions.** Use the HTAN/CELLxGENE T-cell object for the primary
SCLC–LUAD–normal comparison, GSE263196 for orthogonal spatial validation, and
OMIX002441 for cross-platform sensitivity. GSE261348 remains blocked for
treatment-response validation until its clinical outcomes are obtained and
joined."""
    ),
    nbf.v4.new_markdown_cell(
        """## Context and methods

The conference abstract claims SCLC-specific T-cell dysfunction relative to
LUAD and normal lung. This notebook reads the machine-generated audit outputs,
checks the main acceptance assertions, and visualizes cohort scale at the donor
level. The full source-integrity checks are implemented in
`../audit_sclc_data.py`."""
    ),
    nbf.v4.new_code_cell(
        """from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt

AUDIT_DIR = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd() / "sclc_validation/audit"
RESULTS = AUDIT_DIR / "results"
summary = json.loads((RESULTS / "audit_summary.json").read_text())
htan = pd.read_csv(RESULTS / "htan_tcell_summary_by_disease.csv")
spatial = pd.read_csv(RESULTS / "gse263196_spatial_file_audit.csv")
omix = pd.read_csv(RESULTS / "omix_tcell_counts_by_patient.csv")
htan"""
    ),
    nbf.v4.new_markdown_cell(
        """## Data

Four open-access resources were audited: one integrated T-cell object spanning
SCLC/LUAD/normal, one independent SCLC single-cell cohort, one Visium SCLC
cohort, and one GeoMx pretreatment cohort."""
    ),
    nbf.v4.new_code_cell(
        """assert summary["htan_tcell_h5ad"]["raw_available"]
assert summary["htan_tcell_h5ad"]["raw_integer_sample"]
assert all(summary["htan_tcell_h5ad"]["signature_genes_in_raw_symbols"].values())
assert summary["gse263196"]["all_matrix_integrity_checks_pass"]
assert summary["gse263196"]["all_signature_genes_present"]
assert summary["gse261348"]["outcome_columns_in_workbook"] == []
print("Primary technical acceptance assertions passed.")"""
    ),
    nbf.v4.new_markdown_cell("## Results"),
    nbf.v4.new_code_cell(
        """plot = htan.sort_values("cells")
fig, ax = plt.subplots(figsize=(8, 4.8))
bars = ax.barh(plot["disease"], plot["cells"], color=["#A8DADC", "#457B9D", "#E76F51"])
ax.bar_label(bars, labels=[f"{v:,}" for v in plot["cells"]], padding=4)
ax.set_xlabel("T cells")
ax.set_title("Audited HTAN T-cell cohort by disease")
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.show()"""
    ),
    nbf.v4.new_code_cell(
        """pd.DataFrame([
    {
        "dataset": "HTAN T cells",
        "scale": f'{summary["htan_tcell_h5ad"]["cells"]:,} cells',
        "key check": "raw counts and donor labels present",
        "decision": "primary discovery",
    },
    {
        "dataset": "GSE263196",
        "scale": f'{summary["gse263196"]["total_in_tissue_spots"]:,} spots / 5 samples',
        "key check": "complete Visium bundles",
        "decision": "spatial validation",
    },
    {
        "dataset": "OMIX002441",
        "scale": f'{summary["omix002441"]["t_cells"]:,} T cells / 11 patients',
        "key check": "only 2 robust paired donors",
        "decision": "sensitivity analysis",
    },
    {
        "dataset": "GSE261348",
        "scale": f'{summary["gse261348"]["expression_aois"]} AOIs',
        "key check": "clinical outcomes absent",
        "decision": "blocked pending join",
    },
])"""
    ),
    nbf.v4.new_markdown_cell(
        """## Takeaways

1. The current LUAD/LUSC repository work is exploratory provenance, not SCLC
   validation.
2. The HTAN object is sufficiently large for donor-held-out SCLC-versus-LUAD
   modeling, although the four-donor normal group limits generalization.
3. Visium spots must be deconvolved or annotated and aggregated to the five
   patient-level replicates; they must not be passed to Geneformer as if they
   were T cells.
4. OMIX002441 is valuable for direction/rank concordance, not as the main
   training cohort.
5. Treatment-response claims are not testable from the downloaded GSE261348
   workbooks without an external clinical outcome join."""
    ),
]

nbf.write(nb, OUTPUT)
print(OUTPUT)
