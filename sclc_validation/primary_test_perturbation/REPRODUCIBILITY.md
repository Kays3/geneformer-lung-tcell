# Methods and reproducibility notes

Covers the post-perturbation analysis chain in this directory: denoising, donor-level
consistency, and the ambient-risk diagnostic. The spatial validation of the resulting
programs is documented separately in
[`../spatial_validation/METHODS_denoised_programs.md`](../spatial_validation/METHODS_denoised_programs.md).

## 1. Pipeline order

Each stage consumes the previous stage's committed table. Run them in this order; none is
resumable from a later stage.

| # | Stage | Script | Key output |
|---|---|---|---|
| 0 | All-gene perturbation (GPU, ~5 days, 2 nodes) | `run_heldout_allgene.py` (in `KD/`) | `stats/{arm}/heldout_allgene_*.csv` |
| 1 | Primary report and concordance | `scripts/build_primary_report.py` | `tables/primary_concordant_hits.csv` |
| 2 | Denoising and program assignment | `scripts/build_denoised_notebook.py` | `tables/immune_cancer_candidates.csv` |
| 3 | Donor-level consistency | `scripts/donor_consistency_allgene.py` | `tables/allgene_concordant_hits_with_donor_robustness.csv` |
| 4 | Ambient-risk diagnostic | `scripts/ambient_risk_diagnostic.py` | `tables/ambient_risk_candidates.csv` |

Stage 4 reads the *donor-annotated* candidate file produced by stage 3, so stage 3 must
complete first.

## 2. Environments

Stages 3 and 4 read the large `KD/` artifacts and run on the workstation. Stage 2 and all
figure regeneration run on the laptop.

**`thinkstation2`** — NVIDIA GB10, `Linux-6.17.0-1029-nvidia-aarch64`, glibc 2.39,
interpreter `~/workspace/geneformer-uv-starter/.venv/bin/python`:

```text
python 3.12.13   numpy 2.4.4     scipy 1.18.0    pandas 3.0.5
sklearn 1.9.0    anndata 0.13.2  scanpy 1.12.3   datasets 5.0.1
torch 2.13.0+cu130
```

**Laptop** — `macOS-26.2-arm64`:

```text
python 3.13.12   numpy 2.4.4     scipy 1.17.1    pandas 2.3.3
matplotlib 3.10.8  anndata 0.12.10  scanpy 1.12.1
```

Stages 3 and 4 are CPU-only; the GPU is used only by stage 0.

## 3. Running the analyses

Run from the repository checkout on the workstation, never from a scratch copy, so the
committed code is provably the code that ran — see [`../../WORKFLOW.md`](../../WORKFLOW.md).

```bash
./tools/sync.sh
```

```bash
ssh ts2 'cd ~/workspace/geneformer-lung-tcell && N_WORKERS=10 \
  ~/workspace/geneformer-uv-starter/.venv/bin/python \
  sclc_validation/primary_test_perturbation/scripts/donor_consistency_allgene.py'
```

```bash
ssh ts2 'cd ~/workspace/geneformer-lung-tcell && \
  ~/workspace/geneformer-uv-starter/.venv/bin/python \
  sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py'
```

Figures can be regenerated from the committed tables alone, without workstation access:

```bash
python sclc_validation/primary_test_perturbation/scripts/donor_consistency_allgene.py --plots-only
```

Runtimes: stage 3 ≈ 12 min with `N_WORKERS=10` (38,180 pickles, 1.3 GB); stage 4 ≈ 1 min.

## 4. Configuration

No machine-specific path is hardcoded in these scripts. Every input is an environment
override with a default:

| Variable | Used by | Default |
|---|---|---|
| `SCLC_PERTURBATION_ROOT` | 3, 4 | `~/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation` |
| `GENEFORMER_TOKEN_DICT` | 3 | `.../Geneformer/geneformer/token_dictionary_gc104M.pkl` |
| `HTAN_H5AD` | 4 | `~/workspace/KD/sclc_luad_normal_htan_finetune/data/htan_..._prepared.h5ad` |
| `DENOISED_CANDIDATES` | 4 | `tables/immune_cancer_candidates_with_donor_robustness.csv` |
| `PRIMARY_HITS`, `DONOR_OUT_DIR`, `AMBIENT_OUT_DIR`, `N_WORKERS` | 3, 4 | see script headers |

## 5. Determinism

None of these stages samples randomly, so no seed is needed and repeated runs on identical
inputs give byte-identical output:

- **Stage 3** is a deterministic accumulation over pickle files. Results are merged by key,
  not by arrival order, so the `ProcessPoolExecutor` worker count changes runtime only.
- **Stage 4** fits `LogisticRegression(max_iter=2000)` on standardized features with
  non-shuffled 5-fold `cross_val_score`. Anchor sets are literals in the script rather
  than data-derived, so they cannot drift as the candidate list changes.

Stage 2's balanced score is recomputed within the denoised subset, so it changes if the
curated gene sets are edited. That is intended, and the gene sets are literals in
`build_denoised_notebook.py`.

## 6. Correctness checks built into the code

Both new stages assert their central assumption instead of trusting it, and both caught
real problems during development.

**Stage 3 — cell-ordering assertion.** Raw output is stored per cell-batch, with the file
name encoding a cell index. The code asserts that index maps to length-descending order
within the shard by requiring the perturbed-gene count to equal `length - offset` for the
correspondingly ranked cell. This surfaced that the offset is **2 for deletion but 3 for
overexpression**, because overexpression skips the already-first gene as a no-op. The
offset is asserted per arm rather than tolerated with a fuzzy comparison.

**Stage 3 — reconstruction validation.** Before any consistency claim, per-cell values are
re-aggregated to gene level and compared against `Shift_to_goal_end` in the Geneformer
stats tables: **165,438 gene comparisons, max absolute difference 1.11e-16, correlation
1.000000, zero detection-count mismatches**. Written to
`tables/allgene_donor_reconstruction_validation.csv`.

**Stage 4 — anchor separation.** Cross-validated AUC on the two anchor sets is reported
before candidates are ranked, and the script warns if AUC < 0.8 so a non-discriminating
score is not used. Achieved 1.000.

**Stage 4 — identifier check.** An early run silently matched zero anchors because the
h5ad is Ensembl-indexed while the anchors are symbols. The mapping is now built from the
stats tables and the script exits with a message if it cannot be constructed.

## 7. Data not under version control

Large artifacts live outside Git under `KD/` (see
[`../../RELATIONSHIP_TO_KD.md`](../../RELATIONSHIP_TO_KD.md)) and are **not** reproducible
from this repository alone:

- `raw/{delete,overexpress}/{source}/*.pickle` — 38,180 files, 1.3 GB, the per-cell
  perturbation output that stage 3 consumes.
- `htan_sclc_luad_normal_tcells_prepared.h5ad` — 442 MB, 46,140 × 24,540.
- `tables/allgene_per_donor_shifts.parquet` — 7.8 MB stage-3 intermediate, excluded by the
  repository's `*.parquet` rule and regenerable by re-running stage 3.

Everything needed to re-derive the reported numbers *from the committed tables* — summary
statistics, classifications and figures — is tracked.

## 8. Provenance

| Commit | Stage |
|---|---|
| `bc0f30d` | Denoising analysis (stage 2) |
| `2f3b8ba`, `6bd1e8f` | Donor consistency (stage 3) and the per-arm offset fix |
| `c0cb67d` | Donor consistency results |
| `55f5776`, `d0cef97`, `12b3b03`, `5f54cf7` | Ambient diagnostic (stage 4) and three fixes |

## 9. Known limitations carried by the report

Stated here so they are not lost between documents:

- **Donor consistency rests on 4 luad and 3 sclc donors**, and the single normal donor
  makes all 408 normal-source hits unassessable. "All donors agree" is a weaker claim than
  it sounds. The sclc base is 75% one sample whose label, `PleuralEffusion`, is a sample
  type rather than a patient identifier.
- **The ambient stage is a diagnostic, not a correction.** CellBender could not be run:
  the available matrix has no empty droplets (minimum 113 UMI/cell) and is already subset
  to T cells, leaving the ambient profile unidentifiable. A true run needs primary
  CellRanger output from the HTAN portal, which is a data-access task and remains open.
- **Nothing here makes the in-silico perturbations causal.** Deletion and overexpression
  are model-input interventions on a rank-value encoding, not experimental perturbations.
