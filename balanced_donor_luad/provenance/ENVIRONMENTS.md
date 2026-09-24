# Compute environments (deletion on thinkstation1, overexpression on thinkstation2)

Recorded 2026-09-24. Both hosts must stay identical in everything listed
here. The cross-host equivalence gate (Amendment 3, 3.5) tests the parts
that cannot be listed, such as driver and numerical effects.

| Item | thinkstation1 | thinkstation2 |
|---|---|---|
| GPU | NVIDIA GB10 | NVIDIA GB10 |
| Driver | **595.84** | **580.173.02** (differs, see below) |
| Python | 3.12.13 (uv) | 3.12.13 (uv) |
| uv project | `geneformer-uv-starter/sclc_analysis` | `balanced_donor_luad/env/sclc_analysis` (copy) |
| `uv.lock` sha256 prefix | 781e42036c5f8873 | 781e42036c5f8873 (identical; `uv sync --frozen`) |
| `pyproject.toml` sha256 prefix | 5e3a9bd9a866822e | 5e3a9bd9a866822e |
| Sorted `uv pip freeze` (excl. editable geneformer), sha256 prefix | 353d6596cbcab16c | 353d6596cbcab16c |
| Key pins | torch 2.13.0+cu130, transformers 4.46.0, datasets 5.0.1, pandas 3.0.5, numpy 2.5.2 | same |
| Geneformer | `geneformer-uv-starter/Geneformer` @ f45a6c7de57f + bf16 export patch | fresh clone of huggingface.co/ctheodoris/Geneformer @ f45a6c7de57f + the same patch (`apply_bf16_patch.sh`) |
| Patch check | reverse dry-run clean | reverse dry-run clean |
| 316M weights sha256 | 965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3 | identical |
| Tokenised cohort (3 files) | provenance/tokenized_dataset_sha256.txt | identical, all 3 files |

- **Driver difference.** The drivers differ (595.84 vs 580.173.02). Both meet
  CUDA 13.0's minimum. Any numerical effect is exactly what the cross-host
  equivalence gate measures before Phase 6. If the gate fails, Phase 6 does
  not start.
- **Not used:** `/srv/lab/geneformer`, where the 316M file is an LFS pointer,
  and `tools/lab_env.sh` / `paths.env`.
