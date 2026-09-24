# Balanced-donor LUAD T-cell workflow

Geneformer classifier and in-silico perturbation (ISP) on a cohort that is
**balanced and paired by construction**:
- LUAD CD4/CD8 T cells, tumour (`tumor_primary`) vs adjacent normal
  (`normal_adjacent`);
- 43 donors, each contributing both tissues;
- exactly 100 analysis cells per donor per tissue.

## Why this exists

Earlier lung analyses had to correct for donor imbalance after the fact. For
example, one SCLC site carried 74.9% of cells. There is also a further problem:
in the source atlas, the `disease == 'normal'` donors come only from studies that
contribute no tumour T cells, so tumour-vs-normal is confounded with study.

Pairing tumour and adjacent-normal T cells **within each donor** removes study,
chemistry, site and between-donor variation from the contrast. Capping every
donor at the same cell count makes donor weighting and cell weighting identical.

## Status

| Phase | What | State |
|---|---|---|
| 0 | Selection rule, registered before screening | `provenance/PHASE0_SELECTION_RULE.md` (append-only; original sha256 `58e9b153…`) |
| 1 | Candidate survey (structural fields only) | `provenance/PHASE1_SURVEY.md` |
| 2 | Download and integrity check | `provenance/MANIFEST.json`: S3 multipart ETag recomputed and matched, sha256 `f0f7f434…` |
| 3 | Cohort draw, tokenisation, QC | `cohort/`, `qc/QC_COHORT.md`, `qc/QC_TOKENIZATION.md` |
| 4 | Classifier fine-tune and held-out evaluation | not started. Needs human GPU approval. |
| 5 | ISP pre-registration | not started |
| 6-7 | ISP run, analysis | not started |

## Source

- **Dataset:** LuCA extended atlas, CELLxGENE dataset version
  `33165751-ae33-4a65-94ee-52d9bc38f97e`, 17,616,543,376 bytes.
- **Where to fetch it:** the URL, S3 version ID, ETag and our sha256 are in
  `provenance/MANIFEST.json`.
- **Why it isn't stored here:** it is re-fetchable and verifiable, so it is
  not committed and not copied to the NAS.

## What is committed and what is regenerated

| Artifact | Where | Committed? |
|---|---|---|
| Cell IDs, draw ranks, analysis flag | `cohort/cohort_cells.csv` | yes |
| Seed, filters, hashes | `cohort/COHORT_DEFINITION.json` | yes |
| Slim h5ad (25,700 cells, raw counts) | `data/slim/`, outside the repo | no. Regenerate with `build_cohort.py`; sha256 is in the cohort definition. |
| Tokenised dataset | `data/tokenized/`, outside the repo | no. Regenerate with `tokenize_cohort.py`; file hashes are in `provenance/tokenized_dataset_sha256.txt`. |

## Reproduce

```bash
# 1. fetch the file (see provenance/MANIFEST.json for the URL), then
python scripts/verify_download.py <raw_dir> provenance/phase1_paired_donors.csv
# 2. cohort + slim h5ad (checks the file's sha256 against MANIFEST.json first)
python scripts/build_cohort.py --h5ad <raw_dir>/33165751-...h5ad \
  --manifest provenance/MANIFEST.json --paired-donors provenance/phase1_paired_donors.csv \
  --out-dir . --slim-dir <data>/slim
python scripts/qc_cohort.py --cohort cohort/cohort_cells.csv --out qc/QC_COHORT.md
# 3. tokenise (Geneformer V2; loom input path, see the script docstring)
PYTHONPATH=/srv/lab/geneformer python scripts/tokenize_cohort.py \
  --slim-dir <data>/slim --out-dir <data>/tokenized --repo-dir .
# tests: each QC check is shown to FAIL on a cohort broken in that way
python -m pytest -q tests
```

## Known limitations (recorded, not corrected)

- **T-cell-rich samples are favoured.** Requiring >= 100 T cells per donor and
  tissue selects for those samples.
- **Pairing excludes 17 donors.** On 10x, 14 LUAD donors qualify on the tumour
  side only and 3 on the normal side only.
- **Two tissue states, not three.** LUSC fails the donor floor (7 paired donors
  on 10x).
- **Tokenisation uses the loom input path.** The installed Geneformer h5ad path
  breaks under pandas 3. The h5ad-to-loom conversion is checked to be exact.
