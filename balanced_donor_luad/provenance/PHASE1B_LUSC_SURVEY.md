# Phase 1b — LUSC survey under Amendment 3.3 (2026-09-24)

- **Fields used:** structural only, meaning labels, donor IDs, per-donor cell
  counts, chemistry and access. No expression.
- **Sources:**
  - CELLxGENE Census 2025-11-08, every LUSC/NSCLC dataset;
  - a repository and literature sweep (GEO, ArrayExpress/BioStudies, HTAN,
    ScienceDB, papers' sample tables): 19 candidates;
  - direct metadata reads (HTTP byte-range, cell metadata only) of the
    top candidate.

## Result: the outcome depends on how one registered phrase is read

Rule F3-LUSC says a LUSC-contributing study must also contribute qualifying
LUAD tumour donors "**on the same chemistry**". Rule F4 (Phase 0) names
chemistry at the family level: "10x 3′ or 10x 5′". The phrase can mean the
**family** or the **exact version**, and the donor count lands on either
side of the floor of 12.

| LUSC tumour donors (>= 100 CD4/CD8-type T cells, 10x, primary tumour) | Family reading (10x 3′ / 5′) | Version reading (exact, e.g. 3′ v2) |
|---|---|---|
| LuCA Leader_Merad_2021 (LUAD donors on 3′ v2: 26) | 5 | 5 |
| LuCA Lambrechts_Thienpont_2018 (LUSC: 2 on 3′ **v1**, 1 on 3′ v2; LUAD: 3 on 3′ v2) | 3 | 1 |
| LuCA Goveia_Carmeliet_2020 (no LUAD donors, so it **fails F3-LUSC**) | 0 | 0 |
| De Zuani 2024, E-MTAB-13526 (10x 3′ v3.1; LUAD: 5 paired tumour+normal donors) | 4 | 4 |
| **Total** | **12, PASS (floor 12; below the preferred 16)** | **10, FAIL (stop rule)** |

**I have not chosen between the readings.** Choosing after seeing which
side of 12 each lands on is the post-hoc freedom the registration exists to
prevent. It goes to god and the human.

## De Zuani 2024 (E-MTAB-13526), verified by me from its cell metadata

- **Access:** open BioStudies; per-sample raw `mtx` plus annotated `.h5ad`
  files; untreated patients (per the paper, as reported by the survey).
- **Fraction used:** only the whole-tissue `CD235a-` fraction, to match
  LuCA. The `Myeloid` / `MDSC`-sorted fractions contain few T cells.
- **Histology:** harmonised by a rule stated here.
  - LUSC = {Squamous cell carcinoma, Squamous carcinoma, Squamous cell lung
    cancer, Squamous cancer}.
  - LUAD = {Adenocarcinoma, TTF1 +ve lung adenocarcinoma, lung
    adenocarcinoma, Mucinous adenocarcinoma}.
  - Excluded: "Lung cancer", "Presumed lung cancer", "NSCLC", NA, and
    **squamous dysplasia** (pre-invasive).
- **T cells:** author labels containing "T cells", excluding Tregs, cycling,
  γδ and NK. That matches LuCA's `T cell CD4/CD8` choice.
  - **Caveat:** De Zuani's labels are functional states (naive, exhausted,
    "downregulated"), not an explicit CD4/CD8 split, so the match to LuCA's
    definition is approximate.
- **Counts:**
  - LUSC tumour donors >= 100: **4** (203, 5,631, 6,824, 11,298), of which 3
    also have qualifying background lung.
  - LUAD paired donors: **5** (tumour 108–5,248).
- **Rule F4 consequence:** 10x 3′ v3.1 must appear in all three classes, so
  De Zuani's **5 LUAD paired donors must join the LUAD↔normal cohort**,
  taking it from 43 to 48.
  - The existing 43 donors' draws are unchanged, because the draw is seeded
    per donor.
  - This extends Stage L2's "LUAD part unchanged": the extension is
    required by F4, not chosen.
- **Download:** only the per-sample `mtx` for these 9 donors' `CD235a-`
  samples (small, to be measured before pulling). The annotated `.h5ad`
  files (55 + 42 GB) are needed only for the cell labels. Their metadata
  table has already been read by byte-range, so no full download is needed.

## Other candidates (first criterion failed)

| Candidate | Reason |
|---|---|
| GSE243013 (Liu 2025, Cell): 180 LUSC, 63 LUAD, open | **F4.** No adjacent normal on its chemistry, so that chemistry is absent from the normal class. Also **all post-neoadjuvant**, which is not a registered criterion but a treatment confound between LUSC/LUAD and the untreated LuCA cohort. Flagged, not used. |
| ScienceDB 02028 (Wang 2022): 9 LUSC, 10 LUAD | **Unverified.** About 1.2 TB (> 100 GB, so god first), chemistry version and matrix availability unverified. Not screened further without approval. |
| Zhang 2022 (West China): 11 LUSC | **F6 unverified.** Access route unclear, probably controlled. |
| GSE207422, GSE229353 | Mostly post-therapy; BD / small. |
| GSE194070 | **F3-LUSC:** no LUAD. |
| GSE99254, GSE162498, GSE139555, GSE117570, GSE200972 | **A1:** fewer than 5 LUSC donors, or sorted T cells / Smart-seq2 (F4). |
| GSE179994, GSE176021 | Histology split unverified; post-therapy / sorted. |
| GSE148071, E-MTAB-6149, GSE154826 | Already in LuCA; no larger public LUSC release found. |
| HLCA (3 LUSC), CxG FFPE study (2 LUSC) | **A1.** |

## What each reading leads to

- **Family reading:** 3-class goes ahead with exactly 12 LUSC donors from 3
  studies.
  - About 2–3 LUSC donors are held out per fold.
  - It rests on a between-donor contrast (Amendment 3.6), and the De Zuani
    label mapping is approximate.
- **Version reading:** the stop rule applies. The human is asked before any
  relaxation, and the default is the registered 2-class run, unchanged.
