# Balanced-donor workstream — Phase 1 survey and recommendation (2026-09-24)

- Rule applied: `PHASE0_SELECTION_RULE.md`, sha256 `58e9b153…1e3501`. It was
  registered and sent to god at 08:47Z, before any screening.
- Only structural fields were read: donor IDs, per-donor cell counts, and the
  origin, study, assay and cell-type labels. No expression values, DE,
  embeddings or published T-cell results were looked at.

## Sources screened, and how

1. **CELLxGENE Census 2025-11-08.** The cell-metadata table for all primary
   lung cells (6,167,731 cells; 748,205 T cells across 35 T-cell labels),
   grouped by dataset, disease and donor. Script: `census_lung.py` and
   `summarize.py`. Output: `census_2025-11-08_lung_tcell_donor_counts.csv`.
2. **Author metadata** for the candidates that the Census's standard columns
   cannot split into tumour vs normal. Read by HTTP byte-range from the
   published `.h5ad` file (metadata table only, no expression download).
   Script: `remote_obs.py`.

## Screen results. Criteria are applied in rule order; the reason given is the first failure.

| Candidate | Tumour donors / normal donors (>=100 T cells) | Result | First failure / note |
|---|---|---|---|
| **LuCA extended atlas, LUAD, tumour_primary vs normal_adjacent, 10x only** | **43 / 43, all paired within donor** | **PASS; R1 (paired) applies** | 5 studies, each contributing both groups: Leader_Merad 23, Kim_Lee 10, He_Fan 5, Lambrechts 3, Laughney 2. Assay is identical on both sides for every donor (10x 3' v2, and 10x 5' v1 in Leader). 41 of 43 have >= 300 T cells on both sides. |
| LuCA LUAD including BD Rhapsody (UKIM-V) | 55 paired | FAIL F4 | BD Rhapsody is not 10x. The 12 UKIM-V donors are dropped, which leaves the row above. |
| LuCA LUAD, tumour donors with no normal pair | 75 tumour / 58 normal | superseded | The 6 Singleron/InDrop tumour-only donors also fail F3 and F4. |
| LuCA LUSC, paired | 12 paired, of which 7 on 10x | FAIL A1 | 7 < 12 once F4 is applied. |
| LuCA NSCLC-NOS | 13 / 2 | FAIL A1 | |
| LuCA disease=='normal' as the normal group | 24 donors (>=100 T) | FAIL F3 | All from 4 non-cancer studies. None of those studies contributes a tumour donor (see below). |
| HLCA full (9f222629) | 0 tumour | FAIL (contrast) | No tumour group. |
| Treg TME LUAD (d41f45c1, HTAN) | 23 LUAD | FAIL (contrast/F2) | Author `Tissue Type` is Primary / Metastasis / Recurrence only. There is no normal tissue. |
| "Combined samples" (576f193c) | 14 / 4 | FAIL A1 | |
| Multi-tissue TME atlas, lung (b617ee1b) | 7 / 6 | FAIL A1 | |
| Kim_Lee_2020 alone (GSE131907) | 10 paired | FAIL A1 | Already inside LuCA. |
| ILD, COVID-19, BPD, COPD lung datasets | — | out of contrast | Not tumour vs normal. |

## Recommendation

- **Dataset:** LuCA extended atlas (CELLxGENE dataset version `33165751-ae33-4a65-94ee-52d9bc38f97e`,
  revised 2026-06-11).
  - Subset: LUAD, CD4/CD8 T cells, singlets, 10x assays.
  - Contrast: `tumor_primary` vs `normal_adjacent`.
  - The 43 paired donors are listed in `luca_luad_paired_10x_donors.csv`.
  - Analysis cohort: 100 cells per donor per side. Training cap: 300.
- **Raw counts:** present in `raw/X` (integer-valued; 17,764 genes).
- **Download size:**
  - **Possibly 0 GB.** The lab's `KD/data/nsclc/nsclc_integrated.h5ad` on
    thinkstation2 looks like this same atlas. *Inferred, not checked:* it has
    the same 17,764 genes, the same 58 `normal` donors and the same
    donor-ID style (`Adams_Kaminski_2020_001C`).
  - Phase 2 should first confirm on thinkstation2 that the file has
    `origin` and raw counts.
  - Otherwise, the full file is **17.6 GB** (17,616,543,376 bytes), well
    under the 100 GB threshold. Only ~0.25 GB of it is needed rows (about
    25,700 cells at the 300 cap). A partial read cannot be hashed against
    the source, so for integrity I prefer the full file plus our own sha256.
- **Where it lives:** a new top-level directory in `geneformer-lung-tcell`
  (for example `balanced_donor_luad/`). The data is lung T cells, so the repo
  name stays true. It also sits next to the s100 analysis layer, so there is
  no cross-repo copy. The alternative is `geneformer-nsclc-tcell`, which is
  where the NSCLC line "lives" since the split. It would need the analysis
  layer vendored in, and cross-repo duplication is what caused this week's
  mess.

## Classifier fine-tuning code: it exists

- **In git:** `sclc_validation/bf16_bench/run_finetune_316m.py` on
  `geneformer-lung-tcell` **origin/main** (added b33a002, 2026-09-18). It calls
  Geneformer `Classifier(...)`, then `.prepare_data(...)`, then
  `.validate(...)`: 1 epoch, lr 5e-5, batch 8, model V2. It is also on the
  `feature/s100-analysis-layer` branch and 5 other remote branches.
- **Why it was probably missed** (inferred): the local `main` checkout is
  9ec518c (2026-09-17), which predates the file. It is not on local `main`.
- **Outside git:** the original 104M runner,
  `KD/tcell_luad_lusc_normal_luscmax7000_finetune/scripts/run_finetune.py` on
  thinkstation2 (named in `current_workflow/METHODS.md` s.4).
- **Consequence:** Phase 4 is adapting existing code, not writing new code.
  It changes the label (tumour/normal), the split and the cap.

## Points in the brief I think are wrong or missing

1. **The biggest one, and it is outside this workstream: the existing lung
   classifier's "normal" class is study-confounded.**
   - In LuCA extended, all 58 `disease=='normal'` donors with T cells come
     from 6 studies: Adams_Kaminski, Habermann_Kropski, Madissoon_Meyer,
     Mayr_Schiller, Reyfman_Misharin and Vieira_Teichmann.
   - Those studies contribute **zero** LUAD or LUSC T cells.
   - The July classifier's normal class has 58 donors (37 + 9 + 12).
     *Inferred:* it is this set. So its LUAD/LUSC-vs-normal axis is partly a
     study, chemistry and site axis, and **donor balancing does not touch
     that**. A within-donor tumour vs adjacent-normal design does.
2. **"Mirroring the lung work" cannot be literal.** LUSC fails A1 on 10x (7
   paired donors). The balanced design is a two-state LUAD
   tumour-vs-adjacent-normal classifier, not LUAD/LUSC/normal.
3. **(a) must hold after the held-out split.**
   - Option 1, a donor-disjoint split: 43 paired donors split roughly
     22 / 6 / 15 train / eval / test. The ISP then runs on 15 paired test
     donors, so a signed test has a minimum p of about 6e-5.
   - Option 2, k-fold donor cross-fitting: all 43 donors enter the ISP, at
     the cost of k fine-tunes.
   - This is a Phase 4 costing choice.
4. **(b) and (c) are fixed by the equal cap, so they do not discriminate
   between datasets.** Imbalance in the raw data is large: max/min is 28x
   in tumour and 31x in normal among paired donors. After capping, max/min
   is 1.0.
5. **Pairing changes the Phase 5 statistics.** The natural test is a signed
   test over donors of the within-donor tumour-minus-normal shift. That is
   different from the two-group Mann–Whitney setup in the brief. It should be
   pre-registered as such.
6. **Selection bias from pairing:** on 10x, 57 LUAD donors qualify on the
   tumour side and 46 on the normal side. 14 are tumour-only and 3 are
   normal-only, and pairing excludes all 17. This is recorded as a limitation.
