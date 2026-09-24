# Calibration conditions (companion to calibration_rates_316m_bf16.json)

A measurement is the number plus the conditions it was taken under. These are
the conditions for the 2026-09-24 calibration.

## Host and precision

- **Host and GPU:** thinkstation1, NVIDIA GB10, driver 595.84.
- **Model and precision:** Geneformer-V2-316M, bf16. Vendored Geneformer
  f45a6c7 + bf16 export patch, `sclc_analysis` venv.

## Timing

- **Run window (UTC):** started 2026-09-24T09:43:39Z, ended 09:46:15Z.
- **Hard stop:** the OS-level `timeout 540` **did NOT fire**.
  - The process exited with code 0.
  - `total_wall_seconds` = 153.8, against a 540 s limit.
  - The log contains no kill or termination.
  - So the rates are **completed measurements, not lower bounds**.

## What else was running

- **GPU contention:**
  - **0 compute processes on the GPU at launch**, checked with
    `nvidia-smi --query-compute-apps` immediately before starting.
  - GPU residency was **not sampled during the run**, so any contention
    during it is **unknown**.
  - 0 compute processes again when re-checked afterwards.
- **CPU contention:**
  - A single-core bash busy-loop (TE-OA, not ours) was observed at about
    09:15Z.
  - Whether it was still running during the window is **unknown**; it was
    gone when re-checked later.
  - It could only affect the CPU-side data steps, not GPU compute.

## Why the rates are conservative (upper bounds for costing)

- **ISP gene:** B2M is expressed in effectively every cell, so deletion
  touched all 300 cells. That is the worst case for cost.
- **Model loading:** each timed ISP call included a model load. Real runs
  load once per gene over more cells.
- **Fine-tune set:** 900 cells, median 745 tokens. The full cohort's
  training pool has a median of 755 tokens.

## Failed starts before the successful run

There were three, all stopped before any training step, for about 16 s in
total:
1. `prepare_data` output name (`_labeled_train/_test`, not `_labeled`);
2. `/srv/lab` 316M weights file is a git-LFS pointer;
3. the retry.

**Total GPU time charged to the approval:** about 0.05 of 0.15 GPU-h.
