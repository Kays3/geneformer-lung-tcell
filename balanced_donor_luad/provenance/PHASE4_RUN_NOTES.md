# Phase 4 run notes (fine-tunes, thinkstation1)

## Attempt 1 — crashed, archived, not used
- Started 2026-09-24T11:02:07Z at code 9097530. Fold 0 trained and scored its held-out donors, then crashed at 11:39:38Z while writing `fold_record.json`: `float()` was applied to `Classifier.validate`'s per-k-split list for `macro_f1`/`acc`. The driver's `set -e` stopped the run before fold 1.
- Archived unchanged at `~/workspace/balanced_donor_luad/phase4_attempt1_crashed/` (thinkstation1). None of its outputs enter the gate or any later phase.
- GPU cost: 0.63 GPU-h, reported separately in the GPU-hour account.
- The watcher missed the crash because `pgrep -f run_phase4.sh` matched the watcher's own SSH command line. The GPU was idle 11:39Z–17:26Z.

## Attempt 2 — the registered run
- Started 2026-09-24T17:26:18Z at code f9ec970; all 5 folds from scratch.
- `git diff 9097530..f9ec970` covers 5 commits and 16 files. Its only removed lines are the two record-writing lines in `finetune_fold.py`; everything else is additions. The diff was reviewed line by line by god (full text in the hive, `phase4_fix_diff.txt`).
- Inputs re-verified before launch: the tokenised dataset passes `sha256sum -c` against `tokenized_dataset_sha256.txt`, and the base weights hash is 965cecce….

## Hugging Face `datasets` cache files in the input dataset folder
`balanced_donor_luad_pool300.dataset/` on thinkstation1 contains `cache-*.arrow` files. They are not covered by `tokenized_dataset_sha256.txt`:
- `3d9b3970…`, `4f05f06b…`, `e31a8766…` (09:41:38Z) and `a386c97f…` (09:45:33Z): written before Phase 4, by earlier filter/map calls on this dataset.
- `81bb3b8f…` (11:02:11Z): written by attempt 1's fold-0 filter.
- No cache file was written after 17:26Z as of 17:45Z. That is consistent with attempt 2 reusing `81bb3b8f…` for fold 0, but was not verified directly.

The cache was left in place: clearing it under a live run is a larger risk than keeping it. A cache file is a fingerprinted, deterministic transform of the verified arrow, and the attempt-1-to-2 code change touches no map or filter function.

**Scope of any reproducibility statement from the attempt-1 vs attempt-2 fold-0 comparison.** Because preprocessed input may be reused, a match shows only that the TRAINING step reproduces given identical preprocessed input. It is NOT evidence that the tokenise/filter stage reproduces from raw data. Deterministic kernels are off (no `use_deterministic_algorithms`, no `full_determinism`); the seeds are 43 (trainer) and 42 (Geneformer shuffle).

## Same-seed reproducibility check: attempt 1 vs attempt 2, fold 0 (read fixed before the comparison)
Pre-fixed reading: bitwise identical would confirm the code claim by measurement. A mismatch is NOT evidence of a code change, because deterministic kernels are off; the code question is settled by the diff (reviewed: only record lines removed). Any donor crossing 0.5 is reported and not investigated before all five folds exist.

Result, on the same 1,800 held-out test cells (9 donors × 100 tumour + 100 normal), with identical cell order and identical labels in both files:
- **Not bitwise identical.** Every cell's logits differ: median |Δ| 0.47, 99th percentile 1.76, max 2.89.
- (a) Same predicted label: 1,727 / 1,800 cells (95.9%).
- (c) Per-donor balanced-accuracy change: -0.045 to +0.015. No donor crosses 0.5.
- Eval-donor metrics: attempt 1 macro-F1 0.894, accuracy 0.894; attempt 2 macro-F1 0.899, accuracy 0.899.

Reading: GPU training on this stack is not bit-reproducible run to run at a fixed seed. The size of that run-to-run variation is recorded here. It is not a code difference (see the diff). Only attempt 2 enters the gate.

**Note on times (appended 2026-09-24T18:13:20Z):** in these notes and in hive messages, a reported time is an OBSERVATION time unless it is marked as taken from the filesystem or a commit. Example: attempt 2 fold 0 was reported as landing at "18:05Z", which is when it was observed. Filesystem: `fold_record.json` mtime 2026-09-24T18:03:46.827Z. For any time difference, use filesystem or commit times.
