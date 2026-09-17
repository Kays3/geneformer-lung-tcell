# BF16 vs FP32 in-silico perturbation: replication on our T-cell pipeline

**Status: SKELETON -- no GPU arms have run yet.** This document is a
placeholder until the sizing/checkout questions on card
`geneformer-bf16-replication-20260917` are resolved and god authorizes GPU
time (checkpoint 2 in `hive/handoffs/bf16-isp-replication-plan-20260917.md`).
Every number below must eventually come from a generated `compare_runs.py`
JSON or a `power_sample.sh` summary -- none are hand-typed.

## Provenance notes

- **Geneformer checkout pin**: runs use `GENEFORMER_ROOT=/home/kaisar/workspace/geneformer-uv-starter/Geneformer` (private copy, pinned exactly at `f45a6c7de57ff07f946f146c254da02a90e2cdf5`, per god's step-0 decision -- the shared `/srv/lab/geneformer` path (what both runners import by default when `GENEFORMER_ROOT` is unset) is deliberately NOT patched or pointed at, since it is a shared host resource. **Discrepancy on record**: `/srv/lab/geneformer`'s live symlink target is currently `04c2b2e`, 4 commits ahead of the documented `f45a6c7` pin. Diffed: only `README.md` and `geneformer/mtl/{model,train}.py` differ between the two commits -- `perturber_utils.py`/`emb_extractor.py`/`in_silico_perturber.py` (everything this task's patch and ISP runs touch) are byte-identical, so past and present ISP results are unaffected either way. The documentation pin and the live default import path still disagree and someone should eventually reconcile them; out of scope for this task.
- **Peak GPU memory**: `nvidia-smi` reports `Memory-Usage: Not Supported` on this GB10 (unified-memory architecture, confirmed live via `nvidia-smi` on thinkstation1) -- `power_sample.sh` therefore does not attempt a memory column. Peak memory instead comes from `torch.cuda.max_memory_allocated()`, called inside `run_gene()`/`run_one()` right after each unit and reset for the next, recorded as `peak_gpu_mem_gib` in every completion marker.

## Sizing amendment 2026-09-17 (dated, before the bf16 arm ran)

Calibration (single-gene PDCD1, fp32, private `f45a6c7` checkout, `nproc=1` -- see the plan's own dated amendment in `hive/handoffs/bf16-isp-replication-plan-20260917.md` for the full derivation):

| perturb_type | source | source cell pool | elapsed_seconds | peak_gpu_mem_gib |
|---|---|---:|---:|---:|
| delete | normal | 566 | 22.35 | 11.05 |
| overexpress | normal | 566 | 67.30 | 39.93 |
| overexpress | sclc | 2,424 | 169.19 | 39.93 |
| overexpress | luad | 6,387 | 375.51 | 39.93 |

The naive uncapped arm (both types, all 3 sources, 50 genes) projects to well over 2.5h just from overexpress on luad's 6,387 cells. Applying god's lever priority: (a) the 50-gene panel is kept in full; (b) a `max_ncells` cap alone (keeping both perturb types) cannot reach the target -- fixed per-unit overhead alone, summed over 300 units, already exceeds the 2.5h budget; so (c) is also needed. **Decision: overexpress only, `--max-ncells 300` applied identically to every source and every arm.** Projected: ≈150 units × ~59.4s ≈ **2.47h/arm**, ≈7.4h for the three 104M arms sequential. Also found and fixed two pre-existing bugs (unrelated to bf16) that blocked the calibration unit from running at all: a bad `ANALYSIS_ROOT` path and `NPROC=4` forking after CUDA init -- see PR commit `b011667`.

## Objective

Replicate the colleague's finding (petadimensionlab/Geneformer fork:
bf16 matches fp32 for Geneformer ISP -- gene-rank Spearman 0.98, top-20
identical -- while being 2.28x faster and ~74% cheaper in energy) on our
own `geneformer-lung-tcell` pipeline (fine-tuned Geneformer-V2-104M
3-class classifier, thinkstation1 GB10).

## Pre-registered acceptance gate

Frozen before the bf16 arm runs (see `compare_runs.py`'s `RHO_MIN`,
`TOPN`, `TOPN_MIN_OVERLAP` constants -- any change after the bf16 arm has
run is a dated amendment below, not a silent edit):

- Per-gene Spearman rho (bf16 vs fp32 baseline) >= 0.95 on both
  `overexpress_shift` and `delete_shift`.
- Top-20 overlap >= 19/20.
- Sign agreement on all genes whose |shift| exceeds the fp32-vs-fp32
  noise floor.
- bf16-vs-fp32 |delta shift| distribution reported against the fp32-vs-fp32
  floor, split into signal-gene and near-zero-effect-gene subsets.
- T4 arm: no matched-null test flips its significance conclusion at
  alpha = 0.05.

## Arms run

Per the sizing amendment above: **overexpress only** (delete dropped for all
three precision arms), **`--max-ncells 300`** applied identically to every
arm.

| Arm | dtype | run-tag | perturb type | max-ncells | Status |
|---|---|---|---|---|---|
| fp32 baseline | fp32 | `fp32_baseline` | overexpress | 300 | **done** (exit 0, wall 6671.61s) |
| fp32 repeat (noise floor) | fp32 | `fp32_repeat` | overexpress | 300 | running |
| bf16 | bf16 | `bf16` | overexpress | 300 | not run |
| T4 program phase (12/273 units), fp32 | fp32 | `t4_fp32` | overexpress | n/a | not run (optional) |
| T4 program phase (12/273 units), bf16 | bf16 | `t4_bf16` | overexpress | n/a | not run (optional) |

## Gate results (targeted panel, `compare_runs.py panel`)

_Pending. Command once the three arms complete:_

```
python3 compare_runs.py panel \
  --a-stats runs/fp32_baseline/targeted_panel/stats \
  --b-stats runs/bf16/targeted_panel/stats \
  --floor-stats runs/fp32_repeat/targeted_panel/stats \
  --out runs/panel_verdict.json
```

| Perturb type | rho | top-20 overlap | sign agreement (signal genes) | verdict |
|---|---|---|---|---|
| overexpress | -- | -- | -- | -- |

(No `delete` row: dropped for all three precision arms per the sizing
amendment above.)

## Gate results (T4 program phase, `compare_runs.py t4`, optional arm)

**Dated amendment 2026-09-17 (god's PR review, item 1): gate re-scoped to
shift-vector correlation only.** The optional T4 arm is program phase only
(12 of 273 units); the matched-null significance-flip check needs the
separate 240-unit null phase, out of scope for this arm. `compare_runs.py
t4` now accepts `--a-null`/`--b-null` as optional (both omitted together)
and falls back to a correlation-only gate: Spearman rho on the scope=="all"
shift vector must be >= 0.95 (the same `RHO_MIN` as the panel gate), no
separate matched-null check for this arm. This still tests something real
-- whether bf16 preserves the set-level shift ranking on a program-level
(multi-gene) perturbation, a different workload shape than the per-gene
panel -- just not the full G5 gate.

_Pending the arm itself. Command once run:_

```
python3 compare_runs.py t4 \
  --a-shift runs/t4_fp32/.../t4_shift_summary.csv \
  --b-shift runs/t4_bf16/.../t4_shift_summary.csv \
  --out runs/t4_verdict.json
```

| Spearman rho (shift vector) | verdict |
|---|---|
| -- | -- |

## Speed / energy / heat

Per arm, from `power_sample.sh`'s `<label>.power.json` plus each unit's
`elapsed_seconds` / `peak_gpu_mem_gib` markers:

| Arm | wall time | mean power (W) | energy (Wh) | peak temp (C) | peak GPU mem (GiB) |
|---|---|---|---|---|---|
| fp32 baseline | 6671.61s (1.85h) | 63.523 | 117.723 | 84.0 | see per-unit markers (39.93 GiB at calibration) |
| fp32 repeat | -- | -- | -- | -- | -- |
| bf16 | -- | -- | -- | -- | -- |

(fp32 baseline figures are from `runs/fp32_baseline/power/fp32_baseline.power.json`, generated by `power_sample.sh`, not hand-typed.)

Speedup = fp32 baseline wall time / bf16 wall time. Colleague-comparable
headline: rho, top-20 overlap, speedup, Wh per arm.

## Noise floor (fp32 vs fp32 repeat)

_Pending -- bounds how tightly bf16 deltas should be expected to track
fp32, since the fp32-vs-fp32 rerun is not perfectly bitwise-identical for
reasons unrelated to precision (e.g. any within-run nondeterminism); the
colleague measured the same floor for the same reason._

## Recommendation

_Pending arms + gate. One of: switch future ISP runs to bf16 (yes/no),
with the null and floor named for every statistic cited, per this
project's standing rule against hand-waved comparisons._
