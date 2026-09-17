# BF16 vs FP32 in-silico perturbation: replication on our T-cell pipeline

**Status: SKELETON -- no GPU arms have run yet.** This document is a
placeholder until the sizing/checkout questions on card
`geneformer-bf16-replication-20260917` are resolved and god authorizes GPU
time (checkpoint 2 in `hive/handoffs/bf16-isp-replication-plan-20260917.md`).
Every number below must eventually come from a generated `compare_runs.py`
JSON or a `power_sample.sh` summary -- none are hand-typed.

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

| Arm | dtype | run-tag | Status |
|---|---|---|---|
| fp32 baseline | fp32 | `fp32_baseline` | not run |
| fp32 repeat (noise floor) | fp32 | `fp32_repeat` | not run |
| bf16 | bf16 | `bf16` | not run |
| T4 program phase (12/273 units), fp32 | fp32 | `t4_fp32` | not run (optional) |
| T4 program phase (12/273 units), bf16 | bf16 | `t4_bf16` | not run (optional) |

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
| delete | -- | -- | -- | -- |
| overexpress | -- | -- | -- | -- |

## Gate results (T4 program phase, `compare_runs.py t4`, optional arm)

_Pending._

## Speed / energy / heat

Per arm, from `power_sample.sh`'s `<label>.power.json` plus each unit's
`elapsed_seconds` / `peak_gpu_mem_gib` markers:

| Arm | wall time | mean power (W) | energy (Wh) | peak temp (C) | peak GPU mem (GiB) |
|---|---|---|---|---|---|
| fp32 baseline | -- | -- | -- | -- | -- |
| fp32 repeat | -- | -- | -- | -- | -- |
| bf16 | -- | -- | -- | -- | -- |

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
