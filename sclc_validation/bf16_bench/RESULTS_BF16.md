# BF16 vs FP32 in-silico perturbation: replication on our T-cell pipeline

**Status: ALL THREE 104M ARMS COMPLETE. Frozen-gate verdict: FAIL** (sign
criterion, overexpress; see Gate results below). Ruling on whether this
FAIL reflects a gate degeneracy (zero floor from a bitwise-identical
fp32-vs-fp32 rerun) rather than a real bf16 defect is pending Pam as
domain reviewer -- **the implementer (Kevin) is not ruling on the gate his
own result needs**, per god's explicit instruction. No bf16-yes/no
recommendation until that ruling lands. Every number below comes from a
generated `compare_runs.py` JSON (`runs/panel_verdict.json`) or a
`power_sample.sh` summary -- none are hand-typed.

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
| fp32 repeat (noise floor) | fp32 | `fp32_repeat` | overexpress | 300 | **done** (exit 0, wall 4944.71s) |
| bf16 | bf16 | `bf16` | overexpress | 300 | **done** (exit 0, wall 3393.09s) |
| T4 program phase (12/273 units), fp32 | fp32 | `t4_fp32` | overexpress | n/a | not run (optional) |
| T4 program phase (12/273 units), bf16 | bf16 | `t4_bf16` | overexpress | n/a | not run (optional) |

## Gate results (targeted panel, `compare_runs.py panel`)

Run to the canonical path (`runs/panel_verdict.json`), gate constants and
signal-gene definition unchanged from freeze:

```
python3 compare_runs.py panel \
  --a-stats runs/fp32_baseline/targeted_panel/stats \
  --b-stats runs/bf16/targeted_panel/stats \
  --floor-stats runs/fp32_repeat/targeted_panel/stats \
  --perturb-types overexpress \
  --out runs/panel_verdict.json
```

| Perturb type | rho | top-20 overlap | sign agreement (signal genes) | verdict |
|---|---|---|---|---|
| overexpress | 0.9998 | 20/20 | 0.9933 (298/300) | **FAIL** |

(No `delete` row: dropped for all three precision arms per the sizing
amendment above.)

**Dated amendment 2026-09-18 (gate degeneracy found, ruling pending
Pam): the sign-agreement criterion failed on a floor technicality, not on
a large disagreement.** The `fp32_baseline` vs `fp32_repeat` noise-floor
comparison came out **bitwise identical** (`max_abs_diff = 0.0`), so
`floor_threshold = 0.0` and every one of the 300 (gene, comparison) rows
counts as a "signal gene" -- including near-zero-shift genes where a sign
flip carries no real magnitude. 2/300 rows flip sign between fp32 and
bf16, both far below the panel's mean |shift| of 0.00917:

| Gene | comparison | fp32 `Shift_to_goal_end` | bf16 `Shift_to_goal_end` | \|value\| vs mean\|shift\| |
|---|---|---:|---:|---|
| TPSB2 | normal_to_sclc | -0.0000760 | +0.0000094 | ~120x smaller |
| MMP12 | sclc_to_normal | -0.0001084 | +0.0000781 | ~10x smaller |

Both rows independently re-verified against the raw stats CSVs (not just
the JSON summary). Per gate rules and god's explicit instruction, this
FAIL is recorded as-is -- **the gate constants and signal-gene definition
were NOT changed to make this pass.** Whether a zero-floor degeneracy like
this should exempt near-zero-shift genes from the sign criterion is a gate
methodology question routed to Pam as domain reviewer, not decided here.

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

**Dated amendment 2026-09-18 (god's review): fp32_repeat, not fp32_baseline,
is the speed/energy reference.** fp32_baseline's wall time (6671.61s) is
~27 minutes longer than fp32_repeat's (4944.71s) on identical work (same
dtype, same flags, same 150 units) -- the difference is a one-time
`_shared/targeted_panel_sources` cache build that only fp32_baseline pays
(fp32_repeat and bf16 both hit the already-warm cache). Using
fp32_baseline as the speed/energy denominator would therefore overstate
bf16's speedup. **fp32_repeat's wall time and Wh are the fp32
speed+energy reference** (warm-cache vs warm-cache, matching how the
colleague's 2.28x claim was measured); fp32_baseline remains the
ACCURACY reference exactly as already gated (`compare_runs.py panel
--a-stats runs/fp32_baseline/...`), since the gate compares per-gene
statistics, not wall time, and is unaffected by the cache-build asymmetry.

Per arm, from `power_sample.sh`'s `<label>.power.json` plus each unit's
`elapsed_seconds` / `peak_gpu_mem_gib` markers:

| Arm | wall time | mean power (W) | energy (Wh) | peak temp (C) | peak GPU mem (GiB) |
|---|---|---|---|---|---|
| fp32 baseline (accuracy ref; includes one-time cache build) | 6671.61s (1.85h) | 63.523 | 117.723 | 84.0 | 39.93 (calibration) |
| fp32 repeat (**speed/energy ref**) | 4944.71s (1.37h) | 80.235 | 110.2056 | 84.0 | 40.42 (max across units) |
| bf16 | 3393.09s (0.94h) | (max 89.95W) | 35.794 | 78.0 | 20.25 (max across units) |

(Figures are from `runs/<tag>/power/<tag>.power.json`, generated by
`power_sample.sh`, and from the max `peak_gpu_mem_gib` across each arm's
completion markers under `runs/<tag>/targeted_panel/raw/`, not hand-typed.)

**Speedup = fp32_repeat wall time / bf16 wall time = 4944.71 / 3393.09 =
1.457x.** **Energy ratio = fp32_repeat Wh / bf16 Wh = 110.2056 / 35.794 =
3.079x (bf16 uses ~67.5% less energy).** Peak GPU memory is also roughly
halved (40.42 -> 20.25 GiB), consistent with the earlier n=1 smoke-test
bonus data point (~2.07x speedup, ~half memory on one unit).

Colleague-comparable headline: **rho 0.9998, top-20 20/20, 1.457x
speedup, 3.079x less energy** vs the colleague's reported 0.98 / 20/20 /
2.28x / ~74% less energy -- same direction and same order of magnitude,
smaller speedup/energy factor on this workload and this GB10 (unified
memory, different bottleneck profile than the colleague's hardware).
Frozen-gate accuracy verdict is FAIL on the sign-agreement criterion (see
above), pending Pam's ruling on the zero-floor degeneracy.

## Noise floor (fp32 vs fp32 repeat)

From `runs/panel_verdict.json`'s `fp32_vs_fp32_floor` block: fp32_baseline
and fp32_repeat are **bitwise identical** on all 300 (gene, comparison)
rows -- spearman 1.0, sign agreement 1.0, `max_abs_diff = 0.0`,
`mean_abs_diff = 0.0`. This is tighter than anticipated at freeze time
(the plan expected some within-run nondeterminism, as the colleague also
found on their hardware); on this GB10, targeted-panel ISP inference
appears fully deterministic run-to-run for identical fp32 inputs. This
zero floor is the direct cause of the gate's sign-criterion FAIL above --
with any nonzero floor, both flipped genes (magnitudes ~1e-5 and ~1e-4)
would very likely have fallen below the near-zero-effect-gene threshold
and been excluded from the sign check entirely.

## Recommendation

**Frozen-gate verdict: FAIL** (overexpress panel, sign-agreement
criterion: 298/300 = 0.9933, below the implicit 1.0 requirement once the
noise floor collapsed to zero -- see Gate results and Noise floor above).
This FAIL is recorded as-is; the gate constants and signal-gene
definition were not touched to change the outcome, per instruction, since
the implementer (Kevin) must not rule on the gate his own result needs.

**No bf16-yes/no recommendation is made here.** A ruling is pending from
Pam (domain reviewer) on whether the two sign flips -- both magnitude
~1e-5/1e-4, ~10-120x smaller than the panel's mean |shift| of 0.00917,
surfaced only because the fp32-vs-fp32 floor came out exactly zero --
represent a real bf16 defect or a gate methodology gap (no allowance for
near-zero-effect genes when the empirical floor is degenerate). Every
other measured criterion passed cleanly (rho 0.9998 >> 0.95, top-20
20/20) and the speed/energy profile replicates the colleague's direction
of effect (1.457x faster, 3.079x less energy, ~half peak GPU memory).
This section will be updated once Pam's ruling lands.
