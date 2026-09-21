# Committed run stats

Raw per-gene `targeted_panel/stats/overexpress/*.csv` for the six GPU-run
arms behind the precision comparisons in [`../RESULTS_BF16.md`](../RESULTS_BF16.md)
and [`../RESULTS_BF16_316M.md`](../RESULTS_BF16_316M.md), added 2026-09-22
(card `bf16-run-stats-provenance-20260922`) to close a reproducibility gap:
these tables previously lived only under thinkstation1's gitignored `runs/`
tree, so a reader with only the repository could read the two results
documents' conclusions but not independently recompute them.

| Directory | Arm | Used in |
|---|---|---|
| `fp32_baseline/` | 104M, fp32, cold cache | RESULTS_BF16.md gate (accuracy reference) |
| `fp32_repeat/` | 104M, fp32, warm cache | RESULTS_BF16.md gate (noise floor + speed/energy reference) |
| `bf16/` | 104M, bf16 | RESULTS_BF16.md gate; RESULTS_BF16_316M.md panel comparison |
| `316m_bf16/` | 316M, bf16, full 50-gene panel | RESULTS_BF16_316M.md panel comparison |
| `316m_canary_fp32/` | 316M, fp32, 10-gene canary | RESULTS_BF16_316M.md precision canary |
| `316m_canary_bf16/` | 316M, bf16, 10-gene canary | RESULTS_BF16_316M.md precision canary |

**Provenance**: copied verbatim (not regenerated) from the live GPU-run
output on thinkstation1 at `sclc_validation/bf16_bench/runs/<arm>/targeted_panel/stats/`
on 2026-09-22, byte-identical (`diff -rq`) to the source at copy time.
`CHECKSUMS.sha256` records a sha256 for every file here.

**Reproduce the published numbers** from these committed copies alone
(no GPU needed), using the already-committed `compare_runs.py`:

```bash
# RESULTS_BF16.md's gate table (0.9998 / 20/20 / 0.9933 / FAIL):
python3 ../compare_runs.py panel \
  --a-stats fp32_baseline/targeted_panel/stats \
  --b-stats bf16/targeted_panel/stats \
  --floor-stats fp32_repeat/targeted_panel/stats \
  --perturb-types overexpress --out /tmp/verify_104m_panel.json

# RESULTS_BF16_316M.md's panel comparison (rho 0.489, top-20 9/20):
# --floor-stats is required by the CLI but its output is not meaningful
# here -- the gate does not apply across architectures (see RESULTS_BF16_316M.md).
python3 ../compare_runs.py panel \
  --a-stats bf16/targeted_panel/stats \
  --b-stats 316m_bf16/targeted_panel/stats \
  --floor-stats fp32_repeat/targeted_panel/stats \
  --perturb-types overexpress --out /tmp/verify_316m_vs_104m.json

# RESULTS_BF16_316M.md's precision canary (rho 0.99867, sign agreement 1.0, max|delta| 0.00072):
python3 ../compare_runs.py panel \
  --a-stats 316m_canary_fp32/targeted_panel/stats \
  --b-stats 316m_canary_bf16/targeted_panel/stats \
  --floor-stats fp32_repeat/targeted_panel/stats \
  --perturb-types overexpress --out /tmp/verify_canary.json
```

All three were re-run against these exact committed files while preparing
this PR and reproduced the published numbers exactly.
