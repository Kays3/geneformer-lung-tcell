# S100 ISP analysis layer

Built 2026-09-23 (card `s100-isp-execution-20260922`) against
`hive/reports/s100-isp-design-20260922.md` and `retained_rows_spec_20260922.md`.
No GPU. No dependency on `geneformer`/`torch` anywhere in this directory --
every module here is plain `pandas`/`numpy`/`scipy` and is fully testable
without a GPU host, a real dataset, or the vendored Geneformer checkout.

## What's built (fully specified, tested against synthetic fixtures)

| Module | What it does |
|---|---|
| `retained_rows.py` | Builds the planned 96-row Cartesian table (12 genes x 2 sources x 2 goals x 2 operations), left-joined against eligibility manifests, run-completion markers, no-op results, and ISP stats CSVs. A missing/not-yet-run/not-estimable row is always kept, never dropped, with a machine-readable `status`/`status_reason`. Computes `donor_balanced_shift` directly from the raw per-cell pickle + the paired-eligible manifest's donor column -- **not** from `InSilicoPerturberStats`' own `Shift_to_goal_end` column, which is a plain cell-weighted mean, not the donor-weighted mean the design requires. |
| `ambient_stats.py` | `E[g,c]`, `Q[g,c]` (midrank percentile within matched controls), per-stratum rollup, the primary Spearman test with an **exact** (not asymptotic) label-permutation p-value enumerating all `n!` permutations, the zero-GPU stability gates (120 leave-one-control-out recomputations, 10,000 within-stratum bootstrap redraws at seed `20260922`), and the `0.00072` numerical effect floor. |
| `matched_controls.py` | `synthetic_control_table()` -- a structurally-valid **fake** 20-control-per-stratum table (fabricated `ENSG9SYNTH...` ids) used only to prove `ambient_stats.py` is correct before the real controls exist. |
| `provenance_utils.py` | `sha256_file()` -- every file this layer reads gets hashed from its own bytes at the moment it's read, same principle as `run_targeted_panel.py`'s self-hash/input-hash (2026-09-23). |

Run the tests (no GPU, ~5s total, the slow part is the honest n=10
factorial-permutation check):

```
python3 test_retained_rows.py
python3 test_ambient_stats.py
```

## What's blocked (a stub, not a guess)

`matched_controls.build_matched_control_table()` raises `NotImplementedError`
unconditionally. Two definitional questions are open with the human/Pam
(gate 3(b), not an implementation gap):

1. Which file/column defines "median token rank" for a gene -- the
   existing ambient-risk table has a detection fraction but no rank column.
2. Whether the matching tolerance is computed per source-state or globally.

Do not implement it against a guessed answer. Everything in
`ambient_stats.py` works identically against the real control table once
built -- that is the entire reason it was built against
`synthetic_control_table()` now rather than after.

## Gaps surfaced while building this -- RULED 2026-09-23

Three gaps were surfaced rather than resolved unilaterally; Michael ruled
on all three. Full text in `retained_rows_spec_20260922.md`'s dated
amendment; summary here:

1. **`status` gets `not_run`, outside the spec's six-value enum** (every
   listed value describes a completed outcome; none describes "eligible,
   not yet run"). **Ruled: stands** -- it labels a pre-run state and
   changes no threshold, statistic, or interpretation.
2. **`donor_balanced_shift` must come from the raw pickle + manifest, never
   from `InSilicoPerturberStats`' `Shift_to_goal_end`** (a cell-weighted
   mean, not the donor-weighted one the design requires). **Ruled: this
   was the most important catch in the whole layer** -- quantified on LUAD
   S100A2 (four donors at 40/4/4/25 cells), cell-weighting would give the
   40-cell donor 2.19x its intended weight, on the only ambient-flagged
   gene surviving that arm at all. The design doc names this exact failure
   mode as one that already happened once in the existing lung analysis.
3. **`midrank_percentile()`'s exact formula, and whether `status` gets
   overwritten when a stratum's control set fails the 20-count gate.**
   **Ruled: both stand as implemented, and they are one decision, not
   two** -- Q is only ever used in rank-based/relative comparisons in the
   design (never an absolute threshold), so any monotone percentile
   formula is numerically inert to the registered rho *as long as* every
   gene is always scored against the same fixed N=20 controls. Refusing to
   compute Q on a short control set (rather than "rescuing" a gene by
   computing it on fewer) is exactly what keeps N fixed and the formula
   choice inert -- do not change that without revisiting the amendment.
