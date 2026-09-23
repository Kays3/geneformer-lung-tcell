# S100 ISP analysis layer

Built 2026-09-23 (card `s100-isp-execution-20260922`) against
`hive/reports/s100-isp-design-20260922.md` and `retained_rows_spec_20260922.md`.
No GPU. No dependency on `geneformer`/`torch` anywhere in this directory --
every module here is plain `pandas`/`numpy`/`scipy` and is fully testable
without a GPU host, a real dataset, or the vendored Geneformer checkout.

## What's built (fully specified, tested against synthetic fixtures)

| Module | What it does |
|---|---|
| `retained_rows.py` | Builds the planned Cartesian table, left-joined against eligibility manifests, run-completion markers, no-op results, and ISP stats CSVs. As registered, 12 genes x 2 sources x 2 goals x 2 operations = 96 core rows; the real Module A run is now 12 genes x 1 source (LUAD only, SCLC source arm dropped 2026-09-23) x 2 goals x 2 operations = **48** (`MODULE_A_SOURCES = ("luad",)`) -- see `retained_rows_spec_20260922.md`'s second dated amendment. The builder itself stays fully generic over `sources`; its own tests still exercise both source branches. A missing/not-yet-run/not-estimable row is always kept, never dropped, with a machine-readable `status`/`status_reason`. Computes `donor_balanced_shift` directly from the raw per-cell pickle + the paired-eligible manifest's donor column -- **not** from `InSilicoPerturberStats`' own `Shift_to_goal_end` column, which is a plain cell-weighted mean, not the donor-weighted mean the design requires. |
| `ambient_stats.py` | `E[g,c]`, `Q[g,c]` (midrank percentile within matched controls), per-stratum rollup, the primary Spearman test with an **exact two-sided** (not asymptotic) label-permutation p-value enumerating all `n!` permutations, the leave-one-gene-out single-gene-carry check (`primary_test_with_single_gene_check`), the zero-GPU stability gates (120 leave-one-control-out recomputations, 10,000 within-stratum bootstrap redraws at seed `20260922`), and the `0.00072` numerical effect floor. |
| `matched_controls.py` | **Unblocked 2026-09-23.** `build_matched_control_table()` builds the six fixed 20-control matched strata from a LUAD-only per-gene detection-fraction/median-token-rank table (`median_token_rank_and_detection()` / `load_and_freeze_luad_gene_stats()`), matched on both axes against every stratum member, deterministic seeded sampling, `not_estimable_control_stratum` reported rather than rescued below 20 candidates. `synthetic_control_table()` remains for `ambient_stats.py`'s own tests. |
| `provenance_utils.py` | `sha256_file()` -- every file this layer reads gets hashed from its own bytes at the moment it's read, same principle as `run_targeted_panel.py`'s self-hash/input-hash (2026-09-23). |

Run the tests (no GPU, ~5s total, the slow part is the honest n=10
factorial-permutation check):

```
python3 test_retained_rows.py
python3 test_ambient_stats.py
python3 test_matched_controls.py
```

## Mandatory CPU precheck before any GPU work (2026-09-23)

Of the eight LUAD-eligible non-anchor genes, exactly one (`S100A2`) is
ambient-flagged. Their `ambient_risk` values are strongly bimodal (four
near zero, four between 0.74 and 0.96, a gap larger than either cluster's
own range) -- see `retained_rows_spec_20260922.md`'s second dated
amendment for the full table and the reasoning. The primary test's real
`n` is 8, not 10; the exact two-sided permutation math for both n was
independently re-verified there too, correcting an n=8 threshold that had
first been circulated using a one-sided convention -- the registered gate
itself is unchanged, but the supporting arithmetic was wrong and is
corrected in the amendment.

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

## Second round of rulings -- 2026-09-23 (human ruling)

Full text in `retained_rows_spec_20260922.md`'s second dated amendment;
summary here:

1. **SCLC source arm dropped.** No ambient-flagged gene survives
   eligibility in SCLC; a pure two-group difference there would pass the
   primary gate 47% of the time regardless of information content. `SCLC`
   remains a valid goal (`LUAD -> SCLC`). Real run: 48 core rows, not 96.
2. **Median token rank ruled: computed from the tokenized held-out
   dataset**, 0-based convention, frozen and hashed before use.
3. **Matching scope ruled: per-source-state (LUAD-only), not global** --
   both matching axes are source-dependent. `matched_controls.py` is
   unblocked.
4. **n=8, not n=10, for the real primary test** (`S100P`/`S100A16` fail
   LUAD eligibility). **The exact two-sided critical rho at n=8 was
   recomputed and corrected**: ~0.7381, not the one-sided ~0.6190 first
   circulated -- at n=8 (unlike n=10) the p<=0.05 condition is tighter
   than the rho>=0.70 gate over part of its range. The registered gate
   itself does not change; the arithmetic offered in support of it did.
5. **Leave-one-gene-out registered before the number exists.**
   **Corrected same day (Amendment 3):** the first floor (0.60, borrowed
   from leave-one-control-out) was wrong -- n=7 is a different null, not
   the same n as leave-one-control-out. A first replacement (0.75) was
   also wrong: Spearman's rho at n=7 is discrete and 0.75 falls in a gap
   between attainable values (its own exact p is 0.066, which fails). The
   correct three-way rule, evaluated two-sided at the real n=7:
   `rho >= 0.785714 AND p<=0.05` -> "survives"; `0 < rho` but `p>0.05` ->
   "gene_sensitive_open" (no claim either way); `rho <= 0` -> "carried by
   a single gene." `ambient_stats.leave_one_gene_out_check()` /
   `primary_test_with_single_gene_check()`.

CPU precheck (mandatory before GPU, done here, no ts1 needed -- the
ambient-risk table already exists locally): the eight LUAD-eligible
non-anchor genes' `ambient_risk` values are strongly bimodal (four near
zero, four at 0.74-0.96). See the amendment for the full table and why
this is a different situation from the one that dropped SCLC.

**Lesson worth keeping (Amendment 3):** a critical value for a discrete
statistic (small-n Spearman rho) must come from scanning attainable
statistic values for the first whose own exact p clears the threshold --
never from indexing a sorted null at an approximate quantile position,
which can land inside a gap and silently pick the wrong side of a
boundary. Two independent numeric errors in one day (Amendment 2's own
n=8 table, then the leave-one-gene-out floor) both came from that same
shortcut.

## GPU HELD -- Amendment 4 (2026-09-23)

The CPU precheck (Amendment 2) found something bigger than first
described: the 8 LUAD-eligible non-anchor genes split 4-vs-4 into an
ambient-high and ambient-low cluster, and a pure between-cluster Q
difference (zero within-cluster rank information) clears the corrected
primary gate **63.4%** of the time (365/576, verified independently) --
worse than the 47% that got SCLC dropped. GPU is held pending human
review. Two new diagnostics registered and implemented, always reported:
`ambient_stats.exact_group_separation_test()` (exact 4-vs-4 Mann-Whitney,
min p=2/70=0.02857 at complete separation) and `within_cluster_spearman()`
(rho within each cluster separately -- if both are near zero when the real
numbers exist, the result must be reported as a two-group contrast, not a
monotone association; no numeric "near zero" cutoff is fixed yet). One
data-correction: `S100PBP`'s ambient_risk (0.744492, by ensembl_id
ENSG00000116497) is confirmed correct against a proposed replacement that
was actually a different gene, `PEBP1`, matched by a similar-looking
symbol. Full text in the amendment.
