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

## Known gap surfaced while building this (not resolved unilaterally)

`retained_rows_spec_20260922.md`'s `status` enum has no value for "eligible,
but the run has not happened yet" -- every one of its six values describes
a *completed* outcome. `retained_rows.py` adds `STATUS_NOT_RUN = "not_run"`,
outside that enum, to describe the state this table is necessarily built in
right now (before Module A has run). If the spec's intent is that this
table is only ever built strictly after a complete run, `STATUS_NOT_RUN`
rows simply never occur in that case and this is inert, not wrong -- but it
was not in the spec, so it's flagged rather than assumed. See
`retained_rows.py`'s module docstring for the full reasoning.

Also flagged, not resolved: `midrank_percentile()`'s exact percentile
formula (the design specifies the tie-handling rule but not the formula
itself), and whether a completed row's `status` should be overwritten to
`not_estimable_control_stratum` when its stratum's control set fails the
20-control gate, versus leaving `status` describing only the gene's own run
outcome and letting `matched_control_percentile_q` stay null. This module
implements the latter (`status` is about the gene's own run; the stub's
absence does not retroactively unmark rows that already completed).
