#!/usr/bin/env python3
"""Why does T5a's complete-population curated-exhaustion score (RESULTS_T5.md
Goal 2) disagree with T2's test-only score on SCLC vs LUAD ordering? Both read
correctly off the same gene-symbol-mapped pipeline (see differential_expression.py's
`gene_symbol` column) -- this checks whether the disagreement is donor composition,
per PLAN.md's own instruction to report T5a per-donor so a single outlier donor in
the larger pool stays visible.

Runs on the laptop; needs pseudobulk_per_donor_{complete,test_only}.csv, already
pulled back from the compute host.

Reproducibility note (2026-09-19): only pseudobulk_per_donor_test_only.csv
(~10 MB) is committed. pseudobulk_per_donor_complete.csv is ~55 MB -- 6x the
largest file already tracked in this repo, and this repo has no git-lfs -- so
it is regenerated on demand instead of committed, same as every other
compute-host artifact this script's sibling scripts depend on (the H5AD
itself is never committed either). See the missing-file error below for the
exact regeneration command.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
EXHAUSTION = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]


def donor_scores(population: str) -> pd.DataFrame:
    pb_path = RESULTS / f"pseudobulk_per_donor_{population}.csv"
    if not pb_path.exists():
        raise SystemExit(
            f"Missing {pb_path}.\n\n"
            f"pseudobulk_per_donor_{population}.csv is a compute-host artifact "
            "and (for the 'complete' population) too large to commit. Regenerate it:\n\n"
            f"  # on the compute host, needs the prepared H5AD:\n"
            f"  POPULATION={population} python3 differential_expression.py\n\n"
            f"then copy sclc_validation/immune_axis_test/results/pseudobulk_per_donor_{population}.csv "
            "back into this results/ directory before re-running this script. "
            "See README.md's 'What runs today' T5 section for the full command pair."
        )
    pb = pd.read_csv(pb_path)
    ex = pb[pb["gene_symbol"].isin(EXHAUSTION)]
    scored = ex.groupby(["state", "donor"]).agg(
        mean_log1p_cp10k=("mean_log1p_cp10k", "mean"), n_cells=("n_cells", "first")
    ).reset_index()
    scored["population"] = population
    scored["cell_share_in_state"] = scored["n_cells"] / scored.groupby("state")["n_cells"].transform("sum")
    return scored


def main() -> None:
    complete = donor_scores("complete")
    test_only = donor_scores("test_only")
    test_donors = set(test_only["donor"])
    complete["in_test_split"] = complete["donor"].isin(test_donors)

    out = complete.sort_values(["state", "mean_log1p_cp10k"], ascending=[True, False])
    out_path = RESULTS / "t5a_donor_composition_check.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows)\n")

    for state in ("sclc", "luad", "normal"):
        sub = complete[complete["state"] == state]
        print(f"--- {state}: {len(sub)} donors in the complete population ---")
        print(sub[["donor", "mean_log1p_cp10k", "n_cells", "in_test_split"]]
              .sort_values("mean_log1p_cp10k", ascending=False).to_string(index=False))
        t = test_only[test_only["state"] == state].sort_values("mean_log1p_cp10k", ascending=False)
        print(f"\n  test-split donors and their share of test-split {state} cells:")
        print(t[["donor", "mean_log1p_cp10k", "n_cells", "cell_share_in_state"]].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
