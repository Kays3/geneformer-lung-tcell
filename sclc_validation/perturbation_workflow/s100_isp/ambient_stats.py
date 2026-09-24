"""E/Q statistics, stratum rollup, the primary ambient-association test
(exact label-permutation Spearman), and the zero-GPU control-stability
gates, per hive/reports/s100-isp-design-20260922.md ("Controls and
analysis"). Built and testable now against matched_controls.py's
synthetic_control_table() -- the real control table is blocked (see that
module); everything here is control-table-agnostic and works identically
once the real one exists.

No dependency on geneformer/torch/GPU. Pure pandas/numpy/scipy.
"""
from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata

# bf16 canary's max absolute delta vs fp32 (design doc, "Model and
# precision plan"): an arm's donor-balanced mean absolute shift at or
# below this is numerically indeterminate, not zero.
NUMERICAL_FLOOR = 0.00072

SEED = 20260922
BOOTSTRAP_N = 10_000
MIN_COMMON_CONTROLS = 20

# ---------------------------------------------------------------------------
# E: the ambient-test magnitude, per (gene, contrast)
# ---------------------------------------------------------------------------

_KEY_COLS = ["gene", "ensembl_id", "role", "stratum", "source_state", "goal_state", "alt_state"]


def compute_E(retained_rows: pd.DataFrame) -> pd.DataFrame:
    """E[g,c] = (abs(effect_delete[g,c]) + abs(effect_overexpress[g,c])) / 2.

    NaN whenever either arm's donor_balanced_shift is missing (not
    estimable, not yet run, or failed) -- skipna=False deliberately: E must
    never silently degrade to "just the one arm that happened to complete".
    """
    delete_df = (
        retained_rows[retained_rows["perturbation_type"] == "delete"]
        [_KEY_COLS + ["donor_balanced_shift"]]
        .rename(columns={"donor_balanced_shift": "effect_delete"})
    )
    oe_df = (
        retained_rows[retained_rows["perturbation_type"] == "overexpress"]
        [_KEY_COLS + ["donor_balanced_shift"]]
        .rename(columns={"donor_balanced_shift": "effect_overexpress"})
    )
    merged = pd.merge(delete_df, oe_df, on=_KEY_COLS, how="outer")
    merged["E"] = merged[["effect_delete", "effect_overexpress"]].abs().mean(axis=1, skipna=False)
    return merged


# ---------------------------------------------------------------------------
# Q: midrank percentile of E within a stratum's matched controls
# ---------------------------------------------------------------------------


def midrank_percentile(value: float, controls: list[float]) -> float:
    """Percentile rank of `value` among `controls`, ties assigned their
    average rank (design doc's exact tie-handling rule). The design does
    NOT specify the percentile formula itself, only the tie rule -- this
    uses the standard "value's rank among N+1 combined points" convention:
    combine `value` with the `controls` (N+1 points), rank with averaged
    ties, map the value's own rank r to (r - 0.5) / (N + 1) * 100.

    RULED 2026-09-23 (see retained_rows_spec_20260922.md's dated amendment):
    this formula stands as-is, but not because it's the uniquely correct
    one -- because Q is used only in rank-based/relative contexts
    everywhere it appears in the design (Spearman rho, its permutation
    test, LOO/bootstrap, and the median-Q-by-stratum comparison), never
    against an absolute threshold, and every gene is scored against the
    SAME fixed N=20 controls. Any monotone rank-to-percentile formula
    applied identically across genes preserves both the cross-gene
    ordering and the median comparison, so the registered rho is
    numerically unchanged by which one is used -- a free choice only stops
    being free if N ever varies gene-to-gene (see compute_Q_for_contrast's
    docstring and the amendment for why that must never happen).
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    clean_controls = [c for c in controls if c is not None and not (isinstance(c, float) and np.isnan(c))]
    if len(clean_controls) == 0:
        return np.nan
    combined = np.asarray(clean_controls + [value], dtype=float)
    ranks = rankdata(combined, method="average")
    n = len(combined)
    return 100.0 * (ranks[-1] - 0.5) / n


def compute_Q_for_contrast(
    E_by_gene: dict[str, float],
    stratum_by_gene: dict[str, str],
    control_E_by_stratum: dict[str, dict[str, float]],
    min_common_controls: int = MIN_COMMON_CONTROLS,
) -> dict[str, float]:
    """Q for every gene in E_by_gene, for one fixed contrast. A stratum
    with fewer than `min_common_controls` valid control E-values yields
    NaN for every gene in it (never a Q computed on a short control set) --
    per the design's "not widened, split, or given private controls" rule.

    RULED 2026-09-23 (retained_rows_spec_20260922.md's dated amendment):
    this all-or-nothing behavior is load-bearing beyond the rule it was
    written for. midrank_percentile()'s exact formula is provably inert to
    the registered rho ONLY because every gene is always scored against
    the same N=20 -- do not change this function to "rescue" a gene by
    computing Q on however many controls survived a short stratum; that
    would make N vary gene-to-gene and the percentile-formula choice would
    stop being inert."""
    Q: dict[str, float] = {}
    for gene, E in E_by_gene.items():
        stratum = stratum_by_gene[gene]
        controls = control_E_by_stratum.get(stratum, {})
        values = [v for v in controls.values() if v is not None and not (isinstance(v, float) and np.isnan(v))]
        Q[gene] = midrank_percentile(E, values) if len(values) >= min_common_controls else np.nan
    return Q


def stratum_summary(Q_by_gene: dict[str, float], stratum_by_gene: dict[str, str]) -> pd.DataFrame:
    """Stratum rollup: per-stratum median/mean Q among its member genes,
    for the outcome table's "flagged genes have a higher median Q than
    clean controls" style comparisons. NaN Q values are excluded from the
    median/mean but n_estimable / n_total are both reported so a stratum
    that is mostly not-estimable is visible, not averaged away."""
    rows = []
    for stratum in sorted(set(stratum_by_gene.values())):
        genes = [g for g, s in stratum_by_gene.items() if s == stratum]
        q_values = [Q_by_gene[g] for g in genes]
        estimable = [q for q in q_values if not (isinstance(q, float) and np.isnan(q))]
        rows.append({
            "stratum": stratum,
            "n_genes": len(genes),
            "n_estimable": len(estimable),
            "median_Q": float(np.median(estimable)) if estimable else np.nan,
            "mean_Q": float(np.mean(estimable)) if estimable else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Primary test: Spearman(Q, ambient risk), exact label-permutation p
# ---------------------------------------------------------------------------


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Plain Spearman rho (rank-Pearson), no p-value. Used inside the
    stability loops (LOO, bootstrap), which the design asks to recompute
    "rho" only -- NOT the expensive exact permutation p-value, which would
    be 3,628,800x more work per call if repeated 10,120 times."""
    rank_x = rankdata(x, method="average")
    rank_y = rankdata(y, method="average")
    if rank_x.std() == 0 or rank_y.std() == 0:
        return np.nan
    return float(np.corrcoef(rank_x, rank_y)[0, 1])


def spearman_exact_permutation(Q_values: np.ndarray, risk_values: np.ndarray) -> tuple[float, float, np.ndarray]:
    """Exact two-sided permutation p-value for Spearman's rho, enumerating
    ALL n! label permutations (design doc: "exact two-sided permutation
    p-value uses all label permutations"; n=10 -> 3,628,800 -- must NOT use
    scipy.stats.spearmanr's asymptotic p at this n).

    Vectorized rather than a 3.6-million-iteration Python loop: ranks are
    computed once (ties via average rank), and since every permutation of
    a fixed vector has the same norm, rho for a given permutation reduces
    to a single dot product of centered rank vectors. All n! permutations'
    dot products are computed at once via one matrix multiply.

    Memory note: for n=10 this holds two (3,628,800 x 10) arrays in memory
    at once (~580MB combined) for a few seconds. Fine for n=10; would need
    chunking for materially larger n, which this design never uses.
    """
    n = len(Q_values)
    assert len(risk_values) == n, "Q_values and risk_values must be the same length"

    rank_Q = rankdata(Q_values, method="average")
    rank_risk = rankdata(risk_values, method="average")
    centered_Q = rank_Q - rank_Q.mean()
    centered_risk = rank_risk - rank_risk.mean()
    norm_Q = np.linalg.norm(centered_Q)
    norm_risk = np.linalg.norm(centered_risk)

    observed_rho = float(np.dot(centered_Q, centered_risk) / (norm_Q * norm_risk))

    perm_indices = np.array(list(itertools.permutations(range(n))))  # (n!, n)
    permuted_centered_risk = centered_risk[perm_indices]  # (n!, n)
    dots = permuted_centered_risk @ centered_Q  # (n!,)
    perm_rhos = dots / (norm_Q * norm_risk)

    p_exact = float(np.mean(np.abs(perm_rhos) >= np.abs(observed_rho) - 1e-9))
    return observed_rho, p_exact, perm_rhos


def primary_test(
    Q_by_gene: dict[str, float],
    risk_by_gene: dict[str, float],
    non_anchor_genes: list[str],
    rho_gate: float = 0.70,
    p_gate: float = 0.05,
) -> dict:
    """The primary ambient-association test: Spearman(Q, risk) across the
    non-anchor genes (ten as originally registered; RULED 2026-09-23 --
    eight for the real LUAD-only run, see retained_rows_spec_20260922.md's
    second dated amendment), exact TWO-SIDED permutation p (design doc,
    "Controls and analysis": "exact two-sided permutation p-value... the
    gate is rho >= 0.70 and p <= 0.05"). Gate: rho >= 0.70 AND p <= 0.05 --
    BOTH conditions, independently checked, not "p follows automatically
    once rho clears 0.70". At n=10 the minimum rho for exact two-sided
    p <= 0.05 is ~0.6485, comfortably below 0.70, so rho>=0.70 alone was
    the binding constraint there. At n=8 that is NOT true: the minimum rho
    for exact two-sided p <= 0.05 is ~0.7381 -- HIGHER than the 0.70
    rho-gate. An observed rho in [0.70, 0.7381) at n=8 passes the rho-gate
    and still fails the combined gate on p alone. (Verified 2026-09-23 by
    direct enumeration of the n=8 and n=10 null distributions -- this
    corrects an earlier one-sided-looking estimate of the n=8/n=10 critical
    rho that was circulated before this function existed; see
    retained_rows_spec_20260922.md's second dated amendment for the full
    correction and the exact table.) Report the exact p for whatever n was
    actually used -- never quote an n=10 p-value for an n=8 test. Returns
    gate_pass=False (not an error) if any non-anchor gene's Q is not
    estimable -- the test needs every gene in `non_anchor_genes` estimable,
    whatever its length."""
    Q_vals = np.array([Q_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
    risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
    if np.isnan(Q_vals).any() or np.isnan(risk_vals).any():
        return {"rho": np.nan, "p_exact": np.nan, "gate_pass": False, "n": len(non_anchor_genes),
                "reason": "one or more non-anchor genes has no estimable Q or risk score"}
    rho, p_exact, _ = spearman_exact_permutation(Q_vals, risk_vals)
    return {"rho": rho, "p_exact": p_exact, "n": len(non_anchor_genes),
            "gate_pass": bool(rho >= rho_gate and p_exact <= p_gate), "reason": None}


# n=7's exact attainable rho values are quantized (d^2 is always an even
# integer): near the boundary they are ..., 0.6786, 0.7143, 0.7500, 0.7857,
# ... and nothing between. 0.7500 is NOT the critical value -- its exact
# two-sided p is 0.0663 (fails p<=0.05); the correct smallest attainable
# rho that clears the gate is 0.7857 (p=0.0480). RULED 2026-09-23 (human
# ruling, second correction of the day, s100-isp-execution-20260922):
# Stanley's own first replacement (0.75) landed in the same gap for the
# same reason the original 0.60 borrow was wrong -- a critical value must
# come from scanning ATTAINABLE statistic values for the first whose exact
# p clears the threshold, never from indexing a sorted null at a quantile
# position. Kept here only as a documented, independently-verified
# reference point (by direct enumeration of all 7! = 5040 permutations,
# cross-checked against the human's Fraction-exact derivation) --
# leave_one_gene_out_check() below does not hardcode it; it calls
# spearman_exact_permutation() at whatever n the leave-one-out set actually
# has, so the classification is correct even if a future change alters
# which/how many genes are non-anchor.
LOO_GENE_N7_REFERENCE_RHO = 11 / 14  # == 0.785714285714... (d^2 = 12, n=7: rho = 1 - 6*12/336)


def leave_one_gene_out_check(
    Q_by_gene: dict[str, float],
    risk_by_gene: dict[str, float],
    non_anchor_genes: list[str],
    watch_gene: str,
    p_gate: float = 0.05,
) -> dict:
    """The leave-one-gene-out amendment, registered 2026-09-23 (human
    ruling, s100-isp-execution-20260922) before the real number existed,
    THEN CORRECTED THE SAME DAY (see retained_rows_spec_20260922.md's third
    dated amendment): of the eight LUAD-eligible non-anchor genes, exactly
    one -- S100A2 -- is ambient-flagged; if the other seven cluster
    together in ambient risk, the primary rho is decided by where that one
    gene lands, "a single comparison wearing the clothes of a rank
    correlation." Recomputes the FULL exact two-sided permutation test
    (rho AND p, "evaluated two-sided like the primary" -- not the cheap
    rho-only spearman_rho() the LOO-control-out/bootstrap loops use, since
    those loops run thousands of times and this runs once) with
    `watch_gene` removed. Reported ALWAYS, not only on request.

    Three-way outcome (matching the design's existing control-draw
    "sensitive / open" vocabulary rather than inventing a new one):
      - "survives":     leave-one-out rho > 0 AND exact two-sided p <= p_gate
                         -> the ambient-risk association is NOT solely
                         attributable to watch_gene; a panel-wide statement
                         is permitted (still subject to the full test's own
                         gate_pass).
      - "gene_sensitive_open": leave-one-out rho > 0 but p > p_gate ->
                         no panel-wide claim AND no denial.
      - "carried_by_single_gene": leave-one-out rho <= 0 -> the association
                         does not survive watch_gene's removal at all; MUST
                         NOT be stated as a panel-wide ambient-risk finding,
                         only as a watch_gene-specific one.
      - "not_estimable": a reduced-set gene has no estimable Q or risk.
    `carried_by_single_gene` (bool) mirrors the "carried_by_single_gene"
    outcome for callers that only need the one flag.
    """
    reduced_genes = [g for g in non_anchor_genes if g != watch_gene]
    Q_vals = np.array([Q_by_gene.get(g, np.nan) for g in reduced_genes], dtype=float)
    risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in reduced_genes], dtype=float)

    if np.isnan(Q_vals).any() or np.isnan(risk_vals).any():
        return {"watch_gene": watch_gene, "leave_one_out_n": len(reduced_genes),
                "leave_one_out_rho": np.nan, "leave_one_out_p_exact": np.nan,
                "outcome": "not_estimable", "carried_by_single_gene": False}

    rho, p_exact, _ = spearman_exact_permutation(Q_vals, risk_vals)
    if rho <= 0:
        outcome = "carried_by_single_gene"
    elif p_exact <= p_gate:
        outcome = "survives"
    else:
        outcome = "gene_sensitive_open"
    return {
        "watch_gene": watch_gene, "leave_one_out_n": len(reduced_genes),
        "leave_one_out_rho": rho, "leave_one_out_p_exact": p_exact,
        "outcome": outcome, "carried_by_single_gene": outcome == "carried_by_single_gene",
    }


def primary_test_with_single_gene_check(
    Q_by_gene: dict[str, float],
    risk_by_gene: dict[str, float],
    non_anchor_genes: list[str],
    watch_gene: str,
    rho_gate: float = 0.70,
    p_gate: float = 0.05,
) -> dict:
    """primary_test() plus leave_one_gene_out_check(), merged into one
    result -- convenience wrapper for callers that want both in one call.
    See leave_one_gene_out_check()'s docstring for the three-way outcome.
    """
    full = primary_test(Q_by_gene, risk_by_gene, non_anchor_genes, rho_gate=rho_gate, p_gate=p_gate)
    loo = leave_one_gene_out_check(Q_by_gene, risk_by_gene, non_anchor_genes, watch_gene, p_gate=p_gate)
    return {**full, **loo}


# ---------------------------------------------------------------------------
# 4-vs-4 group-separation diagnostics (registered 2026-09-23, human ruling,
# fourth dated amendment): the 8 LUAD-eligible non-anchor genes' ambient
# risk splits cleanly into a 4-gene high cluster and a 4-gene low cluster.
# A pure between-cluster Q difference, with within-cluster order carrying
# NO rank information at all, clears the corrected primary gate
# (rho >= 0.7381 AND p <= 0.05) 63.4% of the time under random
# within-cluster ordering (365/576 arrangements -- verified independently
# by direct enumeration, matching the human's by-hand figure exactly) --
# worse than the 47% that got the SCLC source arm dropped. The primary
# test's 8-point rank correlation has, in the worst case, roughly one
# degree of freedom (which cluster a gene is in), not eight. These two
# diagnostics measure directly whether the real data has more than that.
# ---------------------------------------------------------------------------


def exact_group_separation_test(
    Q_by_gene: dict[str, float],
    high_ambient_genes: list[str],
    low_ambient_genes: list[str],
) -> dict:
    """Exact two-sided Mann-Whitney U test of Q between the ambient-high
    and ambient-low gene groups -- the test the data's actual 4-vs-4
    structure supports ("is Q higher in the high-ambient group than the
    low-ambient group"), as opposed to the primary test's implicit claim
    of a monotone association across all 8 points. With 4-vs-4 groups and
    no ties, the minimum attainable two-sided p (complete separation) is
    2 / C(8,4) = 2/70 = 0.02857 -- so perfect separation is honestly
    significant by this test, which is exactly the point: this measures
    what the primary rho conflates. Uses scipy's own exact enumeration
    (method="exact"), cross-checked directly against 2/70 for the
    complete-separation case. Reported always."""
    Q_high = np.array([Q_by_gene.get(g, np.nan) for g in high_ambient_genes], dtype=float)
    Q_low = np.array([Q_by_gene.get(g, np.nan) for g in low_ambient_genes], dtype=float)
    if np.isnan(Q_high).any() or np.isnan(Q_low).any():
        return {"statistic": np.nan, "p_exact": np.nan,
                "reason": "one or more genes in the two groups has no estimable Q"}
    result = mannwhitneyu(Q_high, Q_low, alternative="two-sided", method="exact")
    return {"statistic": float(result.statistic), "p_exact": float(result.pvalue), "reason": None}



# RULED 2026-09-23 (human ruling, fifth dated amendment,
# s100-isp-execution-20260922): at n=4, the exact null distribution has
# only 4! = 24 permutations and its minimum attainable two-sided p (at
# rho=+/-1.0, PERFECT correlation) is 2/24 = 0.08333 -- verified directly
# by enumeration below. No within-cluster rho, however extreme, can ever
# clear p <= 0.05 at this n. A numeric "near zero" interpretation
# threshold on this quantity would therefore be theatre, not rigour: it
# cannot rescue a result (nothing here can ever be significant) and it
# cannot condemn one either (a rho of exactly 0 is not more "real" than a
# rho of 0.8 at this n -- both merely fail to reach a floor nothing can
# reach). within_cluster_spearman() is DESCRIPTIVE ONLY for this reason --
# it must never gate an interpretation, only illustrate the shape of the
# data alongside this fixed disclaimer.
WITHIN_CLUSTER_N4_MIN_ATTAINABLE_P = 2 / 24  # == 0.08333...


def within_cluster_spearman(
    Q_by_gene: dict[str, float],
    risk_by_gene: dict[str, float],
    gene_group: dict[str, str],
) -> dict[str, dict]:
    """Spearman rho (AND its own exact two-sided p, cheap at n=4 -- 24
    permutations) computed SEPARATELY within each group in `gene_group`
    (e.g. gene -> "high"/"low" cluster membership) -- measures whether Q
    carries any rank information beyond which cluster a gene is in.

    DESCRIPTIVE ONLY -- never a gate. At the real group size (n=4), the
    minimum attainable two-sided p is `WITHIN_CLUSTER_N4_MIN_ATTAINABLE_P`
    (0.0833), reached only by a PERFECT correlation; no observed value can
    ever be distinguished from chance in either direction at this n. The
    claim that survives or fails is decided by exact_group_separation_test()
    (n=8-vs-C(8,4)=70, genuinely reachable at p<=0.05), never by this
    function -- these rhos exist so a reader can see the shape of the data,
    not to rescue or condemn the group-separation result.

    Reported always, one entry per group present in `gene_group`, each a
    dict with `rho`, `p_exact`, `n`, and the fixed `min_attainable_p_note`
    string (present even when n != 4, computed for whatever n the group
    actually has, so this stays correct if group sizes ever change).
    NaN rho/p for a group with a missing Q/risk value or fewer than 2
    estimable members (a single point has no rank correlation to report).
    """
    result: dict[str, dict] = {}
    for group in sorted(set(gene_group.values())):
        genes = [g for g, grp in gene_group.items() if grp == group]
        n = len(genes)
        Q_vals = np.array([Q_by_gene.get(g, np.nan) for g in genes], dtype=float)
        risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in genes], dtype=float)
        if n < 2 or np.isnan(Q_vals).any() or np.isnan(risk_vals).any():
            result[group] = {"rho": np.nan, "p_exact": np.nan, "n": n,
                              "min_attainable_p_note": None}
            continue
        rho, p_exact, _ = spearman_exact_permutation(Q_vals, risk_vals)
        min_attainable_p = 2.0 / math.factorial(n)
        result[group] = {
            "rho": rho, "p_exact": p_exact, "n": n,
            "min_attainable_p_note": (
                f"descriptive only -- at n={n}, the minimum attainable two-sided "
                f"p is {min_attainable_p:.5g} (reached only by a perfect "
                f"correlation); this rho cannot be distinguished from chance "
                f"unless {min_attainable_p:.5g} <= 0.05, and cannot gate the "
                f"result either way -- see exact_group_separation_test() for "
                f"the test that actually can reach significance."
            ),
        }
    return result


# ---------------------------------------------------------------------------
# Stability: leave-one-control-out + within-stratum bootstrap
# ---------------------------------------------------------------------------


def leave_one_control_out(
    E_by_gene: dict[str, float],
    stratum_by_gene: dict[str, str],
    control_E_by_stratum: dict[str, dict[str, float]],
    risk_by_gene: dict[str, float],
    non_anchor_genes: list[str],
    min_common_controls: int = MIN_COMMON_CONTROLS,
) -> pd.DataFrame:
    """All 120 stratum-control-pair leave-one-out recomputations (design
    doc: "for each of the 120 stratum-control pairs, omit that one
    control, recompute all affected Q values and the ten-gene primary
    rho"). Returns one row per (stratum, dropped_control) with the
    resulting rho -- rho only, not the exact p (see spearman_rho's
    docstring)."""
    rows = []
    for stratum, controls in control_E_by_stratum.items():
        for dropped_id in controls:
            reduced = {cid: v for cid, v in controls.items() if cid != dropped_id}
            adjusted = dict(control_E_by_stratum)
            adjusted[stratum] = reduced
            Q_adj = compute_Q_for_contrast(E_by_gene, stratum_by_gene, adjusted, min_common_controls)
            Q_vals = np.array([Q_adj.get(g, np.nan) for g in non_anchor_genes], dtype=float)
            risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
            rho = np.nan if (np.isnan(Q_vals).any() or np.isnan(risk_vals).any()) else spearman_rho(Q_vals, risk_vals)
            rows.append({"stratum": stratum, "dropped_control": dropped_id, "rho": rho})
    return pd.DataFrame(rows)


def loo_gate_passes(loo_df: pd.DataFrame, min_rho: float = 0.60) -> bool:
    """Design doc: "every leave-one-control-out rho remains positive and
    at least 0.60". A NaN rho (e.g. a stratum that became not-estimable
    after dropping its one control) fails the gate -- it is not skipped."""
    if loo_df["rho"].isna().any():
        return False
    return bool((loo_df["rho"] > 0).all() and (loo_df["rho"] >= min_rho).all())


def bootstrap_stability(
    E_by_gene: dict[str, float],
    stratum_by_gene: dict[str, str],
    control_E_by_stratum: dict[str, dict[str, float]],
    risk_by_gene: dict[str, float],
    non_anchor_genes: list[str],
    seed: int = SEED,
    n_boot: int = BOOTSTRAP_N,
    min_common_controls: int = MIN_COMMON_CONTROLS,
) -> np.ndarray:
    """10,000 within-stratum bootstrap redraws (design doc: "10,000
    bootstrap redraws of 20 controls with replacement within every
    stratum simultaneously"), seed 20260922. Every stratum is resampled on
    every draw (simultaneously, not independently per stratum per draw
    order) using ONE shared Generator so the sequence is fully
    reproducible from the seed alone. Returns the array of 10,000 rhos
    (rho only, not the exact p)."""
    rng = np.random.default_rng(seed)
    strata = list(control_E_by_stratum.keys())
    control_value_lists = {s: list(control_E_by_stratum[s].values()) for s in strata}
    rhos = np.empty(n_boot, dtype=float)

    for i in range(n_boot):
        resampled: dict[str, dict[str, float]] = {}
        for s in strata:
            vals = control_value_lists[s]
            draw = rng.choice(vals, size=len(vals), replace=True)
            resampled[s] = {f"boot_{j}": v for j, v in enumerate(draw)}
        Q_boot = compute_Q_for_contrast(E_by_gene, stratum_by_gene, resampled, min_common_controls)
        Q_vals = np.array([Q_boot.get(g, np.nan) for g in non_anchor_genes], dtype=float)
        risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
        rhos[i] = np.nan if (np.isnan(Q_vals).any() or np.isnan(risk_vals).any()) else spearman_rho(Q_vals, risk_vals)

    return rhos


def bootstrap_gate_passes(rhos: np.ndarray, min_rho: float = 0.70, min_fraction: float = 0.95) -> bool:
    """Design doc: "at least 95% of bootstrap rhos are >= 0.70". A NaN rho
    (control set became not-estimable on that draw) counts against the
    gate, not excluded from the denominator -- it is not "no information",
    it is "this draw did not clear 0.70"."""
    return bool(np.mean(np.nan_to_num(rhos, nan=-np.inf) >= min_rho) >= min_fraction)


# ---------------------------------------------------------------------------
# Numerical effect floor
# ---------------------------------------------------------------------------


def above_numerical_floor(donor_balanced_shift) -> bool | None:
    """None (not True/False) when the shift itself is missing -- 'not
    estimable' must never silently read as 'below floor'."""
    if donor_balanced_shift is None or (isinstance(donor_balanced_shift, float) and np.isnan(donor_balanced_shift)):
        return None
    return bool(abs(donor_balanced_shift) > NUMERICAL_FLOOR)
