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

import numpy as np
import pandas as pd
from scipy.stats import rankdata

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

    FLAGGED, not silently assumed: confirm this formula with Pam/Michael.
    It is isolated to this one function -- changing the convention later
    touches nothing else in this module.
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
    NaN for every gene in it -- not_estimable_control_stratum, per the
    design's "not widened, split, or given private controls" rule."""
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
    ten non-anchor genes, exact permutation p. Gate: rho >= 0.70 and
    p <= 0.05. Returns gate_pass=False (not an error) if any non-anchor
    gene's Q is not estimable -- the ten-gene test needs all ten."""
    Q_vals = np.array([Q_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
    risk_vals = np.array([risk_by_gene.get(g, np.nan) for g in non_anchor_genes], dtype=float)
    if np.isnan(Q_vals).any() or np.isnan(risk_vals).any():
        return {"rho": np.nan, "p_exact": np.nan, "gate_pass": False,
                "reason": "one or more non-anchor genes has no estimable Q or risk score"}
    rho, p_exact, _ = spearman_exact_permutation(Q_vals, risk_vals)
    return {"rho": rho, "p_exact": p_exact, "gate_pass": bool(rho >= rho_gate and p_exact <= p_gate), "reason": None}


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
