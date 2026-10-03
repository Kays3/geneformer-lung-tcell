"""ISP-BAL-1: executable checks for the balanced Geneformer ISP criteria.

Each function implements one criterion of protocol_criteria_e3/BALANCED_CRITERIA.md
and returns a plain dict (or a DataFrame for gene tables), so a check can be
applied to any run's outputs and its result stored as JSON. The module needs
only numpy, pandas and scipy; it never imports geneformer, so it runs on CPU
against saved outputs.

Criterion IDs (A1 ... E1) match BALANCED_CRITERIA.md. Sources are the
Geneformer protocol (P1-P10) and this repository's criteria (O1-O19), listed in
protocol_criteria_e3/CRITERIA_INVENTORY.md.
"""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 0.05
CLASSIFIER_BA_MIN = 0.60       # B2 (O4)
NOISE_BAND = 0.045             # B3 (O5)
N_CONTROLS_MIN = 20            # D4 (O8)
PRECISION_PEARSON_MIN = 0.99   # C2 (P7)
PRECISION_SPEARMAN_MIN = 0.95  # C2 (O16)
PRECISION_TOP_K = 20           # C2 (O16)
PRECISION_TOP_MIN = 19         # C2 (O16)
HOST_MAX_ABS = 1e-3            # C3 (O15)
HOST_SPEARMAN_MIN = 0.999      # C3 (O15)
SCREEN_FPR_MAX = 0.05          # D2 (new)
NULL_QUANTILE = 0.95           # D7 (O13)
FOLD_FLIPS_MAX = 1             # D8 (O14)


def _result(criterion: str, passed: bool, **details) -> dict:
    return {"criterion": criterion, "pass": bool(passed), **details}


# ----------------------------------------------------------------- A: design

def check_donor_disjoint(splits: pd.DataFrame) -> dict:
    """A1 (P1, O3). Every donor sits in exactly one partition.

    `splits` has columns donor and partition, one row per sample or cell.
    """
    per_donor = splits.groupby("donor")["partition"].nunique()
    leaking = sorted(per_donor[per_donor > 1].index)
    return _result("A1", not leaking, n_donors=int(per_donor.size), leaking_donors=leaking)


def check_class_confounding(meta: pd.DataFrame, source: str = "study") -> dict:
    """A2 (O1). No pair of classes may be separated by source alone.

    A pair of classes is "bridged" if at least one source (study, site, batch)
    contributes cells of both. An unbridged pair means every difference between
    those classes is also a difference between sources. `meta` has columns
    label and `source`.
    """
    classes = sorted(meta["label"].unique())
    by_source = meta.groupby(source)["label"].agg(lambda s: set(s))
    unbridged = [
        (a, b) for a, b in combinations(classes, 2)
        if not any({a, b} <= labels for labels in by_source)
    ]
    return _result("A2", not unbridged, classes=classes, unbridged_pairs=unbridged,
                   n_sources=int(by_source.size))


def check_confounder_balance(meta: pd.DataFrame, covariate: str) -> dict:
    """A3 (P1 `attr_to_balance`). A covariate does not differ across partitions.

    One row per donor, with columns partition and `covariate`. Categorical
    covariates use a chi-square test of independence, numeric ones a
    Kruskal-Wallis test. Pass if p > ALPHA.
    """
    col = meta[covariate]
    if pd.api.types.is_numeric_dtype(col):
        groups = [g[covariate].to_numpy() for _, g in meta.groupby("partition") if len(g)]
        p = float(stats.kruskal(*groups).pvalue) if len(groups) > 1 else 1.0
        test = "kruskal"
    else:
        table = pd.crosstab(meta["partition"], col)
        p = float(stats.chi2_contingency(table).pvalue) if table.shape[0] > 1 and table.shape[1] > 1 else 1.0
        test = "chi2"
    return _result("A3", p > ALPHA, covariate=covariate, test=test, p=p)


def check_equal_donor_weight(cells_per_donor: pd.Series, cap: int) -> dict:
    """A4 (O2). Every donor contributes exactly `cap` analysis cells."""
    off = cells_per_donor[cells_per_donor != cap]
    return _result("A4", off.empty, cap=cap, donors_off_cap=off.to_dict())


# ------------------------------------------------------------- B: classifier

def _balanced_accuracy(y_true, y_pred) -> float:
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    recalls = [np.mean(y_pred[y_true == c] == c) for c in np.unique(y_true)]
    return float(np.mean(recalls))


def _macro_f1(y_true, y_pred, labels) -> float:
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    f1s = []
    for c in labels:
        tp = np.sum((y_pred == c) & (y_true == c))
        fp = np.sum((y_pred == c) & (y_true != c))
        fn = np.sum((y_pred != c) & (y_true == c))
        denom = 2 * tp + fp + fn
        f1s.append(0.0 if denom == 0 else 2 * tp / denom)
    return float(np.mean(f1s))


def classifier_gate(pred: pd.DataFrame, n_boot: int = 2000, seed: int = 20261003) -> dict:
    """B1-B3 (P3, O4, O5). Held-out classification, scored two ways.

    `pred` has one row per held-out cell: donor, y_true, y_pred. All donors must
    have been scored by a model that never saw them (A1).

    - B1 (P3): confusion matrix and macro F1, with a 95% CI from a bootstrap over
      donors. Pass if the lower bound exceeds 1/K, the expected macro F1 of
      uniform guessing over K classes.
    - B2 (O4): pooled balanced accuracy >= 0.60, and an exact two-sided sign
      test over donors (per-donor balanced accuracy above 1/K) with p <= 0.05
      and more donors above chance than below.
    - B3 (O5): `near_threshold` is True when pooled balanced accuracy lies
      within NOISE_BAND of 0.60. It never changes the decision.
    """
    labels = sorted(pred["y_true"].unique())
    k = len(labels)
    chance = 1.0 / k
    cm = pd.crosstab(pred["y_true"], pred["y_pred"]).reindex(index=labels, columns=labels, fill_value=0)
    macro_f1 = _macro_f1(pred["y_true"], pred["y_pred"], labels)
    donors = pred["donor"].unique()
    rng = np.random.default_rng(seed)
    groups = {d: g for d, g in pred.groupby("donor")}
    boots = []
    for _ in range(n_boot):
        sample = pd.concat([groups[d] for d in rng.choice(donors, size=len(donors), replace=True)])
        boots.append(_macro_f1(sample["y_true"], sample["y_pred"], labels))
    ci = (float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975)))
    pooled_ba = _balanced_accuracy(pred["y_true"], pred["y_pred"])
    per_donor = {d: _balanced_accuracy(g["y_true"], g["y_pred"]) for d, g in groups.items()}
    above = sum(v > chance for v in per_donor.values())
    below = sum(v < chance for v in per_donor.values())
    p_sign = float(stats.binomtest(above, above + below, 0.5).pvalue) if above + below else 1.0
    b1 = ci[0] > chance
    b2 = pooled_ba >= CLASSIFIER_BA_MIN and p_sign <= ALPHA and above > below
    return {
        "criterion": "B1-B3",
        "B1_pass": bool(b1), "macro_f1": macro_f1, "macro_f1_ci95": ci, "chance": chance,
        "confusion_matrix": cm.to_dict(),
        "B2_pass": bool(b2), "pooled_balanced_accuracy": pooled_ba,
        "donors_above_chance": int(above), "donors_below_chance": int(below), "sign_test_p": p_sign,
        "per_donor_balanced_accuracy": per_donor,
        "B3_near_threshold": bool(abs(pooled_ba - CLASSIFIER_BA_MIN) <= NOISE_BAND),
        "pass": bool(b1 and b2),
    }


def permutation_calibration(gate_on_permuted_labels: dict) -> dict:
    """B4 (new). A classifier trained on donor-permuted labels must fail the gate."""
    return _result("B4", not gate_on_permuted_labels["pass"])


# ---------------------------------------------------- C: perturbation integrity

def noop_gate(shifts) -> dict:
    """C1 (O6). A perturbation with no edit gives a shift of exactly 0 in every cell."""
    shifts = np.asarray(shifts, dtype=float)
    return _result("C1", bool(np.all(shifts == 0.0)), n_cells=int(shifts.size),
                   max_abs=float(np.max(np.abs(shifts))) if shifts.size else 0.0)


def precision_agreement(reference, other) -> dict:
    """C2 (P7, O16). Per-gene shifts under two precisions (or two model files) agree.

    Pass needs Pearson r >= 0.99 (the protocol's quantization check), Spearman
    rho >= 0.95 and a top-20 overlap of at least 19 by absolute shift (this
    repository's bf16 canary). Inputs are aligned 1-d arrays of per-gene shifts.
    """
    a, b = np.asarray(reference, float), np.asarray(other, float)
    r = float(stats.pearsonr(a, b).statistic)
    rho = float(stats.spearmanr(a, b).statistic)
    k = min(PRECISION_TOP_K, a.size)
    top = len(set(np.argsort(-np.abs(a))[:k]) & set(np.argsort(-np.abs(b))[:k]))
    need = PRECISION_TOP_MIN if k == PRECISION_TOP_K else k - 1
    passed = r >= PRECISION_PEARSON_MIN and rho >= PRECISION_SPEARMAN_MIN and top >= need
    return _result("C2", passed, pearson_r=r, spearman_rho=rho, top_k=k, top_overlap=top)


def cross_host_equivalence(host_a, host_b) -> dict:
    """C3 (O15). The same probe on two hosts: per-cell max |difference| <= 1e-3
    and Spearman rho >= 0.999."""
    a, b = np.asarray(host_a, float), np.asarray(host_b, float)
    max_abs = float(np.max(np.abs(a - b)))
    rho = float(stats.spearmanr(a, b).statistic)
    return _result("C3", max_abs <= HOST_MAX_ABS and rho >= HOST_SPEARMAN_MIN, max_abs=max_abs, spearman_rho=rho)


def config_consistency(state_cfg: dict, isp_cfg: dict, keys=("emb_layer", "emb_mode", "filter_data")) -> dict:
    """C4 (P8). State embeddings and perturbation use the same layer, mode and cell filter."""
    diff = {k: (state_cfg.get(k), isp_cfg.get(k)) for k in keys if state_cfg.get(k) != isp_cfg.get(k)}
    return _result("C4", not diff, mismatched=diff)


# ------------------------------------------------------ multiple testing helpers

def holm(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    m = p.size
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def bh(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p)
    adj = np.empty(m)
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * m / (rank + 1))
        adj[i] = running
    return adj


def attainable_min_p(n_donors: int) -> float:
    """The smallest two-sided p an exact Wilcoxon signed-rank test can give with n donors."""
    return 2.0 / 2 ** n_donors


def max_holm_family(n_donors: int, alpha: float = ALPHA) -> int:
    """D5 (O9). The largest family in which Holm's strictest step is attainable."""
    return int(math.floor(alpha / attainable_min_p(n_donors) + 1e-12))


def family_correction(pvals, n_donors: int, alpha: float = ALPHA) -> tuple[np.ndarray, str]:
    """D5 (O9, O10, balanced). Holm when the family can reach significance with this
    many donors; otherwise BH, labelled as the power adaptation."""
    m = len(pvals)
    if m <= max_holm_family(n_donors, alpha):
        return holm(pvals), "HOLM"
    return bh(pvals), "BH_POWER_ADAPTED"


# ------------------------------------------------- D: gene-level evidence

def screen_goal_state_shift(cell_shifts: dict, random_pool) -> pd.DataFrame:
    """D1 (P4, P6). The protocol's screen: each gene's per-cell shifts toward the
    goal against the shifts from perturbing random genes, two-sided Wilcoxon
    rank-sum, BH across genes, Sig at FDR < 0.05.

    `cell_shifts` maps gene -> 1-d array of per-cell shifts; `random_pool` is
    the pooled per-cell shifts of random genes.
    """
    pool = np.asarray(random_pool, float)
    rows = []
    for gene, x in cell_shifts.items():
        x = np.asarray(x, float)
        p = float(stats.mannwhitneyu(x, pool, alternative="two-sided").pvalue)
        rows.append({"gene": gene, "Shift_to_goal_end": float(np.mean(x)), "N_Detections": int(x.size), "p": p})
    out = pd.DataFrame(rows)
    out["FDR"] = bh(out["p"].to_numpy()) if len(out) else []
    out["Sig"] = (out["FDR"] < ALPHA).astype(int)
    return out


def screen_calibration(null_cell_shifts: dict) -> dict:
    """D2 (new). Calibration of the D1 screen on random null genes.

    Each null gene is screened against the pooled cells of the other null genes,
    and the share with p < 0.05 is the screen's false-positive rate. The screen
    fails calibration only on evidence of excess: when the exact binomial lower
    95% bound of that rate exceeds 0.05 (with 100 null genes, 11 or more hits).
    A point estimate a little above 0.05 is expected by chance and does not fail.
    When it fails, D1 is UNCALIBRATED for this run and may rank genes but never
    promote one.
    """
    genes = list(null_cell_shifts)
    sig = 0
    for g in genes:
        pool = np.concatenate([np.asarray(null_cell_shifts[h], float) for h in genes if h != g])
        p = float(stats.mannwhitneyu(np.asarray(null_cell_shifts[g], float), pool, alternative="two-sided").pvalue)
        sig += p < ALPHA  # one test per gene against its own pool; BH over one test is the raw p
    n = len(genes)
    fpr = sig / n if n else float("nan")
    ci = stats.binomtest(sig, n).proportion_ci(confidence_level=0.95) if n else None
    return _result("D2", n > 0 and ci.low <= SCREEN_FPR_MAX, n_null=n, n_sig=int(sig), fpr=fpr,
                   fpr_ci95=(float(ci.low), float(ci.high)) if ci else None)


def donor_values(cells: pd.DataFrame) -> pd.Series:
    """D3 (O7). Per-donor mean shift from per-cell output (columns donor, shift)."""
    return cells.groupby("donor")["shift"].mean()


def control_adjust(gene_donor: pd.Series, controls_donor: pd.DataFrame) -> tuple[pd.Series | None, str]:
    """D4 (O8). a(g, d) = s(g, d) minus the median over matched controls in donor d.

    `controls_donor` has one column per control gene and one row per donor.
    Fewer than N_CONTROLS_MIN controls gives NOT_ESTIMABLE_CONTROLS.
    """
    if controls_donor.shape[1] < N_CONTROLS_MIN:
        return None, "NOT_ESTIMABLE_CONTROLS"
    donors = gene_donor.index.intersection(controls_donor.index)
    return gene_donor.loc[donors] - controls_donor.loc[donors].median(axis=1), "OK"


def donor_test(adjusted: pd.Series) -> dict:
    """D5 (O10). Exact two-sided Wilcoxon signed-rank test over donors."""
    x = np.asarray(adjusted, float)
    x = x[x != 0]
    if x.size == 0:
        return {"n": 0, "median": 0.0, "p": 1.0}
    p = float(stats.wilcoxon(x, zero_method="wilcox", method="exact").pvalue)
    return {"n": int(x.size), "median": float(np.median(adjusted)), "p": p}


def dose_concordance(del_sig: bool, del_median: float, ovx_sig: bool, ovx_median: float) -> str:
    """D6 (O11). Deletion and overexpression must move in opposite directions."""
    if del_sig and ovx_sig:
        return "COHERENT" if np.sign(del_median) == -np.sign(ovx_median) else "INCOHERENT"
    if del_sig or ovx_sig:
        return "UNRESOLVED"
    return "NONE"


def alt_state_ok(shift_goal: float, shift_alt: float | None) -> bool | None:
    """D6 (P5). With an alternative state, a candidate moves toward the goal and
    not toward the alternative. None when the design has no alternative state."""
    if shift_alt is None:
        return None
    return bool(shift_goal > 0 and shift_alt <= 0)


def beats_null(effect: float, null_effects) -> dict:
    """D7 (O13). The gene's control-adjusted effect exceeds the 95th percentile of
    the absolute effects of random null genes in the same run."""
    q = float(np.quantile(np.abs(np.asarray(null_effects, float)), NULL_QUANTILE))
    return _result("D7", abs(effect) > q, effect=float(effect), null_q95=q)


def null_pair_correlation(null_del, null_ovx, n_perm: int = 10000, seed: int = 20261003) -> dict:
    """D7 (O13). Rank correlation of deletion and overexpression medians across
    random null genes, with a one-sided permutation p for rho < 0. A negative
    rho means opposed effects are the baseline, which is why D6 alone is not
    enough and D7 is required."""
    a, b = np.asarray(null_del, float), np.asarray(null_ovx, float)
    rho = float(stats.spearmanr(a, b).statistic)
    rng = np.random.default_rng(seed)
    perm = np.array([stats.spearmanr(a, rng.permutation(b)).statistic for _ in range(n_perm)])
    p = float((1 + np.sum(perm <= rho)) / (1 + n_perm))
    return {"criterion": "D7_null", "rho": rho, "p_lower": p, "n_null": int(a.size)}


def fold_stable(fold_medians, overall_median: float) -> dict:
    """D8 (O14). A positive whose sign flips in more than one fold is not stable."""
    flips = int(sum(np.sign(m) == -np.sign(overall_median) for m in fold_medians if m != 0))
    return _result("D8", flips <= FOLD_FLIPS_MAX, flips=flips)


# ------------------------------------------------------- E: evidence levels

LEVELS = ("NOT_ESTIMABLE", "OPEN", "SCREEN", "SUPPORTED", "REPLICATED", "VALIDATED")


def evidence_level(row: dict) -> str:
    """E1 (balanced). The evidence level of one gene in one run.

    Keys used: screen_sig, screen_calibrated, estimable, donor_sig, concordance,
    alt_ok, ambient_clear, beats_null, fold_stable, replicated, validated.

    - SCREEN: D1 Sig = 1 in a run whose screen passed D2. A ranking, not a claim.
    - SUPPORTED: donor-level test significant after D5 correction, COHERENT,
      beats the null (D7), fold-stable (D8), ambient-clear, and toward the goal
      and away from the alternative state when one exists. It does not require
      SCREEN; the screen ranks, the donor test decides.
    - REPLICATED: SUPPORTED, in the same direction, in an independent cohort.
    - VALIDATED: REPLICATED or SUPPORTED, and confirmed experimentally (P10).
    """
    if not row.get("estimable", True):
        return "NOT_ESTIMABLE"
    supported = (
        row.get("donor_sig", False)
        and row.get("concordance") == "COHERENT"
        and row.get("beats_null", False)
        and row.get("fold_stable", False)
        and row.get("ambient_clear", True)
        and row.get("alt_ok") in (True, None)
    )
    if supported and row.get("validated", False):
        return "VALIDATED"
    if supported and row.get("replicated", False):
        return "REPLICATED"
    if supported:
        return "SUPPORTED"
    if row.get("screen_sig", False) and row.get("screen_calibrated", False):
        return "SCREEN"
    return "OPEN"


def evaluate_genes(cells: pd.DataFrame, genes, n_donors: int | None = None) -> pd.DataFrame:
    """Apply D1-D8 and E1 to one run.

    `cells` has one row per perturbed cell with columns gene, role (test,
    control or null), stratum, op (delete or overexpress), donor, shift, and
    optionally fold. Controls are matched to a test gene by stratum. `genes` is
    the family of test genes. Returns one row per test gene with every check and
    its evidence level.
    """
    out = []
    null = cells[cells["role"] == "null"]
    calibration = screen_calibration(
        {g: d["shift"].to_numpy() for g, d in null[null["op"] == "delete"].groupby("gene")}
    ) if len(null) else {"pass": False}
    pool = null.loc[null["op"] == "delete", "shift"].to_numpy()
    screen = screen_goal_state_shift(
        {g: cells.loc[(cells["gene"] == g) & (cells["op"] == "delete"), "shift"].to_numpy() for g in genes}, pool
    ).set_index("gene") if len(pool) else None
    per_op = {}
    null_effects = {}
    for op in ("delete", "overexpress"):
        sub = cells[cells["op"] == op]
        dv = sub.groupby(["gene", "donor"])["shift"].mean().unstack("donor")
        controls = {s: dv.loc[sorted(set(sub.loc[(sub["role"] == "control") & (sub["stratum"] == s), "gene"]))]
                    for s in sub["stratum"].dropna().unique()}
        res = {}
        for g in genes:
            stratum = sub.loc[sub["gene"] == g, "stratum"].iloc[0]
            adj, status = control_adjust(dv.loc[g], controls.get(stratum, pd.DataFrame()).T)
            res[g] = (adj, status)
        per_op[op] = res
        null_genes = sorted(set(sub.loc[sub["role"] == "null", "gene"]))
        null_effects[op] = [float(np.median(dv.loc[n].dropna())) for n in null_genes]
    tests = {op: {g: donor_test(adj) if adj is not None else None for g, (adj, _) in per_op[op].items()}
             for op in per_op}
    nd = n_donors or cells["donor"].nunique()
    adjp = {}
    for op in tests:
        names = [g for g in genes if tests[op][g] is not None]
        corrected, method = family_correction([tests[op][g]["p"] for g in names], nd) if names else ([], "NONE")
        adjp[op] = (dict(zip(names, corrected)), method)
    for g in genes:
        row = {"gene": g}
        estimable = all(per_op[op][g][1] == "OK" for op in per_op)
        row["estimable"] = estimable
        if screen is not None:
            row["screen_sig"] = bool(screen.loc[g, "Sig"])
            row["screen_p"] = float(screen.loc[g, "p"])
        row["screen_calibrated"] = bool(calibration["pass"])
        if estimable:
            d, o = tests["delete"][g], tests["overexpress"][g]
            dp, op_ = adjp["delete"][0][g], adjp["overexpress"][0][g]
            row.update(del_median=d["median"], del_p_adj=dp, ovx_median=o["median"], ovx_p_adj=op_,
                       correction=adjp["delete"][1])
            row["concordance"] = dose_concordance(dp <= ALPHA, d["median"], op_ <= ALPHA, o["median"])
            row["donor_sig"] = dp <= ALPHA
            row["beats_null"] = beats_null(d["median"], null_effects["delete"])["pass"] if null_effects["delete"] else False
            if "fold" in cells.columns:
                # O14 works on control-adjusted values: median a over the donors of each fold.
                adj = per_op["delete"][g][0]
                donor_fold = cells.drop_duplicates("donor").set_index("donor")["fold"]
                folds = [float(np.median(adj[adj.index.isin(donor_fold[donor_fold == f].index)]))
                         for f in sorted(donor_fold.unique()) if adj.index.isin(donor_fold[donor_fold == f].index).any()]
                row["fold_stable"] = fold_stable(folds, d["median"])["pass"]
            else:
                row["fold_stable"] = True
        row["level"] = evidence_level(row)
        out.append(row)
    return pd.DataFrame(out)
