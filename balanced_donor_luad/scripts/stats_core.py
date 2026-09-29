"""Registered statistics for the balanced-donor ISP (PHASE5_ISP_REGISTRATION.md s.5, Amendment 3).

Everything here is exact and CPU-only:
  - exact two-sided Wilcoxon signed-rank p over donors, from the full null
    distribution built by dynamic programming over integer (doubled) ranks,
    with zeros dropped and ties given midranks. Exact rational arithmetic
    (Fraction), no normal approximation;
  - an exact distribution-free CI for the median (Walsh averages, critical
    value found by SCANNING attainable values, never by indexing a sorted null);
  - Holm step-down adjustment;
  - dose concordance and the registered outcome-row status tables, including
    the mandatory `claim_qualifier` (Amendment 3, 3.8) and `design_warrant`
    (3.6) fields, and a validator that refuses rows missing them.
"""
from fractions import Fraction
from itertools import combinations_with_replacement

import numpy as np

PANEL_A_QUALIFIER = "non-replication on 316M; model-vs-design attribution not established"
WARRANT_PAIRED = "within-donor paired"
WARRANT_LUSC = "between-donor (study-overlap guarded)"
NON_REPLICATION = {"OPEN", "DELETION_ONLY", "DOSE_INCOHERENT", "REVERSED"}


def _doubled_midranks(absvals):
    """Midranks of |x|, doubled so they are integers (ties -> shared midrank)."""
    order = np.argsort(absvals, kind="mergesort")
    ranks2 = np.empty(len(absvals), dtype=np.int64)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and absvals[order[j + 1]] == absvals[order[i]]:
            j += 1
        # positions i..j (0-based) share ranks i+1..j+1; doubled midrank = i+1 + j+1
        for k in range(i, j + 1):
            ranks2[order[k]] = (i + 1) + (j + 1)
        i = j + 1
    return ranks2


def signed_rank_null(ranks2):
    """Exact null distribution of W2 = sum of doubled ranks with positive sign.

    Returns (support, counts) with sum(counts) == 2**n.
    """
    total = int(sum(ranks2))
    counts = np.zeros(total + 1, dtype=object)
    counts[0] = 1
    for r in ranks2:
        r = int(r)
        new = counts.copy()
        new[r:] = new[r:] + counts[: total + 1 - r]
        counts = new
    return np.arange(total + 1), counts


def wilcoxon_exact(x):
    """Exact two-sided Wilcoxon signed-rank test.

    Returns dict(n_used, w2, p) where p is a Fraction:
    P(|W2 - mean| >= |w2_obs - mean|) under the exact sign-flip null.
    Zeros are dropped (registered). n_used == 0 -> p = 1.
    """
    x = np.asarray(x, dtype=float)
    x = x[x != 0]
    n = len(x)
    if n == 0:
        return {"n_used": 0, "w2": 0, "p": Fraction(1)}
    ranks2 = _doubled_midranks(np.abs(x))
    w2 = int(ranks2[x > 0].sum())
    support, counts = signed_rank_null(ranks2)
    total2 = int(ranks2.sum())            # mean of W2 is total2/2
    dev_obs = abs(2 * w2 - total2)        # compare 2*|W2-mean| as integers
    extreme = sum(int(c) for s, c in zip(support, counts) if c and abs(2 * int(s) - total2) >= dev_obs)
    return {"n_used": n, "w2": w2, "p": Fraction(extreme, 2 ** n)}


def min_attainable_p(n):
    """Smallest attainable exact two-sided p with n non-zero, untied values: 2/2**n."""
    return Fraction(2, 2 ** n)


def smallest_d_passing(alpha):
    """Scan d = 1,2,... for the first whose minimum attainable p is <= alpha."""
    d = 1
    while min_attainable_p(d) > Fraction(alpha).limit_denominator(10 ** 12):
        d += 1
    return d


def walsh_ci(x, level=Fraction(95, 100)):
    """Exact distribution-free CI for the median via Walsh averages.

    The rank of the Walsh average bounding the CI is found by scanning the
    exact null (untied, n values) for the largest k with P(W <= k-1) <= a/2.
    Returns (median_estimate, lo, hi, achieved_coverage) or Nones if n < 2.
    """
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n < 2:
        return None, None, None, None
    walsh = np.sort([(a + b) / 2 for a, b in combinations_with_replacement(x, 2)])
    _, counts = signed_rank_null(np.arange(1, n + 1) * 2)   # untied ranks, doubled
    probs = [Fraction(int(c), 2 ** n) for c in counts]
    # W (undoubled) = W2/2; only even W2 values are attainable here.
    cdf, acc = {}, Fraction(0)
    for w2, pr in enumerate(probs):
        acc += pr
        if w2 % 2 == 0:
            cdf[w2 // 2] = acc
    half = (1 - level) / 2
    k = 0
    for w in sorted(cdf):
        if cdf[w] <= half:
            k = w + 1
        else:
            break
    if k == 0:
        return float(np.median(walsh)), None, None, Fraction(1)
    lo, hi = walsh[k - 1], walsh[len(walsh) - k]
    coverage = 1 - 2 * cdf[k - 1]
    return float(np.median(walsh)), float(lo), float(hi), coverage


def holm(pvals):
    """Holm step-down adjusted p-values (same order as input). Accepts Fractions or floats."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [None] * m
    running = 0
    for rank, i in enumerate(order):
        val = min(1, (m - rank) * pvals[i])
        running = max(running, val)
        adj[i] = running
    return adj


def concordance(del_sig, del_median, ovx_sig, ovx_median):
    """Registered s.5b dose concordance."""
    if del_sig and ovx_sig:
        return "COHERENT" if np.sign(del_median) == -np.sign(ovx_median) and del_median != 0 else "INCOHERENT"
    if del_sig or ovx_sig:
        return "UNRESOLVED"
    return "NONE"


def panel_a_status(eligible, controls_ok, del_sig, del_median, conc, ambient_flagged):
    """Registered s.5c Panel A status. ambient_flagged covers flagged, CURATED_LINEAGE_FOREIGN and UNDETERMINED."""
    if not controls_ok:
        return "NOT_ESTIMABLE_CONTROLS"
    if not eligible:
        return "NOT_RUN"
    if not del_sig:
        return "OPEN"
    if conc == "INCOHERENT":
        return "DOSE_INCOHERENT"
    if conc == "UNRESOLVED":
        return "DELETION_ONLY"
    if conc == "COHERENT":
        if del_median > 0:
            return "REPLICATED_AMBIENT" if ambient_flagged else "REPLICATED"
        return "REVERSED"
    raise ValueError(f"unreachable: del_sig with concordance {conc!r}")


def panel_b_status(eligible, controls_ok, del_sig, del_median, conc):
    if not controls_ok:
        return "NOT_ESTIMABLE_CONTROLS"
    if not eligible:
        return "NOT_RUN"
    if not del_sig:
        return "OPEN"
    if conc == "INCOHERENT":
        return "DOSE_INCOHERENT"
    if conc == "UNRESOLVED":
        return "DELETION_ONLY"
    return "T_CELL_SIGNAL_TOWARD" if del_median > 0 else "T_CELL_SIGNAL_AWAY"


def make_row(panel, gene, status, involves_lusc=False, **fields):
    row = {"panel": panel, "gene": gene, "status": status,
           "design_warrant": WARRANT_LUSC if involves_lusc else WARRANT_PAIRED, **fields}
    if panel == "A" and status in NON_REPLICATION:
        row["claim_qualifier"] = PANEL_A_QUALIFIER
    return row


def validate_rows(rows):
    """Refuse to render rows that drop a registered qualifier (Amendment 3, 3.6 and 3.8)."""
    for r in rows:
        if "design_warrant" not in r or r["design_warrant"] not in (WARRANT_PAIRED, WARRANT_LUSC):
            raise ValueError(f"row {r.get('gene')}: missing/invalid design_warrant")
        if r["panel"] == "A" and r["status"] in NON_REPLICATION and r.get("claim_qualifier") != PANEL_A_QUALIFIER:
            raise ValueError(f"row {r.get('gene')}: Panel A {r['status']} without the registered claim_qualifier")
    return True
