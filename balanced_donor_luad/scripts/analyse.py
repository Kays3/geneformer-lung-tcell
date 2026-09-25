"""Phase 7: registered analysis of the balanced-donor ISP (PHASE5_ISP_REGISTRATION.md s.3-s.7,
Amendments 3, 3c, 3e, 3f).

Reads ONLY completed Phase 6 output: <root>/<op>/<gene>/<donor>.complete.json plus the raw
per-cell pickle written beside it. Refuses to run on a partial arm (registration s.8).

Per (gene g, operation op, donor d), from the raw per-cell pickle (never Shift_to_goal_end):
  s_op(g,d) = mean per-cell cosine shift toward the goal state over d's perturbed cells
  a_op(g,d) = s_op(g,d) - median over controls c in C(g) of s_op(c,d)
Estimable (g,d): gene token in >= MIN_CELLS of d's 100 cells (marker n_token_cells).
Eligible (g): estimable in >= d_min donors (A: 10, B: 11).
Test: exact two-sided Wilcoxon over estimable donors; Holm within each of 4 families;
dose concordance; registered status tables; S3 downgrade; panel-level secondaries;
S1 (global goal), S2 (Leader_Merad vs rest), S4 (cell-weighted), and the Amendment 3f
accuracy/effect correlations (descriptive only).

Two rules the registration leaves open are PARAMETERS (see RULES below). They must be fixed
by a dated amendment before any real Phase 6 output is read.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
from dataclasses import dataclass, field

import numpy as np

import stats_core as sc

OPS = ("delete", "overexpress")
GOAL = "normal_adjacent"                 # registered primary goal (donor's own centroid)
GOAL_S1 = "global_normal_adjacent"       # S1 sensitivity
ALPHA = 0.05
MIN_CELLS = 10
D_MIN = {"A": 10, "B": 11}
SECONDARY_FLOOR = 12
POSITIVE = {"REPLICATED", "REPLICATED_AMBIENT", "T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY"}

# Open rules (to be fixed by amendment before unblinding). Defaults are the PROPOSAL.
RULES = {
    # Minimum token-positive cells for a CONTROL to count in donor d (Amendment 3g: same as the gene).
    "control_min_cells": 10,
    # A control contributes to donor d's median only if it is estimable in d (same MIN_CELLS rule
    # as the gene). a(g,d) is undefined if fewer than this many of the 20 controls qualify.
    "min_controls_per_donor": 10,
    # "panel": Holm over the full registered family (15 / 36), untested members entered as p = 1,
    # which is the m the registered d_min table was derived from. "tested": only genes tested.
    "holm_m": "panel",
}
# Amendment 3g: the alternatives are reported as a SENSITIVITY LINE only; no status ever changes on them.
ALT_RULES = {
    "rule1_any_control": {"control_min_cells": 1, "min_controls_per_donor": 1, "holm_m": "panel"},
    "rule2_holm_tested": {"control_min_cells": 10, "min_controls_per_donor": 10, "holm_m": "tested"},
}


# ----------------------------------------------------------------------------- loading
class IncompleteRun(RuntimeError):
    pass


def load_call(root, op, gene, donor):
    """Return dict(status, n_token_cells, shifts={state: np.array}) or raise IncompleteRun."""
    safe = donor.replace("/", "_")
    marker = os.path.join(root, op, gene, f"{safe}.complete.json")
    if not os.path.exists(marker):
        raise IncompleteRun(f"missing marker {op}/{gene}/{donor}")
    try:
        rec = json.load(open(marker))
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise IncompleteRun(f"unparseable marker {op}/{gene}/{donor}: {e}")
    if rec.get("status") == "no_token_cells":
        return {"status": "no_token_cells", "n_token_cells": 0, "shifts": {}}
    if rec.get("status") != "done":
        raise IncompleteRun(f"marker status {rec.get('status')!r} for {op}/{gene}/{donor}")
    fs = glob.glob(os.path.join(root, op, gene, f"in_silico_{op}_{safe}_cell_embs_dict_*_raw.pickle"))
    if len(fs) != 1:
        raise IncompleteRun(f"expected 1 pickle for {op}/{gene}/{donor}, found {len(fs)}")
    d = pickle.load(open(fs[0], "rb"))
    shifts = {}
    for state, per_key in d.items():
        vals = list(per_key.values())
        if len(vals) != 1:
            raise IncompleteRun(f"{op}/{gene}/{donor} state {state}: {len(vals)} keys, expected 1")
        shifts[state] = np.asarray(vals[0], dtype=float)
    n = rec.get("n_token_cells")
    if n is None or any(len(v) != n for v in shifts.values()):
        raise IncompleteRun(f"{op}/{gene}/{donor}: per-cell count != marker n_token_cells ({n})")
    return {"status": "done", "n_token_cells": int(n), "shifts": shifts}


# ----------------------------------------------------------------------------- config
@dataclass
class Design:
    panel_a: list                 # [{"ensembl_id","symbol"}], registered order
    panel_b: list                 # perturbable genes
    panel_b_not_run: list         # [{"ensembl_id","symbol","reason"}] e.g. TRAC/TRBC1/TRBC2
    strata: dict                  # stratum -> {"status","controls","members"}
    gene_stratum: dict            # panel gene -> stratum
    ambient_flagged: dict         # Panel A gene -> bool (flagged / curated-foreign / undetermined)
    ambient_label: dict
    donors: list
    donor_fold: dict
    donor_study: dict
    donor_ba: dict = field(default_factory=dict)
    expected_estimable: dict = field(default_factory=dict)   # from the pre-GPU eligibility table
    not_run_genes: set = field(default_factory=set)         # genes deliberately not run (BTG1)


# ----------------------------------------------------------------------------- core
def donor_means(calls, gene, op, state, min_cells=MIN_CELLS):
    """{donor: (s, n_cells)} for donors where the gene has >= min_cells token-positive cells."""
    out = {}
    for d, c in calls[(op, gene)].items():
        if c["status"] == "done" and c["n_token_cells"] >= min_cells:
            out[d] = (float(np.mean(c["shifts"][state])), c["n_token_cells"])
    return out


def adjusted_values(calls, design, gene, op, state, rules):
    """a_op(g,d) for estimable donors with enough estimable controls. Returns {d: a}, info."""
    st = design.strata[design.gene_stratum[gene]]
    g_vals = donor_means(calls, gene, op, state)
    ctrl_vals = {c: donor_means(calls, c, op, state, rules["control_min_cells"]) for c in st["controls"]}
    a, n_ctrl_used, dropped = {}, {}, []
    for d, (s, _) in g_vals.items():
        cs = [ctrl_vals[c][d][0] for c in st["controls"] if d in ctrl_vals[c]]
        if len(cs) < rules["min_controls_per_donor"]:
            dropped.append(d)
            continue
        a[d] = s - float(np.median(cs))
        n_ctrl_used[d] = len(cs)
    n_ctrl_all = {d: sum(d in ctrl_vals[c] for c in st["controls"]) for d in g_vals}
    return a, {"n_estimable": len(g_vals), "n_controls_used": n_ctrl_used, "n_controls_qualifying": n_ctrl_all,
               "dropped_few_controls": dropped}


def test_arm(a):
    x = np.array(list(a.values()), dtype=float)
    w = sc.wilcoxon_exact(x)
    med, lo, hi, cov = sc.walsh_ci(x) if len(x) >= 2 else (float(np.median(x)) if len(x) else None, None, None, None)
    return {"n": len(x), "n_used": w["n_used"], "p": w["p"], "median": med, "ci_lo": lo, "ci_hi": hi,
            "ci_coverage": None if cov is None else float(cov), "n_pos": int((x > 0).sum()), "n_neg": int((x < 0).sum())}


def holm_family(pvals_by_gene, family_genes, rules):
    """Holm-adjust. pvals_by_gene only holds tested genes; 'panel' pads untested members with p = 1."""
    genes = list(pvals_by_gene)
    ps = [pvals_by_gene[g] for g in genes]
    if rules["holm_m"] == "panel":
        pad = [g for g in family_genes if g not in pvals_by_gene]
        adj = sc.holm(ps + [sc.Fraction(1)] * len(pad))[: len(genes)]
    elif rules["holm_m"] == "tested":
        adj = sc.holm(ps) if ps else []
    else:
        raise ValueError(rules["holm_m"])
    return dict(zip(genes, adj))


def sign_test(n_pos, n_neg):
    from math import comb
    n = n_pos + n_neg
    if n == 0:
        return sc.Fraction(1)
    k = min(n_pos, n_neg)
    return min(sc.Fraction(1), sc.Fraction(2 * sum(comb(n, i) for i in range(k + 1)), 2 ** n))


def spearman(x, y):
    from scipy.stats import spearmanr
    if len(x) < 3:
        return None
    r = spearmanr(x, y).correlation
    return None if r is None or np.isnan(r) else float(r)


def analyse(calls, design, rules=RULES, state=GOAL):
    """Full registered analysis for one goal state. Returns dict with rows and secondaries."""
    panels = {"A": [g["ensembl_id"] for g in design.panel_a], "B": [g["ensembl_id"] for g in design.panel_b]}
    symbol = {g["ensembl_id"]: g["symbol"] for g in design.panel_a + design.panel_b}
    arms, info, eligible, controls_ok, not_run_reason = {}, {}, {}, {}, {}
    for P, genes in panels.items():
        for g in genes:
            stratum = design.gene_stratum.get(g)
            st = design.strata.get(stratum) if stratum else None
            if st is None:
                # ineligible before any GPU run (s.4): never perturbed, kept as a NOT_RUN row
                k = design.expected_estimable.get(g)
                not_run_reason[g] = f"estimable donors = {k} < d_min = {D_MIN[P]}"
                eligible[g], controls_ok[g] = False, True
                continue
            controls_ok[g] = st["status"] == "eligible" and len(st.get("controls") or []) == 20
            if not controls_ok[g] or g in design.not_run_genes:
                eligible[g] = False
                continue
            for op in OPS:
                a, inf = adjusted_values(calls, design, g, op, state, rules)
                arms[(g, op)], info[(g, op)] = a, inf
            n_est = info[(g, "delete")]["n_estimable"]
            if info[(g, "overexpress")]["n_estimable"] != n_est:
                raise ValueError(f"{g}: estimable donors differ between operations")
            if g in design.expected_estimable and design.expected_estimable[g] != n_est:
                raise ValueError(f"{g}: estimable donors {n_est} != pre-GPU table {design.expected_estimable[g]}")
            eligible[g] = n_est >= D_MIN[P]            # registered s.4 rule: estimable donors
            if not eligible[g]:
                not_run_reason[g] = f"estimable donors = {n_est} < d_min = {D_MIN[P]}"
    results = {}
    for P, genes in panels.items():
        for op in OPS:
            tested = {g: test_arm(arms[(g, op)]) for g in genes if eligible[g]}
            adj = holm_family({g: r["p"] for g, r in tested.items()}, genes, rules)
            m = len(genes) if rules["holm_m"] == "panel" else len(tested)
            for g, r in tested.items():
                r["p_holm"] = adj[g]
                r["holm_m"] = m
                r["sig"] = adj[g] <= sc.Fraction(ALPHA).limit_denominator(1000)
                results[(g, op)] = r
    rows = []
    for P, genes in panels.items():
        for g in genes:
            fields = {"symbol": symbol[g], "stratum": design.gene_stratum.get(g)}
            if P == "A":
                fields["ambient"] = design.ambient_label.get(g)
            else:
                fields["ambient"] = "CURATED_T_CELL"
            if not eligible[g]:
                if not controls_ok[g]:
                    status, reason = "NOT_ESTIMABLE_CONTROLS", "fewer than 20 matched controls in its stratum"
                else:
                    status, reason = "NOT_RUN", not_run_reason.get(g, "not run")
                rows.append(sc.make_row(P, g, status, reason=reason,
                                        estimable_donors=design.expected_estimable.get(g), **fields))
                continue
            dr, orr = results[(g, "delete")], results[(g, "overexpress")]
            conc = sc.concordance(dr["sig"], dr["median"], orr["sig"], orr["median"])
            if P == "A":
                status = sc.panel_a_status(True, True, dr["sig"], dr["median"], conc, design.ambient_flagged[g])
            else:
                status = sc.panel_b_status(True, True, dr["sig"], dr["median"], conc)
            # S3: fold heterogeneity can only downgrade a positive status to OPEN
            fold_med = {}
            for d, v in arms[(g, "delete")].items():
                fold_med.setdefault(design.donor_fold[d], []).append(v)
            fold_med = {k: float(np.median(v)) for k, v in sorted(fold_med.items())}
            opposite = sum(1 for v in fold_med.values() if dr["median"] and np.sign(v) == -np.sign(dr["median"]))
            s3 = None
            if status in POSITIVE and opposite > 1:
                s3, status = f"fold-heterogeneous (downgraded from {status})", "OPEN"
            fields.update({
                "concordance": conc, "s3": s3, "fold_medians_delete": fold_med,
                **{f"{tag}_{k}": (float(v) if isinstance(v, sc.Fraction) else v)
                   for tag, r in (("del", dr), ("ovx", orr)) for k, v in r.items()},
                "dropped_few_controls": {op: info[(g, op)]["dropped_few_controls"] for op in OPS},
                "controls_qualifying_per_donor": {op: info[(g, op)]["n_controls_qualifying"] for op in OPS},
            })
            rows.append(sc.make_row(P, g, status, **fields))
    for nr in design.panel_b_not_run:
        rows.append(sc.make_row("B", nr["ensembl_id"], "NOT_RUN", reason=nr["reason"], symbol=nr["symbol"],
                                ambient="CURATED_T_CELL"))
    sc.validate_rows(rows)
    return {"rows": rows, "arms": arms, "results": results, "eligible": eligible,
            "secondaries": secondaries(rows, arms, design, panels)}


def secondaries(rows, arms, design, panels):
    out = {}
    for P, genes in panels.items():
        pr = [r for r in rows if r["panel"] == P and r["gene"] in genes and "del_median" in r]
        subsets = {"all": pr}
        if P == "A":
            subsets["ambient_flagged"] = [r for r in pr if design.ambient_flagged.get(r["gene"])]
            subsets["not_flagged"] = [r for r in pr if not design.ambient_flagged.get(r["gene"])]
        for name, sub in subsets.items():
            n = len(sub)
            key = f"{P}_{name}"
            if n < SECONDARY_FLOOR:
                out[key] = {"status": f"not run: n_eligible = {n} < {SECONDARY_FLOOR}"}
                continue
            pos = sum(r["del_median"] > 0 for r in sub); neg = sum(r["del_median"] < 0 for r in sub)
            p = sign_test(pos, neg)
            out[key] = {"n_eligible": n, "n_toward": pos, "n_away": neg, "sign_test_p": float(p),
                        "caveat": "genes in one stratum share controls (Amendment 3c.2)"}
        out[f"{P}_n_coherent"] = sum(r.get("concordance") == "COHERENT" for r in pr)
    return out


def sensitivities(calls, design, primary, rules=RULES):
    """S1 (global goal), S2 (Leader_Merad vs rest), S4 (cell-weighted), 3f (BA vs effect). Never change a status."""
    s1 = analyse(calls, design, rules, state=GOAL_S1)
    out = {"S1_global_goal": {r["gene"]: {"status": r["status"], "del_median": r.get("del_median")} for r in s1["rows"]}}
    s2, s4, f3 = {}, {}, {}
    for (g, op), a in primary["arms"].items():
        if not primary["eligible"].get(g):
            continue
        lm = [v for d, v in a.items() if design.donor_study[d].startswith("Leader_Merad")]
        ot = [v for d, v in a.items() if not design.donor_study[d].startswith("Leader_Merad")]
        s2[f"{g}|{op}"] = {"leader_merad_median": float(np.median(lm)) if lm else None, "n_lm": len(lm),
                           "others_median": float(np.median(ot)) if ot else None, "n_others": len(ot)}
        cells = np.concatenate([calls[(op, g)][d]["shifts"][GOAL] for d in a]) if a else np.array([])
        s4[f"{g}|{op}"] = {"cell_weighted_mean_shift": float(cells.mean()) if cells.size else None,
                           "n_cells": int(cells.size), "donor_median_a": primary["results"][(g, op)]["median"],
                           "note": "pooled per-cell raw shift (not control-adjusted); InSilicoPerturberStats not run"}
        if design.donor_ba and primary["results"][(g, op)]["sig"]:
            ds = sorted(a)
            f3[f"{g}|{op}"] = {"spearman_ba_vs_a": spearman([design.donor_ba[d] for d in ds], [a[d] for d in ds]),
                               "conditioned_on_significance": True}
    if design.donor_ba:
        for op in OPS:
            for P in ("A", "B"):
                genes = [g for g in ([x["ensembl_id"] for x in (design.panel_a if P == "A" else design.panel_b)])
                         if primary["eligible"].get(g)]
                per_d = {}
                for g in genes:
                    for d, v in primary["arms"][(g, op)].items():
                        per_d.setdefault(d, []).append(v)
                ds = sorted(per_d)
                f3[f"panel_{P}|{op}"] = {"spearman_ba_vs_panel_median_a":
                                         spearman([design.donor_ba[d] for d in ds], [float(np.median(per_d[d])) for d in ds]),
                                         "n_donors": len(ds), "conditioned_on_significance": False,
                                         "n_genes_per_donor": {d: len(per_d[d]) for d in ds},
                                         "note": "per-donor value = median of a over a VARYING number of genes"}
    alt = {}
    base = {r["gene"]: r["status"] for r in primary["rows"]}
    for name, ar in ALT_RULES.items():
        other = {r["gene"]: r["status"] for r in analyse(calls, design, ar)["rows"]}
        alt[name] = {"rules": ar, "status_changes": {g: {"registered": base[g], "alternative": other[g]}
                                                    for g in base if other.get(g) != base[g]},
                     "note": "sensitivity line only (Amendment 3g); no registered status changes on it"}
    out.update({"A3g_alternative_rules": alt, "S2_leader_merad_vs_rest": s2, "S4_cell_weighted": s4, "A3f_accuracy_vs_effect": f3,
                "A3f_low_ba_donors": sorted(d for d, b in design.donor_ba.items() if b < 0.65)})
    return out


# ----------------------------------------------------------------------------- driver
def genes_to_load(design):
    """Genes Phase 6 ran: panel genes in an eligible stratum (not in not_run_genes) plus every control.
    Panel genes without a stratum (ineligible before GPU, s.4) were never perturbed and have no output."""
    run_panel = {g for g, s in design.gene_stratum.items()
                 if design.strata[s]["status"] == "eligible" and g not in set(design.not_run_genes)}
    controls = {c for s in design.strata.values() for c in (s.get("controls") or [])}
    return sorted(run_panel | controls)


def load_all(root, genes, donors):
    calls, missing = {}, []
    for op in OPS:
        for g in genes:
            calls[(op, g)] = {}
            for d in donors:
                try:
                    calls[(op, g)][d] = load_call(root, op, g, d)
                except IncompleteRun as e:
                    missing.append(str(e))
    if missing:
        raise IncompleteRun(f"{len(missing)} calls incomplete (a partial arm is never analysed); first: {missing[:5]}")
    return calls


def jsonable(o):
    if isinstance(o, sc.Fraction):
        return float(o)
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(type(o))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--phase6-root", required=True)
    p.add_argument("--design", required=True, help="design.json built by build_design.py")
    p.add_argument("--rules", required=True, help="rules JSON fixed by amendment (min_controls_per_donor, holm_m)")
    p.add_argument("--out", required=True)
    p.add_argument("--host-drift", choices=("true", "false"), required=True,
                   help="from the end-of-run equivalence gate (s.3.5): true stamps every row")
    a = p.parse_args()
    rules = json.load(open(a.rules))
    assert set(rules) == set(RULES), f"rules must set exactly {sorted(RULES)}"
    design = Design(**json.load(open(a.design)))
    design.not_run_genes = set(design.not_run_genes)
    calls = load_all(a.phase6_root, genes_to_load(design), design.donors)
    primary = analyse(calls, design, rules)
    for r in primary["rows"]:
        r["host_drift"] = a.host_drift == "true"
    sens = sensitivities(calls, design, primary, rules)
    os.makedirs(a.out, exist_ok=True)
    json.dump(primary["rows"], open(os.path.join(a.out, "outcome_rows.json"), "w"), indent=1, default=jsonable)
    json.dump(primary["secondaries"], open(os.path.join(a.out, "secondaries.json"), "w"), indent=1, default=jsonable)
    json.dump(sens, open(os.path.join(a.out, "sensitivities.json"), "w"), indent=1, default=jsonable)
    per_donor = {f"{g}|{op}": v for (g, op), v in primary["arms"].items()}
    json.dump(per_donor, open(os.path.join(a.out, "per_donor_adjusted.json"), "w"), indent=1, default=jsonable)
    json.dump({"rules": rules, "n_rows": len(primary["rows"])}, open(os.path.join(a.out, "run_record.json"), "w"), indent=1)
    print(json.dumps({r["gene"]: r["status"] for r in primary["rows"]}, indent=1))


if __name__ == "__main__":
    main()
