"""Step 0.3: ambient classification for the panel genes (registration s.1c), non-circular.

Features are computed exactly as in
sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py
(log1p-CPM; detect_frac, mean_expr, log subtype-F, log donor-F, library-size
Spearman, breadth/depth, subtype/donor), on this cohort's training-pool T cells.

What differs, and why (registered s.1c):
  * VALIDITY by leave-one-anchor-out (LOAO): each anchor is scored by a model
    fitted without it. valid iff LOAO AUC >= 0.8. The original CV-AUC trains on
    the anchors it scores (recorded 1.0 +/- 0.0), so it can essentially never fail.
  * THRESHOLD = 25th percentile of the ambient anchors' LOAO scores.
  * NON-ANCHOR genes are scored by the model fitted on all anchors (they are
    not in its training set, so that is not circular) and flagged if >= threshold.
  * ANCHOR genes are labelled by curation: ambient anchors -> CURATED_LINEAGE_FOREIGN
    (treated as flagged), T-cell anchors -> CURATED_T_CELL (not scored).
  * If invalid, every non-anchor gene is UNDETERMINED (treated as flagged).
  * A non-anchor panel gene below the 0.5% detection filter cannot be scored:
    it is UNDETERMINED (treated as flagged). This case was not named in s.1c;
    it is resolved in the conservative direction (it can only prevent a plain
    REPLICATED) and reported as such.
"""
import argparse
import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

COLUMNS = ["log_subtype_f", "log_donor_f", "libsize_corr", "breadth_over_depth", "subtype_over_donor", "mean_expr"]
MIN_DETECT = 0.005
VALID_AUC = 0.8


def anchor_lists(diagnostic_path):
    s = open(diagnostic_path).read()
    get = lambda name: re.findall(r'"([^"]+)"', re.search(name + r" = \[(.*?)\]", s, re.S).group(1))
    return get("KNOWN_AMBIENT"), get("KNOWN_TCELL")


def group_f_statistic(counts, labels):
    """Verbatim logic of the original diagnostic's group_f_statistic."""
    codes = pd.Categorical(labels).codes
    n_groups = int(codes.max()) + 1
    n_cells, _ = counts.shape
    ind = sp.csr_matrix((np.ones(n_cells), (codes, np.arange(n_cells))), shape=(n_groups, n_cells))
    gn = np.asarray(ind.sum(axis=1)).ravel()
    gs = np.asarray(sp.csr_matrix(ind @ counts).todense())
    gss = np.asarray(sp.csr_matrix(ind @ counts.multiply(counts)).todense())
    gm = gs / gn[:, None]
    grand = np.asarray(counts.mean(axis=0)).ravel()
    between = (gn[:, None] * (gm - grand) ** 2).sum(axis=0) / max(n_groups - 1, 1)
    within = (gss.sum(axis=0) - (gn[:, None] * gm ** 2).sum(axis=0)) / max(n_cells - n_groups, 1)
    return between / np.maximum(within, 1e-12)


def features(adata):
    raw = sp.csr_matrix(adata.X)
    tot = np.asarray(raw.sum(axis=1)).ravel()
    norm = sp.csr_matrix(sp.diags(1e4 / np.maximum(tot, 1)) @ raw)
    norm.data = np.log1p(norm.data)
    n = adata.n_obs
    detect = np.asarray((raw > 0).sum(axis=0)).ravel() / n
    mean_expr = np.asarray(norm.mean(axis=0)).ravel()
    sub_f = group_f_statistic(norm, adata.obs["celltype"])
    don_f = group_f_statistic(norm, adata.obs["individual"])
    rl = pd.Series(tot).rank().to_numpy(); rl = (rl - rl.mean()) / rl.std()
    cm = np.asarray(norm.mean(axis=0)).ravel()
    num = np.asarray(norm.T @ rl).ravel() - n * cm * rl.mean()
    ss = np.asarray(norm.multiply(norm).sum(axis=0)).ravel()
    lib = num / np.maximum(np.sqrt(np.maximum(ss - n * cm ** 2, 1e-12)) * np.sqrt(n), 1e-12)
    f = pd.DataFrame({"ensembl_id": adata.var["ensembl_id"].astype(str).to_numpy(), "detect_frac": detect,
                      "mean_expr": mean_expr, "log_subtype_f": np.log1p(sub_f), "log_donor_f": np.log1p(don_f),
                      "libsize_corr": lib})
    f["breadth_over_depth"] = f.detect_frac / np.maximum(f.mean_expr, 1e-6)
    f["subtype_over_donor"] = f.log_subtype_f / np.maximum(f.log_donor_f, 1e-6)
    return f


def fit(X, y):
    sc = StandardScaler().fit(X)
    m = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(X), y)
    return lambda Z: m.predict_proba(sc.transform(Z))[:, 1]


def classify(feat, amb_ens, tc_ens, panel_a_ens):
    scored = feat[feat.detect_frac >= MIN_DETECT].reset_index(drop=True)
    lab = scored[scored.ensembl_id.isin(amb_ens | tc_ens)].reset_index(drop=True)
    lab["is_ambient"] = lab.ensembl_id.isin(amb_ens).astype(int)
    loao = np.empty(len(lab))
    for i in range(len(lab)):
        keep = np.arange(len(lab)) != i
        loao[i] = fit(lab.loc[keep, COLUMNS], lab.loc[keep, "is_ambient"])(lab.loc[[i], COLUMNS])[0]
    lab["loao_score"] = loao
    auc = float(roc_auc_score(lab.is_ambient, lab.loao_score))
    valid = auc >= VALID_AUC
    threshold = float(np.percentile(lab.loc[lab.is_ambient == 1, "loao_score"], 25))
    full = fit(lab[COLUMNS], lab.is_ambient)
    scored["ambient_risk"] = full(scored[COLUMNS])
    rows = []
    for e in panel_a_ens:
        if e in amb_ens:
            label, risk = "CURATED_LINEAGE_FOREIGN", None
        elif e in tc_ens:
            label, risk = "CURATED_T_CELL", None
        elif not valid:
            label, risk = "UNDETERMINED (diagnostic invalid: LOAO AUC < 0.8)", None
        elif e not in set(scored.ensembl_id):
            label, risk = "UNDETERMINED (detection < 0.5%, not scorable)", None
        else:
            risk = float(scored.loc[scored.ensembl_id == e, "ambient_risk"].iloc[0])
            label = "AMBIENT_FLAGGED" if risk >= threshold else "NOT_FLAGGED"
        rows.append({"ensembl_id": e, "ambient_label": label, "ambient_risk": risk,
                     "treated_as_flagged": label not in ("NOT_FLAGGED", "CURATED_T_CELL")})
    summary = {"n_cells": int(0), "n_genes_scored": int(len(scored)),
               "n_ambient_anchors_scored": int(lab.is_ambient.sum()), "n_tcell_anchors_scored": int((1 - lab.is_ambient).sum()),
               "loao_auc": auc, "valid": bool(valid), "threshold_p25_ambient_loao": threshold,
               "median_loao_ambient": float(lab.loc[lab.is_ambient == 1, "loao_score"].median()),
               "median_loao_tcell": float(lab.loc[lab.is_ambient == 0, "loao_score"].median())}
    return pd.DataFrame(rows), summary, lab


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--slim", required=True)
    p.add_argument("--symbols", required=True, help="csv ensembl_id,symbol from the source atlas var")
    p.add_argument("--diagnostic", required=True)
    p.add_argument("--panel-a", required=True)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()
    adata = ad.read_h5ad(a.slim)
    sym = pd.read_csv(a.symbols)
    s2e = sym.groupby("symbol").ensembl_id.agg(list)
    ka, kt = anchor_lists(a.diagnostic)
    amb = {s2e[g][0] for g in ka if g in s2e.index and len(s2e[g]) == 1}
    tc = {s2e[g][0] for g in kt if g in s2e.index and len(s2e[g]) == 1}
    panel = [g["ensembl_id"] for g in json.load(open(a.panel_a))["genes"]]
    feat = features(adata)
    table, summary, lab = classify(feat, amb, tc, panel)
    summary["n_cells"] = int(adata.n_obs)
    table = table.merge(sym.drop_duplicates("ensembl_id"), on="ensembl_id", how="left")
    table.to_csv(a.out_prefix + "_panelA.csv", index=False)
    lab.to_csv(a.out_prefix + "_anchors_loao.csv", index=False)
    json.dump(summary, open(a.out_prefix + "_summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1)); print(table.to_string(index=False))


if __name__ == "__main__":
    main()
