"""Registered S100 LUAD analysis (Amendments 6-8), run once after
Stanley's post-run PASS (2026-09-28 ~17:30 JST). CPU only.

Run from sclc_validation/perturbation_workflow/targeted_panel/ with
S100_ISP_DIR set and the same explicit roots as the GPU run.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.argv = ["run_targeted_panel.py", "--run-tag", "s100_luad_analysis_20260928", "--dtype", "bf16"]
sys.path.insert(0, ".")
import run_targeted_panel as rtp  # noqa: E402

sys.path.insert(0, os.environ["S100_ISP_DIR"])
import ambient_stats as ast  # noqa: E402
import retained_rows as rr  # noqa: E402

RUNS = rtp.BF16_BENCH_ROOT / "runs"
MAIN = RUNS / "s100_luad_20260928" / "targeted_panel"
NOOP = RUNS / "s100_luad_noop_20260928" / "targeted_panel"
PRE = RUNS / "s100_luad_preflight_pergene_20260928"
OUT = RUNS / "s100_luad_analysis_20260928"
OUT.mkdir(parents=True, exist_ok=True)
EXPECT = {PRE / "s100_luad_run_panel_pergene_20260928.json": "8a0b36683542e51fa81a655505a348f9286669df9adac22859d5f6406f34afba",
          PRE / "matched_control_sets_pergene.json": "607ac6c9",
          Path("run_targeted_panel.py"): "974535b7b9ad6b49e33e4dfaf2a84e7070f0cc63e64bd8b271233cf77b60f754"}
HIGH = {"S100A2": "ENSG00000196754", "S100B": "ENSG00000160307", "S100A13": "ENSG00000189171", "S100PBP": "ENSG00000116497"}
LOW = {"S100A4": "ENSG00000196154", "S100A6": "ENSG00000197956", "S100A10": "ENSG00000197747", "S100A11": "ENSG00000163191"}
RISK = {"S100A4": 0.0000011, "S100A6": 0.000034, "S100A11": 0.000144, "S100A10": 0.000312,
        "S100PBP": 0.744492, "S100A13": 0.795706, "S100B": 0.859066, "S100A2": 0.955296}
PRIMARY = {**HIGH, **LOW}
LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def as_json(x):
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    raise TypeError(type(x))


for p, want in EXPECT.items():
    got = sha(p)
    say(f"{p.name} sha256={got}")
    assert got.startswith(want), (p, got)

run_panel = json.loads((PRE / "s100_luad_run_panel_pergene_20260928.json").read_text())
sets = json.loads((PRE / "matched_control_sets_pergene.json").read_text())["sets"]
for g, eid in PRIMARY.items():
    assert next(x for x in run_panel["genes"] if x["gene"] == g)["ensembl_id"] == eid, g

rows = rr.build_retained_rows(
    "s100_luad_run_panel_pergene_20260928", run_panel["genes"],
    raw_root=MAIN / "raw", noop_raw_root=NOOP / "raw", paired_eligible_dir=rtp.PAIRED_ELIGIBLE_DIR,
    stats_root=MAIN / "stats", gene_token_dict=rtp.gene_token_dict(), run_id="s100_luad_20260928",
    runner_sha256=EXPECT[Path("run_targeted_panel.py")], panel_sha256=EXPECT[PRE / "s100_luad_run_panel_pergene_20260928.json"],
    sources=("luad",))
rows.to_csv(OUT / "retained_rows.csv", index=False)
say(f"\nretained rows: {len(rows)} (154 genes x 2 goals x 2 operations); sha256={sha(OUT / 'retained_rows.csv')}")
say(rr.status_counts(rows).to_string())
say("no_op_status: " + rows["no_op_status"].value_counts(dropna=False).to_string().replace("\n", "; "))
panel_rows = rows[rows["role"] != "matched_control"]
say("\npanel rows (status / reason / no-op):")
for _, r in panel_rows.iterrows():
    say(f"  {r['gene']:>8s} {r['perturbation_type']:>11s} ->{rr.SLUGS[r['goal_state']]:<6s} {r['status']:<26s} "
        f"shift={r['donor_balanced_shift']} no_op={r['no_op_status']} {r['no_op_score']} reason={r['status_reason']}")

E = ast.compute_E(rows)
E.to_csv(OUT / "E_by_gene_contrast.csv", index=False)
floor = rows[rows["donor_balanced_shift"].notna()]
n_below = int((floor["donor_balanced_shift"].abs() <= ast.NUMERICAL_FLOOR).sum())
say(f"\narms with |donor_balanced_shift| <= numerical floor {ast.NUMERICAL_FLOOR}: {n_below} of {len(floor)}")

results = {}
stratum_by_gene = {g: g for g in PRIMARY}
for goal in (rr.NORMAL, rr.SCLC):
    c = f"luad_to_{rr.SLUGS[goal]}"
    Ec = E[E["goal_state"] == goal].set_index("ensembl_id")["E"]
    E_by_gene = {g: float(Ec.get(eid, np.nan)) for g, eid in PRIMARY.items()}
    controls = {g: {cid: float(Ec.get(cid, np.nan)) for cid in sets[g]["controls"]} for g in PRIMARY}
    n_valid = {g: int(np.sum(~np.isnan(list(v.values())))) for g, v in controls.items()}
    Q = ast.compute_Q_for_contrast(E_by_gene, stratum_by_gene, controls)
    test = ast.exact_group_separation_test(Q, list(HIGH), list(LOW))
    loo = ast.group_separation_leave_one_control_out(E_by_gene, stratum_by_gene, controls, list(HIGH), list(LOW))
    boot = ast.group_separation_bootstrap(E_by_gene, stratum_by_gene, controls, list(HIGH), list(LOW))
    gate = ast.group_separation_stability_gate(test, loo, boot)
    loo.to_csv(OUT / f"loo_{c}.csv", index=False)
    if np.isnan(test["p_exact"]):
        status = "not_estimable"
    elif gate["outcome"] == "positive_stable":
        status = "positive"
    elif gate["outcome"] == "control_draw_sensitive_open":
        status = "control_draw_sensitive_open"
    elif test["p_exact"] <= ast.GROUP_SEPARATION_ALPHA and test["direction"] == "low_above_high":
        status = "opposite_direction"
    else:
        status = "negative"
    rho = ast.primary_test(Q, RISK, list(PRIMARY))
    logo = ast.leave_one_gene_out_check(Q, RISK, list(PRIMARY), "S100A2")
    results[c] = {"status": status, "E": E_by_gene, "Q": Q, "n_valid_controls": n_valid,
                  "primary_test": test, "stability_gate": gate,
                  "loo_n_positive": int(loo["positive"].sum()),
                  "loo_direction_counts": loo["direction"].value_counts().to_dict(),
                  "loo_p_range": [float(loo["p_exact"].min()), float(loo["p_exact"].max())],
                  "descriptive_spearman_Q_vs_risk": {k: v for k, v in rho.items() if not isinstance(v, np.ndarray)},
                  "descriptive_leave_one_gene_out": {k: v for k, v in logo.items() if not isinstance(v, (np.ndarray, list))}}
    say(f"\n== {c}: STATUS {status}")
    for g in PRIMARY:
        say(f"  {g:>8s} {'high' if g in HIGH else 'low ':>4s} E={E_by_gene[g]:.6f} Q={Q[g]:.2f} valid_controls={n_valid[g]}")
    say(f"  exact MWU: U={test['statistic']} p={test['p_exact']:.4f} ({round(test['p_exact'] * 70)}/70) "
        f"direction={test['direction']} complete_separation={test['complete_separation']} positive={test['positive']}")
    say(f"  stability: LOO positive {int(loo['positive'].sum())}/{len(loo)}, bootstrap positive fraction "
        f"{gate['bootstrap_positive_fraction']:.4f}, outcome={gate['outcome']}")
    say(f"  descriptive rho(Q, risk)={rho.get('rho')} p={rho.get('p_exact')}")

(OUT / "primary_results.json").write_text(json.dumps(results, indent=2, default=as_json) + "\n")
say(f"\nprimary_results.json sha256={sha(OUT / 'primary_results.json')}")
(OUT / "analysis_log.txt").write_text("\n".join(LOG) + "\n")
