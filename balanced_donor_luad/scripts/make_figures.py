"""Figures for RESULTS.md and the slides, from COMMITTED Phase 7 outputs only.

Every registered row (15 Panel A + 39 Panel B = 54) must appear in the outcome table and in the panel
figures; the script refuses (exit 1) otherwise. Palette: Okabe-Ito (colour-blind safe).
Usage: python make_figures.py --base balanced_donor_luad --out balanced_donor_luad/figures
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

STATUS_COLOUR = {
    "T_CELL_SIGNAL_TOWARD": "#0072B2", "T_CELL_SIGNAL_AWAY": "#D55E00", "DELETION_ONLY": "#E69F00",
    "DOSE_INCOHERENT": "#CC79A7", "REPLICATED": "#009E73", "REPLICATED_AMBIENT": "#56B4E9", "REVERSED": "#D55E00",
    "OPEN": "#7F7F7F", "NOT_RUN": "#D9D9D9", "NOT_ESTIMABLE_CONTROLS": "#BDBDBD",
}
STATUS_LABEL = {"T_CELL_SIGNAL_TOWARD": "T-cell signal, toward normal", "T_CELL_SIGNAL_AWAY": "T-cell signal, away from normal",
                "DELETION_ONLY": "deletion only", "DOSE_INCOHERENT": "dose-incoherent", "OPEN": "open (no call)",
                "NOT_RUN": "not run", "NOT_ESTIMABLE_CONTROLS": "no controls"}
NAVY, MUTED = "#10243C", "#4A5768"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight", "axes.titlecolor": NAVY})


def load(base):
    j = lambda p: json.load(open(os.path.join(base, p)))
    d = {"rows": j("phase7_results/outcome_rows.json"), "sens": j("phase7_results/sensitivities.json"),
         "sec": j("phase7_results/secondaries.json"), "per_donor": j("phase7_results/per_donor_adjusted.json"),
         "ann": j("phase7_results/row_annotations.json"), "fig": j("phase7_results/figure_data.json"),
         "gate": j("phase4_results/classifier_gate.json"), "cohort": j("cohort/COHORT_DEFINITION.json"),
         "folds": j("cohort/folds.json"), "pa": j("registration/panel_A.json"), "pb": j("registration/panel_B.json"),
         "design": j("phase7_results/design.json")}
    return d


def check_rows(d):
    want = [("A", g["ensembl_id"]) for g in d["pa"]["genes"]] + [("B", g["ensembl_id"]) for g in d["pb"]["genes"]]
    got = [(r["panel"], r["gene"]) for r in d["rows"]]
    missing, extra = set(want) - set(got), set(got) - set(want)
    if missing or extra or len(got) != len(want) or len(want) != 54:
        sys.exit(f"REFUSED: rows do not match the 54 registered ({len(got)} rows; missing {missing}; extra {extra})")


def estimable(d, r):
    return r.get("estimable_donors") or d["design"]["expected_estimable"].get(r["gene"])


def legend_for(ax, statuses, **kw):
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=STATUS_COLOUR[s], label=STATUS_LABEL.get(s, s)) for s in statuses],
              frameon=False, fontsize=7.5, **kw)


# ----------------------------------------------------------------------------------------------- 1
def fig_design(d, out):
    counts = d["cohort"]["eligible_cells_per_donor_origin"]
    donors = sorted({k.split("|")[0] for k in counts}, key=lambda x: -counts[f"{x}|tumor_primary"])
    tu = np.array([counts[f"{x}|tumor_primary"] for x in donors]); no = np.array([counts[f"{x}|normal_adjacent"] for x in donors])
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw={"width_ratios": [3, 1], "wspace": 0.55})
    x = np.arange(len(donors))
    ax[0].bar(x - 0.2, tu, 0.4, color="#B8465E", label="tumour (primary)")
    ax[0].bar(x + 0.2, no, 0.4, color="#0B6F6A", label="adjacent normal")
    ax[0].axhline(100, color=NAVY, lw=1, ls="--"); ax[0].text(-0.5, 72, "analysis cap: 100 per tissue", fontsize=7.5, color=NAVY); ax[0].set_ylim(60, None)
    ax[0].set_yscale("log"); ax[0].set_xticks([]); ax[0].set_xlabel(f"{len(donors)} donors (sorted by tumour T cells)")
    ax[0].set_ylabel("eligible T cells (log)")
    ax[0].legend(frameon=False, fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    ax[0].set_title(f"Before the cap: {tu.max() / tu.min():.1f}x (tumour) and {no.max() / no.min():.1f}x (normal) between donors")
    studies = d["folds"]["donor_study"]; sc = {}
    for s in studies.values(): sc[s] = sc.get(s, 0) + 1
    names = sorted(sc, key=lambda s: -sc[s])
    ax[1].barh(range(len(names)), [sc[s] for s in names], color=MUTED)
    ax[1].set_yticks(range(len(names))); ax[1].set_yticklabels([s.replace("_", " ") for s in names], fontsize=7.5)
    ax[1].invert_yaxis(); ax[1].set_xlabel("donors"); ax[1].set_title("5 studies, 10x only")
    fig.suptitle("Design: every donor contributes both tissues; after the cap each donor weighs exactly 1/43", color=NAVY, fontsize=10.5, y=1.03)
    fig.savefig(os.path.join(out, "fig1_design.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 2
def fig_gate(d, out):
    g = d["gate"]; fold = d["folds"]["donor_test_fold"]; ba = g["per_donor_balanced_accuracy"]
    fig, ax = plt.subplots(figsize=(7, 3.2))
    rng = np.random.default_rng(0)
    for k in range(5):
        v = [ba[x] for x in ba if fold[x] == k]
        ax.scatter(k + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color="#2C6E9B", alpha=0.8)
        ax.hlines(g["per_fold_balanced_accuracy"][str(k)], k - 0.3, k + 0.3, color=NAVY, lw=2)
    v = list(ba.values())
    ax.scatter(5 + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color="#0B6F6A", alpha=0.8)
    ax.hlines(g["pooled_balanced_accuracy"], 4.7, 5.3, color=NAVY, lw=2)
    ax.axhline(0.60, color="#B8465E", ls="--", lw=1); ax.text(-0.45, 0.61, "gate 0.60", color="#B8465E", fontsize=7.5)
    ax.axhline(0.5, color=MUTED, ls=":", lw=1); ax.text(-0.45, 0.51, "chance", color=MUTED, fontsize=7.5)
    ax.set_xticks(range(6)); ax.set_xticklabels([f"fold {k}" for k in range(5)] + ["pooled"])
    ax.set_ylim(0.45, 1.02); ax.set_ylabel("held-out balanced accuracy (per donor)")
    ax.set_title(f"Classifier gate PASS: pooled {g['pooled_balanced_accuracy']:.3f}; {g['donors_above_0.5']}/{g['n_donors']} donors above chance "
                 f"(sign test p = {g['sign_test_p_float']:.1e})")
    fig.savefig(os.path.join(out, "fig2_classifier_gate.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 3
def fig_panel_a(d, out):
    rows = [r for r in d["rows"] if r["panel"] == "A"]
    rows.sort(key=lambda r: -(estimable(d, r) or 0))
    fig, ax = plt.subplots(figsize=(8, 3.6))
    y = np.arange(len(rows))
    est = [estimable(d, r) or 0 for r in rows]
    ax.barh(y, est, color=[STATUS_COLOUR[r["status"]] for r in rows], edgecolor="white")
    ax.axvline(10, color="#B8465E", ls="--", lw=1); ax.text(10.4, -0.9, "d_min = 10 donors", color="#B8465E", fontsize=7.5)
    ax.axvline(43, color=MUTED, ls=":", lw=1); ax.text(43.4, -0.9, "all 43", color=MUTED, fontsize=7.5)
    for yi, r, e in zip(y, rows, est):
        tag = {"AMBIENT_FLAGGED": " [ambient]", "CURATED_LINEAGE_FOREIGN": " [lineage-foreign]"}.get(r.get("ambient"), "")
        ax.text(e + 0.4, yi, f"{STATUS_LABEL[r['status']]}{tag}", va="center", fontsize=7.5, color=MUTED)
    ax.set_yticks(y); ax.set_yticklabels([r["symbol"] for r in rows]); ax.invert_yaxis()
    ax.set_xlim(0, 60); ax.set_xlabel("donors whose tumour T cells express the gene (>= 10 of 100 cells)")
    pol = next(r for r in rows if r["status"] == "OPEN")
    ax.set_title("Panel A (July hits): 13 of 15 are not expressed in enough donors' T cells to test")
    fig.text(0.01, -0.06, f"{pol['symbol']}, the one testable gene: delete median {pol['del_median']:+.4f} (p = {pol['del_p']:.2f}), "
             f"overexpress {pol['ovx_median']:+.4f} (p = {pol['ovx_p']:.2f}), n = {pol['del_n']}, Holm m = 15 -> OPEN. "
             "NOT_RUN is not non-replication.", fontsize=7.5, color=NAVY)
    fig.savefig(os.path.join(out, "fig3_panel_a_testability.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 4
def fig_concordance(d, out):
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    cloud = [(v["delete"]["median"], v["overexpress"]["median"]) for v in d["fig"]["control_cloud"].values()
             if v["delete"]["median"] is not None and v["overexpress"]["median"] is not None]
    cx, cy = zip(*cloud)
    ax.scatter(cx, cy, s=8, color="#BDBDBD", alpha=0.7, label=f"matched control genes (n = {len(cloud)}), null cloud")
    tested = [r for r in d["rows"] if "del_median" in r]
    for r in tested:
        ax.scatter(r["del_median"], r["ovx_median"], s=34, color=STATUS_COLOUR[r["status"]], edgecolor="white", lw=0.5, zorder=3,
                   marker="D" if r["panel"] == "A" else "o")
        if r["status"] != "OPEN":
            ax.annotate(r["symbol"], (r["del_median"], r["ovx_median"]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.axhline(0, color=MUTED, lw=0.6); ax.axvline(0, color=MUTED, lw=0.6)
    allx = list(cx) + [r["del_median"] for r in tested]; ally = list(cy) + [r["ovx_median"] for r in tested]
    px, py = (max(allx) - min(allx)) * 0.05, (max(ally) - min(ally)) * 0.05
    ax.set_xlim(min(allx) - px, max(allx) + px); ax.set_ylim(min(ally) - py, max(ally) + 3 * py)
    ax.text(0.98, 0.02, "toward-normal quadrant\n(delete +, overexpress -)", transform=ax.transAxes, ha="right", fontsize=7.5, color=STATUS_COLOUR["T_CELL_SIGNAL_TOWARD"])
    ax.text(0.02, 0.02, "away-from-normal quadrant\n(delete -, overexpress +) is top-left", transform=ax.transAxes, fontsize=7.5, color=STATUS_COLOUR["T_CELL_SIGNAL_AWAY"])
    ax.set_xlabel("delete: median control-adjusted shift toward donor's own normal")
    ax.set_ylabel("overexpress: median control-adjusted shift")
    ax.set_title("Dose concordance vs the control-gene null cloud\n(both arms on the same token-positive cells, 3h)")
    legend_for(ax, ["T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY", "DELETION_ONLY", "OPEN"], loc="upper right")
    fig.savefig(os.path.join(out, "fig4_dose_concordance.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 5
def fig_per_donor(d, out):
    sig = [r for r in d["rows"] if "del_sig" in r and (r["del_sig"] or r["ovx_sig"])]
    sig.sort(key=lambda r: (list(STATUS_COLOUR).index(r["status"]), r["symbol"]))
    n = len(sig); cols = 6; rws = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rws, cols, figsize=(12, 2.0 * rws), sharey=False)
    for ax, r in zip(axes.flat, sig):
        for i, op in enumerate(("delete", "overexpress")):
            v = np.array(list(d["per_donor"][f"{r['gene']}|{op}"].values()))
            ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.18, 0.18, len(v)), v, s=6,
                       color=STATUS_COLOUR[r["status"]], alpha=0.7)
            ax.hlines(np.median(v), i - 0.3, i + 0.3, color=NAVY, lw=1.6)
        ax.axhline(0, color=MUTED, lw=0.5); ax.set_xticks([0, 1]); ax.set_xticklabels(["del", "ovx"], fontsize=7)
        ax.set_title(f"{r['symbol']} (n={r['del_n']})", fontsize=8, color=STATUS_COLOUR[r["status"]]); ax.tick_params(labelsize=6.5)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.suptitle(f"Per-donor control-adjusted shifts for every gene with a Holm-significant arm ({n} genes); bar = median",
                 color=NAVY, fontsize=10.5, y=1.01)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig5_per_donor.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 6
def fig_panel_b(d, out):
    rows = [r for r in d["rows"] if r["panel"] == "B"]
    order = {g["ensembl_id"]: i for i, g in enumerate(d["pb"]["genes"])}
    rows.sort(key=lambda r: order[r["gene"]])
    cannot = set(d["ann"]["cannot_attain"])
    fig, axes = plt.subplots(1, 2, figsize=(10, 9.5), sharey=True)
    y = np.arange(len(rows))
    for ax, op, title in ((axes[0], "del", "delete"), (axes[1], "ovx", "overexpress")):
        for yi, r in zip(y, rows):
            if f"{op}_median" not in r:
                ax.text(0, yi, r["status"].replace("_", " ").lower(), va="center", ha="center", fontsize=7, color=MUTED)
                continue
            c = STATUS_COLOUR[r["status"]]
            ax.hlines(yi, r[f"{op}_ci_lo"], r[f"{op}_ci_hi"], color=c, lw=1.4)
            ax.scatter(r[f"{op}_median"], yi, s=26 if r[f"{op}_sig"] else 14, color=c if r[f"{op}_sig"] else "white", edgecolor=c, zorder=3)
        ax.axvline(0, color=MUTED, lw=0.6); ax.set_title(f"{title}: median a, exact 95% CI")
        ax.set_xlabel("control-adjusted shift toward donor's own normal")
    labels = []
    for r in rows:
        s = r["symbol"]
        if r["symbol"] in cannot:
            s += "  (cannot attain)"
        labels.append(s)
    axes[0].set_yticks(y); axes[0].set_yticklabels(labels, fontsize=7.5); axes[0].invert_yaxis()
    for t, r in zip(axes[0].get_yticklabels(), rows):
        t.set_color(STATUS_COLOUR[r["status"]] if r["status"] not in ("OPEN", "NOT_RUN") else NAVY)
    legend_for(axes[1], ["T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY", "DELETION_ONLY", "OPEN", "NOT_RUN"],
               loc="upper center", bbox_to_anchor=(-0.1, -0.06), ncol=5)
    fig.text(0.5, -0.05, "filled marker = Holm-adjusted p <= 0.05 (m = 36); '(cannot attain)' = no effect size could reach significance at that n",
             ha="center", fontsize=7.5, color=MUTED)
    fig.suptitle("Panel B (39 curated T-cell genes): 11 toward, 2 away, 2 deletion-only, 19 open, 5 not run",
                 color=NAVY, fontsize=10.5, y=1.0)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig6_panel_b_forest.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 7
def fig_sensitivities(d, out):
    rows = {r["gene"]: r for r in d["rows"]}
    tested = [g for g, r in rows.items() if "del_median" in r]
    s = d["sens"]
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.5))
    # S1
    x = [rows[g]["del_median"] for g in tested]; yv = [s["S1_global_goal"][g]["del_median"] for g in tested]
    ax[0].scatter(x, yv, s=14, c=[STATUS_COLOUR[rows[g]["status"]] for g in tested])
    m = max(map(abs, x + yv)) * 1.1; ax[0].plot([-m, m], [-m, m], color=MUTED, lw=0.6, ls=":")
    ch = [rows[g]["symbol"] for g in tested if s["S1_global_goal"][g]["status"] != rows[g]["status"]]
    ax[0].set_title(f"S1 own-donor vs global goal\nwould change: {', '.join(ch) or 'none'}", fontsize=9)
    ax[0].set_xlabel("delete median, own goal"); ax[0].set_ylabel("delete median, global goal")
    # S2
    for g in tested:
        v = s["S2_leader_merad_vs_rest"][f"{g}|delete"]
        if v["leader_merad_median"] is not None and v["others_median"] is not None:
            ax[1].scatter(v["others_median"], v["leader_merad_median"], s=14, color=STATUS_COLOUR[rows[g]["status"]])
    ax[1].axhline(0, color=MUTED, lw=0.5); ax[1].axvline(0, color=MUTED, lw=0.5)
    ax[1].set_title("S2 Leader_Merad (23) vs other studies (20)\ndelete medians", fontsize=9)
    ax[1].set_xlabel("other studies"); ax[1].set_ylabel("Leader_Merad")
    # S3
    coh = [g for g in tested if rows[g]["status"] not in ("OPEN",)]
    mat = np.array([[rows[g]["fold_medians_delete"].get(str(k), np.nan) for k in range(5)] for g in coh])
    lim = np.nanmax(np.abs(mat))
    im = ax[2].imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax[2].set_yticks(range(len(coh))); ax[2].set_yticklabels([rows[g]["symbol"] for g in coh], fontsize=6.5)
    ax[2].set_xticks(range(5)); ax[2].set_xticklabels([f"f{k}" for k in range(5)])
    ax[2].set_title("S3 per-fold delete medians (called genes)\ndowngrades: none", fontsize=9)
    fig.colorbar(im, ax=ax[2], fraction=0.05)
    # S4
    for g in tested:
        v = s["S4_cell_weighted"][f"{g}|delete"]
        ax[3].scatter(v["donor_median_a"], v["cell_weighted_mean_shift"], s=14, color=STATUS_COLOUR[rows[g]["status"]])
    ax[3].axhline(0, color=MUTED, lw=0.5); ax[3].axvline(0, color=MUTED, lw=0.5)
    ax[3].set_title("S4 donor-weighted (adjusted) vs\ncell-weighted (raw) delete", fontsize=9)
    ax[3].set_xlabel("donor median a"); ax[3].set_ylabel("cell-weighted raw mean")
    fig.suptitle("Sensitivities: reported alongside, never used to change a registered status", color=NAVY, fontsize=10.5, y=1.04)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig7_sensitivities.png")); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 8
def fig_table(d, out):
    rows = d["rows"]; cannot = set(d["ann"]["cannot_attain"])
    fr = {}
    for name, v in d["sens"]["A3g_alternative_rules"].items():
        for g in v["status_changes"]:
            fr.setdefault(g, []).append("3g " + name.split("_")[0])
    for g in d["sens"]["A3h_B_all_cells_overexpress"]["status_changes"]:
        fr.setdefault(g, []).append("B all-cells")
    hdr = ["panel", "gene", "status", "concordance", "ambient", "estimable", "n del/ovx", "Holm del", "Holm ovx", "note"]
    cells = []
    for r in rows:
        est = estimable(d, r)
        note = []
        if r["symbol"] in cannot: note.append("cannot attain")
        if r["gene"] in fr: note.append("fragile: " + ", ".join(fr[r["gene"]]))
        if r.get("reason"): note.append(r["reason"][:48])
        f = lambda k: f"{r[k]:.2g}" if isinstance(r.get(k), (int, float)) and not isinstance(r.get(k), bool) else "—"
        cells.append([r["panel"], r["symbol"], r["status"], r.get("concordance") or "—", (r.get("ambient") or "—").replace("CURATED_", "").lower(),
                      str(est if est is not None else "—"), f"{r.get('del_n', '—')}/{r.get('ovx_n', '—')}", f("del_p_holm"), f("ovx_p_holm"), "; ".join(note)])
    fig, ax = plt.subplots(figsize=(13, 0.22 * len(cells) + 0.8)); ax.axis("off")
    t = ax.table(cellText=cells, colLabels=hdr, loc="center", cellLoc="left", colLoc="left",
                 colWidths=[0.04, 0.07, 0.17, 0.08, 0.1, 0.06, 0.07, 0.06, 0.06, 0.29])
    t.auto_set_font_size(False); t.set_fontsize(6.8); t.scale(1, 1.18)
    for (i, j), c in t.get_celld().items():
        c.set_edgecolor("#E3E8EF")
        if i == 0:
            c.set_facecolor("#E3E8EF"); c.set_text_props(weight="bold", color=NAVY)
        elif j == 2:
            c.set_text_props(color=STATUS_COLOUR[cells[i - 1][2]] if cells[i - 1][2] not in ("NOT_RUN", "NOT_ESTIMABLE_CONTROLS") else MUTED, weight="bold")
    ax.set_title(f"All {len(cells)} registered outcome rows (15 Panel A + 39 Panel B); host_drift = false on every row", loc="left", pad=18)
    fig.savefig(os.path.join(out, "fig8_outcome_table.png")); plt.close(fig)
    return len(cells)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    d = load(a.base)
    check_rows(d)
    os.makedirs(a.out, exist_ok=True)
    for f in (fig_design, fig_gate, fig_panel_a, fig_concordance, fig_per_donor, fig_panel_b, fig_sensitivities):
        f(d, a.out)
    n = fig_table(d, a.out)
    assert n == 54
    print("8 figures;", n, "rows in the outcome table")
    return 0


if __name__ == "__main__":
    sys.exit(main())
