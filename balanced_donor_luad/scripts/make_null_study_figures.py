"""Null study figures for the IMRaD report (Amendments 4-6), from the two files the gatekeeper signed.

Reads ONLY --a5-result (the A5 N=100 primary output) and --n200-result (the N=200 combined output),
pinned by sha256 to the exact hashes Stanley PASSED (2026-09-30T11-06-43Z). Refuses (exit 3) if either
file does not match, before making any plot -- these figures must never be drawn from a different run
than the one the report's numbers already cite. No new statistics are computed from raw data; every
number plotted is read directly from the signed files' own fields.

The one exception is figure (b), the reference null distribution: the signed files do not store the
100,000 individual permutation draws, only the observed rho and that p sits at the permutation floor
(no draw reached it). Figure (b) therefore draws an illustrative reference distribution of Spearman rho
under independence at the same n, from a fresh simulation of independent random data -- never from the
study's own gene data -- exactly the same construction the registration's own power-by-simulation
methodology (Amendment 6 s.6.4) already uses, and the caption says so explicitly.

Captions follow Stanley's 2026-09-30T11-06-43Z wording rules: "random eligible genes detectable in LUAD
tumour T cells", never "genome-wide"; p reported as "at the permutation floor", never as an exact
value; N=200 as "confirming" the N=100 primary, never "replicated"; the R3 validity comparison to -0.62
by sign only.
"""
import argparse
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EXPECTED = {
    "a5": "a00ccb2f62a7a888386208eb63a2f349d1c2ee8717208f742a563bc508588267",
    "n200": "c9e56026fac99d30a89da91ab557b09e7889717163712242a961d33295a3b847",
}
SEED = 20260929  # the study's own registered seed, reused here only for the illustrative reference curve


def load_pinned(path, key):
    raw = open(path, "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    if h != EXPECTED[key]:
        raise SystemExit(f"REFUSING: {path} sha256 {h} != signed {EXPECTED[key]}")
    return json.loads(raw)


def fig_scatter(a5, n200, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, rows, n, title in (
        (axes[0], [r for r in n200["rows"][:100] if r["status"] == "estimable"], 100, "N=100 (registered primary)"),
        (axes[1], [r for r in n200["rows"] if r["status"] == "estimable"], 200, "N=200 (nested confirmatory extension)"),
    ):
        x = np.array([r["delete_median"] for r in rows])
        y = np.array([r["overexpress_median"] for r in rows])
        ax.scatter(x, y, s=16, alpha=0.6, color="#0B6F6A", edgecolor="none")
        ax.axhline(0, color="0.8", lw=0.8, zorder=0)
        ax.axvline(0, color="0.8", lw=0.8, zorder=0)
        # symlog: one gene's shift is roughly an order of magnitude larger than the rest on this axis
        # (real data, not an artifact -- see the report text); linear axes would compress every other
        # point into an unreadable cluster.
        ax.set_xscale("symlog", linthresh=0.01)
        ax.set_yscale("symlog", linthresh=0.01)
        rho = a5["primary"]["rho"] if n == 100 else n200["extension_n200"]["rho"]
        ax.set_title(f"{title}\nSpearman rho = {rho:.3f}, n = {len(x)}", fontsize=10)
        ax.set_xlabel("raw donor-median deletion shift (symlog)")
        ax.set_ylabel("raw donor-median overexpression shift (symlog)")
    fig.suptitle("Delete vs. overexpress shift, random eligible genes detectable in LUAD tumour T cells",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_null_reference(a5, n200, out_path):
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, n, rho_obs, label in (
        (axes[0], 100, a5["primary"]["rho"], "N=100 primary"),
        (axes[1], 200, n200["extension_n200"]["rho"], "N=200 extension"),
    ):
        # illustrative reference only: independent random data at the same n, not the study's genes
        sims = np.array([
            np.corrcoef(np.argsort(np.argsort(rng.standard_normal(n))),
                        np.argsort(np.argsort(rng.standard_normal(n))))[0, 1]
            for _ in range(20000)
        ])
        ax.hist(sims, bins=60, color="0.75", edgecolor="none")
        ax.axvline(rho_obs, color="#8B1E3F", lw=2)
        ax.text(rho_obs, ax.get_ylim()[1] * 0.92, f"  observed rho = {rho_obs:.3f}\n  p at the permutation\n  floor (0/100,000 draws)",
                color="#8B1E3F", fontsize=8, va="top")
        ax.set_title(f"{label}: reference null (independent data, illustrative)", fontsize=9.5)
        ax.set_xlabel("Spearman rho")
        ax.set_ylabel("count (20,000 simulated draws)")
    fig.suptitle("Reference null distribution of Spearman rho under independence -- illustrative only;\n"
                 "the registered test permuted the study's own 100,000 times and found no draw at or beyond "
                 "the observed value", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_loo_stability(a5, n200, out_path):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    labels = ["N=100 primary", "N=200 extension"]
    pass_n = [a5["primary"]["loo_n_pass"], n200["extension_n200"]["loo_n_pass"]]
    total_n = [a5["primary"]["loo_n"], n200["extension_n200"]["loo_n"]]
    fail_n = [t - p for t, p in zip(total_n, pass_n)]
    x = np.arange(len(labels))
    ax.bar(x, pass_n, color="#0B6F6A", label="leave-one-out recomputation stable", width=0.6)
    ax.bar(x, fail_n, bottom=pass_n, color="#8B1E3F", label="not stable", width=0.6)
    for i, (p, t) in enumerate(zip(pass_n, total_n)):
        ax.text(i, t + max(total_n) * 0.03, f"{p}/{t}", ha="center", fontsize=10)
    boot = [a5["primary"]["bootstrap_fraction_stable"], n200["extension_n200"]["bootstrap_fraction_stable"]]
    for i, b in enumerate(boot):
        ax.text(i, total_n[i] * 0.4, f"bootstrap stable\nfraction: {b:.3f}\n(10,000 draws)",
                ha="center", va="center", fontsize=8.5, color="white")
    ax.set_xticks(x, labels)
    ax.set_ylabel("genes")
    ax.set_title("Leave-one-out and bootstrap stability, both runs")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2)
    ax.set_ylim(0, max(total_n) * 1.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_amendment_flow(out_path):
    events = [
        ("Amendment 4", "622efaf", "200-gene draw registered"),
        ("cost corrections", "3d26d4a, f945374", "wall-clock -> per-call seconds; division error fixed"),
        ("estimable-count finding", "-", "only 12/200 drawn genes estimable"),
        ("Amendment 5", "3c4a5a4 .. 36476f7", "N=100 estimable redesign; 3 gate cycles, PASSED"),
        ("A5 launch -> complete", "-", "6.6153 GPU-h, no-op 5/5"),
        ("Amendment 6", "1b2f2fe, 037683c", "N=200 extension + s.6.7 pair rule, PASSED"),
        ("A6 run-package gates", "94fc5b6, c287178, 2fe4e31", "B1-B6, m1, m2, pip-freeze fix; 2 launch attempts"),
        ("A6 launch -> complete", "-", "6.6555 GPU-h, no-op 5/5"),
        ("V1 fix + PASS", "-", "merged 318-control tree; gatekeeper PASS 2026-09-30"),
    ]
    fig, ax = plt.subplots(figsize=(7, 9))
    ax.axis("off")
    n = len(events)
    for i, (title, commit, note) in enumerate(events):
        y = n - i
        ax.add_patch(plt.Rectangle((0.05, y - 0.38), 0.9, 0.7, fc="#EEF1F5", ec="#0B6F6A", lw=1.2))
        ax.text(0.1, y + 0.12, title, fontsize=10, weight="bold", va="center")
        ax.text(0.1, y - 0.16, f"{commit}  --  {note}", fontsize=7.5, va="center", color="0.25")
        if i < n - 1:
            ax.annotate("", xy=(0.5, y - 0.62), xytext=(0.5, y - 0.42),
                        arrowprops=dict(arrowstyle="-|>", color="#0B6F6A"))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n + 1)
    ax.set_title("Null study: amendments and gates, in order", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a5-result", required=True)
    ap.add_argument("--n200-result", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    a5 = load_pinned(a.a5_result, "a5")
    n200 = load_pinned(a.n200_result, "n200")
    fig_scatter(a5, n200, f"{a.out_dir}/fig9_null_scatter.png")
    fig_null_reference(a5, n200, f"{a.out_dir}/fig10_null_reference_distribution.png")
    fig_loo_stability(a5, n200, f"{a.out_dir}/fig11_null_loo_stability.png")
    fig_amendment_flow(f"{a.out_dir}/fig12_amendment_flow.png")
    print("wrote fig9-fig12 to", a.out_dir)


if __name__ == "__main__":
    main()
