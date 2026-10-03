#!/usr/bin/env python3
"""Build the lung-to-colon graphical abstract as an SVG scene inside one HTML page.

Every value is read from committed result files and checked against README.md
("Latest studies") before anything is drawn:

- balanced_donor_luad/phase4_results/classifier_gate.json and
  pelka_crc_e2/phase4_results/classifier_gate.json (per-donor balanced accuracy);
- pelka_crc_e2/results/h2b_null_result.json (100 null genes, deletion and
  overexpression medians, the H2b test);
- balanced_donor_luad/phase7_results/outcome_rows.json and
  pelka_crc_e2/results/panel_b/outcome_rows.json with h2c_result.json (Panel B).

The background is drawn from a fixed seed. Render to PNG with render.mjs.
"""

from __future__ import annotations

import json
import math
import random
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
README = re.sub(r"\s+", " ", (ROOT / "README.md").read_text())
W, H = 1600, 800
rng = random.Random(20261003)


def load(rel: str):
    return json.loads((ROOT / rel).read_text())


def in_readme(text: str) -> None:
    assert text in README, f"not in README: {text!r}"


def ranks(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    for k, i in enumerate(order):
        r[i] = k
    return r


def spearman(a: list[float], b: list[float]) -> float:
    ra, rb = ranks(a), ranks(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    return cov / math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))


# ------------------------------------------------------------------ numbers
gate = {"lung": load("balanced_donor_luad/phase4_results/classifier_gate.json"),
        "colon": load("pelka_crc_e2/phase4_results/classifier_gate.json")}
for g in gate.values():
    assert g["donors_above_0.5"] == g["n_donors"]
in_readme(f"pooled held-out balanced accuracy {gate['lung']['pooled_balanced_accuracy']:.3f}, 43 of 43 donors")
in_readme(f"balanced accuracy {gate['colon']['pooled_balanced_accuracy']:.3f}, 19 of 19 donors")

h2b = load("pelka_crc_e2/results/h2b_null_result.json")
prim = h2b["primary"]
dele, ovx = h2b["null_delete_medians"], h2b["null_overexpress_medians"]
assert len(dele) == len(ovx) == prim["n_estimable"] == 100
assert abs(spearman(dele, ovx) - prim["rho"]) < 1e-9
in_readme(f"rho {prim['rho']:.3f} (p = {prim['p_lower_tail']:.4f})".replace("-", "\u2212"))
in_readme(f"only {round(100 * prim['bootstrap_fraction_stable'])}% of gene bootstraps")
in_readme(f"`{prim['status']}`")
in_readme("(rho −0.593)")

h2c = load("pelka_crc_e2/results/h2c_result.json")
luad = {x["symbol"]: x for x in load("balanced_donor_luad/phase7_results/outcome_rows.json")
        if x.get("panel") == "B" and x.get("del_median") is not None}
crc_rows = load("pelka_crc_e2/results/panel_b/outcome_rows.json")
crc = {x["symbol"]: x for x in crc_rows if x.get("del_median") is not None}
shared = sorted(set(luad) & set(crc))
ref = [g for g, x in luad.items() if x["status"] in ("T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY")]
tested = [g for g in ref if g in crc]
agree = [g for g in tested if (luad[g]["del_median"] > 0) == (crc[g]["del_median"] > 0)]
assert (len(ref), len(tested), len(agree)) == (h2c["reference_genes"], h2c["n_tested"], h2c["del_agree"])
in_readme(f"{len(agree)} of {len(tested)} LUAD reference genes keep their deletion sign (p = {h2c['p_one_sided']:.2f})")
crc_hits = [x["symbol"] for x in crc_rows if x["status"] in ("T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY")]
assert crc_hits == ["PRF1"]
rho_b = spearman([luad[g]["del_median"] for g in shared], [crc[g]["del_median"] for g in shared])

# ------------------------------------------------------------------ palette
INK, MUTED, GOLD, LUNG, COLON, GREY, CORAL = "#f1eef8", "#b3abc8", "#ffc857", "#7fb7ff", "#ff8fb1", "#8f8aa3", "#ff8a5b"


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=16, weight=400, fill=INK, anchor="start", extra=""):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-weight="{weight}" '
            f'fill="{fill}" text-anchor="{anchor}" {extra}>{esc(s)}</text>')


def backdrop() -> str:
    out = [f'<rect width="{W}" height="{H}" fill="url(#bg)"/>']
    for _ in range(70):
        x, y, r = rng.uniform(0, W), rng.uniform(0, H), rng.uniform(10, 34)
        col = rng.choice([LUNG, COLON, "#b39ddb"])
        out.append(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="{r:.0f}" fill="none" stroke="{col}" '
                   f'stroke-opacity="{rng.uniform(0.05, 0.13):.2f}" stroke-width="{rng.uniform(1, 3):.1f}"/>')
        out.append(f'<circle cx="{x + rng.uniform(-4, 4):.0f}" cy="{y + rng.uniform(-4, 4):.0f}" r="{r * 0.35:.0f}" '
                   f'fill="{col}" fill-opacity="{rng.uniform(0.03, 0.07):.2f}"/>')
    return "\n".join(out)


def frame(x0: float, num: str, kicker: str, title: str) -> str:
    return "\n".join([
        f'<rect x="{x0}" y="96" width="490" height="600" rx="26" fill="#0d0a18" fill-opacity="0.74" '
        f'stroke="#c9b8f0" stroke-opacity="0.22"/>',
        f'<circle cx="{x0 + 52}" cy="140" r="20" fill="none" stroke="{GOLD}" stroke-width="2"/>',
        text(x0 + 52, 147, num, 20, 800, GOLD, "middle"),
        text(x0 + 84, 135, kicker, 14, 800, GOLD, extra='letter-spacing="3"'),
        text(x0 + 84, 158, title, 19, 700),
    ])


def callout(x0, y, line1, line2, col=GOLD):
    return "\n".join([
        f'<rect x="{x0 + 40}" y="{y}" width="410" height="58" rx="14" fill="#170f2b" stroke="{col}" stroke-opacity="0.55"/>',
        text(x0 + 245, y + 25, line1, 16, 700, col, "middle"),
        text(x0 + 245, y + 46, line2, 13, 400, MUTED, "middle"),
    ])


# ------------------------------------------------------------------ panel 1
def panel_classifier(x0: float) -> str:
    out = [frame(x0, "1", "THE CLASSIFIER", "Tumour vs normal T cells, same donor")]
    lo, hi = 0.4, 1.0
    sy = lambda v: 560 - 330 * (v - lo) / (hi - lo)
    for v in (0.5, 0.6, 0.8, 1.0):
        dash = 'stroke-dasharray="4 5"' if v != 0.6 else ""
        out.append(f'<line x1="{x0 + 70}" y1="{sy(v):.1f}" x2="{x0 + 450}" y2="{sy(v):.1f}" stroke="{MUTED}" '
                   f'stroke-opacity="{0.55 if v in (0.5, 0.6) else 0.18}" {dash}/>')
        out.append(text(x0 + 62, sy(v) + 4, f"{v:.1f}", 12, 400, MUTED, "end"))
    out.append(text(x0 + 452, sy(0.6) - 6, "gate 0.60", 12, 600, MUTED, "end"))
    out.append(text(x0 + 452, sy(0.5) + 15, "chance", 12, 600, MUTED, "end"))
    for i, (key, col, label) in enumerate((("lung", LUNG, "Lung (LUAD)"), ("colon", COLON, "Colon (E2)"))):
        cx = x0 + 165 + i * 190
        vals = list(gate[key]["per_donor_balanced_accuracy"].values())
        for v in vals:
            out.append(f'<circle cx="{cx + rng.uniform(-38, 38):.1f}" cy="{sy(v):.1f}" r="5" fill="{col}" fill-opacity="0.7"/>')
        pooled = gate[key]["pooled_balanced_accuracy"]
        out.append(f'<line x1="{cx - 50}" y1="{sy(pooled):.1f}" x2="{cx + 50}" y2="{sy(pooled):.1f}" stroke="{GOLD}" stroke-width="4"/>')
        out.append(text(cx, 590, label, 16, 700, col, "middle"))
        out.append(text(cx, 610, f"{len(vals)} donors · pooled {pooled:.3f}", 13, 400, MUTED, "middle"))
    out.append(text(x0 + 30, sy(0.7), "balanced accuracy per held-out donor", 12, 600, MUTED, "middle",
                    f'transform="rotate(-90 {x0 + 30} {sy(0.7):.1f})"'))
    out.append(callout(x0, 624, "Transfers: every donor above chance", "fresh fold models, same recipe, no LUAD weights reused"))
    return "\n".join(out)


# ------------------------------------------------------------------ panel 2
def panel_null(x0: float) -> str:
    out = [frame(x0, "2", "RANDOM GENES", "Built-in anti-correlation, weaker in colon")]
    pad = 0.08
    xs, ys = dele, ovx
    xr = (min(xs), max(xs)); yr = (min(ys), max(ys))
    sx = lambda v: x0 + 80 + 350 * (v - xr[0]) / (xr[1] - xr[0])
    sy = lambda v: 540 - 320 * (v - yr[0]) / (yr[1] - yr[0])
    out.append(f'<rect x="{x0 + 70}" y="210" width="370" height="340" fill="none" stroke="{MUTED}" stroke-opacity="0.25"/>')
    if xr[0] < 0 < xr[1]:
        out.append(f'<line x1="{sx(0):.1f}" y1="210" x2="{sx(0):.1f}" y2="550" stroke="{MUTED}" stroke-opacity="0.35"/>')
    if yr[0] < 0 < yr[1]:
        out.append(f'<line x1="{x0 + 70}" y1="{sy(0):.1f}" x2="{x0 + 440}" y2="{sy(0):.1f}" stroke="{MUTED}" stroke-opacity="0.35"/>')
    for a, b in zip(xs, ys):
        out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="4.5" fill="{COLON}" fill-opacity="0.7"/>')
    # least-squares line on ranks, drawn in data space for orientation only
    n = len(xs); mx, my = sum(xs) / n, sum(ys) / n
    slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
    x1, x2 = xr
    out.append(f'<line x1="{sx(x1):.1f}" y1="{sy(my + slope * (x1 - mx)):.1f}" x2="{sx(x2):.1f}" '
               f'y2="{sy(my + slope * (x2 - mx)):.1f}" stroke="{GOLD}" stroke-width="2.5" stroke-dasharray="7 5"/>')
    out.append(text(x0 + 255, 572, "deletion shift toward normal (median over donors)", 12, 600, MUTED, "middle"))
    out.append(text(x0 + 56, 380, "overexpression shift", 12, 600, MUTED, "middle", f'transform="rotate(-90 {x0 + 56} 380)"'))
    out.append(text(x0 + 430, 232, f"100 random genes, colon", 13, 700, COLON, "end"))
    out.append(text(x0 + 40, 600, f"Colon rho {prim['rho']:.3f} (p = {prim['p_lower_tail']:.4f})".replace("-", "\u2212"), 17, 800, COLON))
    out.append(text(x0 + 450, 600, "Lung rho −0.593", 17, 800, LUNG, "end"))
    out.append(callout(x0, 624, "control_draw_sensitive_open",
                       f"same direction as lung; {round(100 * prim['bootstrap_fraction_stable'])}% of bootstraps stable vs a 95% bar"))
    return "\n".join(out)


# ------------------------------------------------------------------ panel 3
def panel_genes(x0: float) -> str:
    out = [frame(x0, "3", "CURATED T-CELL GENES", "The lung gene pattern does not replicate")]
    xs = [luad[g]["del_median"] for g in shared]
    ys = [crc[g]["del_median"] for g in shared]
    lim = max(abs(v) for v in xs + ys) * 1.08
    sx = lambda v: x0 + 255 + 175 * v / lim
    sy = lambda v: 380 - 165 * v / lim
    out.append(f'<rect x="{x0 + 80}" y="215" width="350" height="330" fill="none" stroke="{MUTED}" stroke-opacity="0.25"/>')
    out.append(f'<rect x="{sx(0):.1f}" y="215" width="{x0 + 430 - sx(0):.1f}" height="{sy(0) - 215:.1f}" fill="{LUNG}" fill-opacity="0.06"/>')
    out.append(f'<rect x="{x0 + 80}" y="{sy(0):.1f}" width="{sx(0) - x0 - 80:.1f}" height="{545 - sy(0):.1f}" fill="{LUNG}" fill-opacity="0.06"/>')
    out.append(f'<line x1="{sx(0):.1f}" y1="215" x2="{sx(0):.1f}" y2="545" stroke="{MUTED}" stroke-opacity="0.4"/>')
    out.append(f'<line x1="{x0 + 80}" y1="{sy(0):.1f}" x2="{x0 + 430}" y2="{sy(0):.1f}" stroke="{MUTED}" stroke-opacity="0.4"/>')
    pos = {g: (sx(luad[g]["del_median"]), sy(crc[g]["del_median"])) for g in tested}
    crowded = sorted(g for g in tested if any(h != g and math.hypot(pos[g][0] - pos[h][0], pos[g][1] - pos[h][1]) < 40 for h in tested))
    for g in shared:
        a, b = luad[g]["del_median"], crc[g]["del_median"]
        if g in tested:
            ok = g in agree
            out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="7" fill="{GOLD if ok else CORAL}" stroke="#0d0a18" stroke-width="1.5"/>')
            if g not in crowded:
                out.append(text(sx(a) + 10, sy(b) + 4, g, 11, 600, INK))
        else:
            fill = COLON if g in crc_hits else GREY
            out.append(f'<circle cx="{sx(a):.1f}" cy="{sy(b):.1f}" r="{6 if g in crc_hits else 4}" fill="{fill}" fill-opacity="0.85"/>')
            if g in crc_hits:
                out.append(text(sx(a) + 10, sy(b) + 4, f"{g} (colon hit)", 11, 700, COLON))
    out.append(text(x0 + 255, 565, "deletion median, lung", 12, 600, LUNG, "middle"))
    out.append(text(x0 + 66, 380, "deletion median, colon", 12, 600, COLON, "middle", f'transform="rotate(-90 {x0 + 66} 380)"'))
    out.append(text(x0 + 425, 233, "shaded: same sign", 11, 600, MUTED, "end"))
    if crowded:
        kept = [g for g in crowded if g in agree]; lost = [g for g in crowded if g not in agree]
        parts = ([f"kept: {', '.join(kept)}"] if kept else []) + ([f"lost: {', '.join(lost)}"] if lost else [])
        out.append(text(x0 + 255, 538, "unlabelled cluster · " + " · ".join(parts), 11, 600, MUTED, "middle"))
    lx, ly = x0 + 40, 598
    for j, (col, lab) in enumerate(((GOLD, "LUAD hit, sign kept"), (CORAL, "LUAD hit, sign lost"), (GREY, "other panel gene"))):
        cx = lx + 8 + j * 145
        out.append(f'<circle cx="{cx}" cy="{ly - 4}" r="6" fill="{col}"/>')
        out.append(text(cx + 12, ly, lab, 12, 600))
    out.append(callout(x0, 624, f"{len(agree)} of {len(tested)} lung hits keep their sign",
                       f"pattern_not_replicated · {len(shared)} shared genes, rho {rho_b:.2f}".replace("-", "\u2212").replace("pattern\u2212not\u2212replicated", "pattern_not_replicated")))
    return "\n".join(out)


DEFS = """
<defs>
  <linearGradient id="bg" x1="0" y1="0" x2="0.4" y2="1">
    <stop offset="0" stop-color="#2a1d4a"/><stop offset="0.55" stop-color="#160f2a"/><stop offset="1" stop-color="#0a0714"/>
  </linearGradient>
</defs>"""


def svg() -> str:
    return "\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        f'font-family="Inter, \'Liberation Sans\', \'DejaVu Sans\', sans-serif">', DEFS, backdrop(),
        text(40, 58, "Does a Geneformer T-cell perturbation screen transfer from lung to colon?", 30, 800),
        text(40, 84, "Balanced-donor design, Geneformer-V2-316M, re-run unchanged on Pelka 2021 colorectal cancer (E2)", 16, 400, MUTED),
        panel_classifier(40), panel_null(555), panel_genes(1070),
        f'<rect x="0" y="{H - 62}" width="{W}" height="62" fill="#07050e" fill-opacity="0.85"/>',
        text(W / 2, H - 25, "Pre-registered readings · donor cross-fitting · study, chemistry and annotation also "
             "differ between cohorts, so tissue is not the only explanation", 16, 600, MUTED, "middle"),
        "</svg>"])


def main() -> None:
    body = svg()
    (HERE / "graphical-abstract.svg").write_text(body + "\n")
    (HERE / "graphical-abstract.html").write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Lung to Colon Abstract</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700;800&display=block">
<style>html,body{{margin:0;background:#0a0714}} svg{{display:block;width:100vw;height:auto}}</style>
</head><body>
{body}
</body></html>
""")
    print("wrote graphical-abstract.svg and graphical-abstract.html")


if __name__ == "__main__":
    main()
