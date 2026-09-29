"""Slide deck for the balanced-donor LUAD ISP results, in the house format of talk/JSDP_P25_talk.pptx
(16:9, #EEF1F5 background, Cambria navy titles, teal small-caps section label, Calibri body, white cards).
Numbers are read from committed Phase 7 outputs; each result slide names its confound on the slide.
Usage: python make_slides.py --base balanced_donor_luad --out balanced_donor_luad/slides/balanced_donor_luad_isp.pptx
"""
import argparse
import json
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Pt

BG, NAVY, BODY, MUTED, TEAL, RED, GREEN, PURPLE = "EEF1F5", "10243C", "16202C", "4A5768", "0B6F6A", "B8465E", "1A6B3C", "3D2F6B"
CARD, CARD_TEAL, CARD_RED, CARD_PURPLE, EDGE = "FFFFFF", "E8F2EC", "F9ECEF", "F0EDF9", "DCE3EC"
W, H = 12192000, 6858000
M = 548640                                                   # left margin used by the template


def rgb(h):
    return RGBColor.from_string(h)


class Deck:
    def __init__(self):
        self.p = Presentation(); self.p.slide_width, self.p.slide_height = W, H
        self.n = 0

    def slide(self, section, title, footer):
        s = self.p.slides.add_slide(self.p.slide_layouts[6]); self.n += 1
        s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(BG)
        if section:
            self.text(s, section.upper(), M, 91440, W - 2 * M, 201168, size=9.5, bold=True, color=TEAL, spacing=200)
        if title:
            self.text(s, title, M, 329184, W - 2 * M, 760000, size=26, bold=True, color=NAVY, font="Cambria")
        if footer:
            self.text(s, footer, M, 6419088, 7315200, 256032, size=9, color=MUTED)
            self.text(s, str(self.n), 11002975, 6419088, 640080, 256032, size=9, color=MUTED, align=PP_ALIGN.RIGHT)
        return s

    def text(self, s, txt, x, y, w, h, size=12.5, bold=False, color=BODY, font="Calibri", align=None, spacing=None, italic=False):
        tb = s.shapes.add_textbox(Emu(x), Emu(y), Emu(w), Emu(h)); tf = tb.text_frame; tf.word_wrap = True
        tf.margin_left = tf.margin_right = 0
        lines = txt if isinstance(txt, list) else [txt]
        for i, line in enumerate(lines):
            para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            if align is not None:
                para.alignment = align
            parts = line if isinstance(line, list) else [(line, bold)]
            for seg, b in parts:
                r = para.add_run(); r.text = seg; f = r.font
                f.size, f.bold, f.name, f.italic = Pt(size), b, font, italic; f.color.rgb = rgb(color)
                if spacing:
                    r._r.get_or_add_rPr().set("spc", str(spacing))
            para.space_after = Pt(4)
        return tb

    def card(self, s, x, y, w, h, fill=CARD, edge=EDGE, head=None, head_color=TEAL, body=None, size=11.5):
        sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Emu(x), Emu(y), Emu(w), Emu(h))
        sh.adjustments[0] = 0.04
        sh.fill.solid(); sh.fill.fore_color.rgb = rgb(fill); sh.line.color.rgb = rgb(edge); sh.line.width = Pt(0.75)
        sh.shadow.inherit = False
        pad = 228600; yy = y + 160000
        if head:
            self.text(s, head, x + pad, yy, w - 2 * pad, 300000, size=12.5, bold=True, color=head_color); yy += 330000
        if body:
            self.text(s, body, x + pad, yy, w - 2 * pad, h - (yy - y) - 100000, size=size, color=MUTED)
        return sh

    def image(self, s, path, x, y, w=None, h=None):
        """Fit inside the w x h box (either may be None), keeping the aspect ratio."""
        from PIL import Image
        iw, ih = Image.open(path).size
        sc = min((w / iw) if w else float("inf"), (h / ih) if h else float("inf"))
        return s.shapes.add_picture(path, Emu(x), Emu(y), width=Emu(int(iw * sc)), height=Emu(int(ih * sc)))

    def save(self, out):
        os.makedirs(os.path.dirname(out), exist_ok=True); self.p.save(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    j = lambda p: json.load(open(os.path.join(a.base, p)))
    rows = j("phase7_results/outcome_rows.json"); gate = j("phase4_results/classifier_gate.json")
    sec = j("phase7_results/secondaries.json"); ann = j("phase7_results/row_annotations.json")
    fig = lambda n: os.path.join(a.base, "figures", n)
    by = {}
    for r in rows:
        by.setdefault((r["panel"], r["status"]), []).append(r["symbol"])
    toward, away = by[("B", "T_CELL_SIGNAL_TOWARD")], by[("B", "T_CELL_SIGNAL_AWAY")]
    assert len(rows) == 54
    d = Deck()
    CW = W - 2 * M

    # 1 title
    s = d.slide(None, None, None)
    d.text(s, "Paired tumour-vs-normal T cells from 43 patients: which genes does a foundation model read as T-cell state?",
           M, 1000000, CW, 1400000, size=30, bold=True, color=NAVY, font="Cambria")
    d.text(s, [[("Kaisar Dauyey", True), ("  ·  ", False), ("Shinji Nakaoka", True)]], M, 2750000, CW, 300000, size=14, color=NAVY)
    d.text(s, "Laboratory of Mathematical Biology, Faculty of Advanced Life Science, Hokkaido University", M, 3080000, CW, 300000, size=11, color=MUTED)
    bw = (CW - 2 * 380000) // 3
    for i, (t, c) in enumerate((("43 donors, both tissues,\n100 cells each", TEAL), ("Delete AND overexpress,\ndonor's own normal as goal", PURPLE),
                                ("Registered before any number;\nall 54 rows reported", RED))):
        d.card(s, M + i * (bw + 380000), 3800000, bw, 1000000, edge=c, body=t.split("\n"), size=12)
    d.text(s, "Geneformer V2-316M (bf16) · in silico perturbation is a hypothesis generator, not a knockout screen",
           M, 5150000, CW, 300000, size=12, bold=True, color=RED, italic=True)

    # 2 question
    s = d.slide("The question", "Two registered questions, one cleaner design", "Question")
    d.card(s, M, 1300000, CW // 2 - 150000, 2300000, fill=CARD_TEAL, head="Panel A — do the July hits replicate?",
           body=["The 15 LUAD→normal top genes of the July screen.", "July compared LUAD with normal tissue from other studies (study-confounded)."], size=12)
    d.card(s, M + CW // 2 + 150000, 1300000, CW // 2 - 150000, 2300000, fill=CARD_PURPLE, head="Panel B — does the model read T-cell genes?",
           head_color=PURPLE, body=["The lab's 39 curated T-cell genes (36 perturbable).", "July never ranked them in its top 120."], size=12)
    d.card(s, M, 3850000, CW, 1900000, head="What 'a hit' means here", body=[
        "Per donor: shift of tumour T cells toward the SAME donor's normal T cells, minus the median of 20 matched control genes.",
        "Exact Wilcoxon over donors, Holm per panel, for deletion and overexpression separately.",
        "A signal needs BOTH arms significant with opposite signs (dose concordance)."], size=12)

    # 3 design
    s = d.slide("Design", "Pairing within a donor cancels study, chemistry and person — not ambient RNA", "Balanced paired design")
    d.image(s, fig("fig1_design.png"), M, 1250000, w=CW, h=3550000)
    d.card(s, M, 4950000, CW, 1150000, fill=CARD_RED, head="What pairing cannot fix", head_color=RED, body=[
        "Tumour and normal tissue differ in ambient contamination within each donor. 20 of 43 donors overlap July's LUAD pool. Selection favours T-cell-rich samples; 17 unpaired donors excluded."], size=11.5)

    # 4 gate
    s = d.slide("Classifier gate", f"The model tells a donor's tumour from normal T cells: {gate['donors_above_0.5']}/{gate['n_donors']} donors above chance",
                "Held-out, 5-fold donor cross-fitting")
    d.image(s, fig("fig2_classifier_gate.png"), M, 1300000, w=int(CW * 0.62))
    d.card(s, M + int(CW * 0.65), 1300000, int(CW * 0.35), 2500000, head="Gate: PASS", body=[
        f"Pooled balanced accuracy {gate['pooled_balanced_accuracy']:.3f} (gate 0.60).",
        f"Sign test p = {gate['sign_test_p_float']:.1e}.",
        "Every donor scored by a model that never saw it.",
        "Fine-tune not bit-reproducible (band 0.045 per donor); cannot flip the gate."], size=11.5)

    # 5 panel A
    s = d.slide("Panel A — July hits", "13 of 15 July hits are not expressed in T cells — not testable, not refuted",
                "NOT_RUN is not non-replication")
    d.image(s, fig("fig3_panel_a_testability.png"), M, 1250000, w=int(CW * 0.64), h=4800000)
    d.card(s, M + int(CW * 0.67), 1250000, int(CW * 0.33), 3300000, head="Reading", body=[
        "Keratins, mucin, a secretoglobin, a haemoglobin chain: epithelial and blood markers, present in ≤ 8 of 43 donors' T cells.",
        "The one testable gene, POLR2J3, is OPEN (p = 0.10 / 0.25).",
        "The design did not refute the July panel; it showed the panel is mostly not T-cell biology.",
        "Confound: 316M replaced 104M (no 104M arm), so a non-replication could not be attributed — moot, none occurred."], size=11)

    # 6 panel B
    s = d.slide("Panel B — T-cell genes", f"{len(toward) + len(away)} of 34 tested T-cell genes give a dose-concordant signal",
                "Curated T-cell genes")
    d.image(s, fig("fig6_panel_b_forest.png"), M, 1150000, w=5600000, h=5150000)
    x0 = M + 5800000; cw = W - M - x0
    d.card(s, x0, 1250000, cw, 1700000, fill=CARD_TEAL, head=f"Toward normal ({len(toward)})", body=[", ".join(toward)], size=11.5)
    d.card(s, x0, 3050000, cw, 800000, head=f"Away from normal ({len(away)})", head_color=RED, body=[", ".join(away)], size=11.5)
    d.card(s, x0, 3950000, cw, 2050000, head="Stable — and what cannot be judged", body=[
        "None of the 13 moves under any registered sensitivity.",
        f"Below the design d_min of 11 (registered): {', '.join(ann['below_registered_d_min'])}. Of these, {', '.join(ann['cannot_attain'])} cannot attain significance at any effect size in this family (post hoc) — no evidence either way.",
        "TRAC/TRBC1/TRBC2: not in the model's vocabulary."], size=11)

    # 7 concordance
    s = d.slide("Dose concordance", "Opposite-signed arms are common even for control genes — effect size separates few",
                "Both arms on the same token-positive cells")
    d.image(s, fig("fig4_dose_concordance.png"), M, 1200000, w=5900000, h=5000000)
    x0 = M + 6100000; cw = W - M - x0
    d.card(s, x0, 1300000, cw, 2300000, fill=CARD_RED, head="Confound", head_color=RED, body=[
        "360 matched control entries are anti-correlated too (ρ = −0.62); 35% sit in the toward-normal quadrant.",
        "Post hoc, descriptive: only CCR7, CD3D, GZMA, CD7 (and CD27 on overexpression) lie beyond ~95% of controls."], size=11.5)
    d.card(s, x0, 3750000, cw, 2200000, head="Measured on the same cells (Amendment 3h)", body=[
        "Geneformer overexpresses into all cells; the token-positive subset was reconstructed and checked.",
        "15,142 / 15,142 count checks; 8 / 8 GPU identity pairs bit-identical."], size=11.5)

    # 8 per donor
    s = d.slide("Per donor", "The medians are not carried by a few donors", "Per-donor control-adjusted shifts")
    d.image(s, fig("fig5_per_donor.png"), M, 1150000, w=CW, h=5150000)

    # 9 sensitivities
    s = d.slide("Sensitivities", "The coherent set does not move under any registered sensitivity",
                "Reported alongside; never used to change a status")
    d.image(s, fig("fig7_sensitivities.png"), M, 1250000, w=CW, h=3100000)
    d.card(s, M, 4450000, CW, 1500000, head="Three fragile rows, outside the coherent set", body=[
        "GZMK: deletion-only → away under all-cells overexpress.  PRF1: open → deletion-only if any control counts.  CD69: deletion-only → dose-incoherent under Holm-over-tested and under the global goal.",
        f"Panel-level: {sec['B_all']['n_toward']} toward vs {sec['B_all']['n_away']} away, sign test p = {sec['B_all']['sign_test_p']:.2f} — genes share controls, so signs are not independent."], size=11)

    # 10 interpretation
    s = d.slide("Interpretation", "What each status licenses — and what it does not", "Interpretation breakdown")
    items = [("T-cell signal (13)", "Model-internal, donor-consistent, dose-concordant response beyond matched controls.", "Not a causal driver; mostly not unusual among genes.", TEAL),
             ("Deletion only (2)", "A consistent deletion effect.", "No direction; no signal claim.", "E69F00"),
             ("Open (20)", "The registered test did not call it.", "Not 'no effect'; for CD28/CTLA4/PDCD1 no evidence either way.", MUTED),
             ("Not run (18) · no controls (1)", "Not expressed in enough donors' T cells, or not in the vocabulary.", "Not non-replication; not no effect.", "93A3B5")]
    cw = (CW - 3 * 200000) // 4
    for i, (h, yes, no, c) in enumerate(items):
        d.card(s, M + i * (cw + 200000), 1350000, cw, 2900000, edge=c, head=h, head_color=NAVY if c == MUTED else c,
               body=["Licenses: " + yes, "", "Does not: " + no], size=11.5)

    # 11 limitations
    s = d.slide("Limitations", "What this design still cannot rule out", "Limitations")
    lim = [("Ambient RNA", "Tumour vs normal tissue differ in contamination within a donor; pairing cannot remove it."),
           ("Model change", "316M bf16 replaced July's 104M fp32; no 104M arm."),
           ("Not independent", "20 of 43 donors overlap July's LUAD pool."),
           ("Selection", "T-cell-rich samples only; 17 unpaired donors excluded."),
           ("Nondeterminism", "Fine-tuning is not bit-reproducible; the ISP inherits these fold models."),
           ("Process", "No-op and wrapper checks ran after Phase 6 started; the overexpress cell set was fixed by Amendment 3h after Phase 6.")]
    cw = (CW - 2 * 200000) // 3
    for i, (h, b) in enumerate(lim):
        d.card(s, M + (i % 3) * (cw + 200000), 1350000 + (i // 3) * 1700000, cw, 1500000, head=h, body=[b], size=11.5)

    # 12 July + next
    s = d.slide("So what", "For the July screen: its LUAD→normal list was mostly not T-cell biology", "Implications and next steps")
    d.card(s, M, 1350000, CW // 2 - 150000, 2800000, fill=CARD_TEAL, head="What this changes", body=[
        "13 of 15 July hits cannot be tested in T cells under a paired, balanced design.",
        "The model does read canonical T-cell genes (CD3 complex, LCK, LAT, CD8A, CCR7, GZMA) in a donor-consistent, dose-concordant way.",
        "Effect size, not significance, separates a handful from arbitrary genes."], size=12)
    d.card(s, M + CW // 2 + 150000, 1350000, CW // 2 - 150000, 2800000, head="Next", body=[
        "Genome-wide control distribution, to rank effect size formally (registered first).",
        "Checkpoint genes (CTLA4, PDCD1, CD28) need more token-positive donors: a larger cohort or a looser cap.",
        "A 104M arm, if attribution to model vs design is wanted.",
        "Independent tissue validation for the top four (CCR7, CD3D, GZMA, CD7)."], size=12)

    d.save(a.out)
    print(d.n, "slides ->", a.out)


if __name__ == "__main__":
    main()
