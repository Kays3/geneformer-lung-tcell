"""IMRaD-restructured oral-presentation deck, geneformer-lung-tcell, 2026-09-25.

18 slides across five explicit sections (Background / Introduction / Methods /
Results / Discussion), restoring the original JSDP_P25_talk / poster_final
results as the Results section's starting point, with the post-poster audit
findings (S100A8/S100A9 concordance failure, the ambient screen-level audit,
the T6 donor-level reanalysis, the bf16 precision control) following.

House standard, unchanged from 2026-09-24: every number is pulled from its
source file at build time via `git show <ref>:<path>` -- nothing is
hardcoded from memory. Run lung_tcell_talk_imrad_figures_20260925.py first
to stage the JSDP image assets; the four reused 2026-09-24 figures are
referenced directly from their original directory, not copied.

Full citations: lung_tcell_talk_imrad_sources_20260925.md.
Byline matches poster/poster_final and talk/JSDP_P25_talk exactly -- no
addition pending the open human ruling on AI co-authorship.

STANDING RULE, sharpened 2026-09-24/25: the confound/caveat is named in the
same sentence as the result it qualifies, at FIRST utterance -- not a later
slide, not notes, not Q&A. Two corrections to how that rule was almost
mis-applied this round, both checked at source before building:
  (1) S100A8/A9 have zero presence in JSDP_P25_talk or poster_final -- they
      are not an "original claim now qualified", they simply enter the
      story for the first time already caveated (Results, R4).
  (2) The T2/T5 SCLC-vs-LUAD ordering was already flagged as an OPEN,
      PENDING test in both original documents (poster: "remains an open
      question rather than a reconciled conclusion"; talk, under a shape
      reading "TESTED SINCE THE POSTER": "Pending T2 and T3") -- R3
      completes a test the original team flagged as open, framed as
      resolution, not reversal.

Timestamps in this deck/sources note: JST (UTC+9) first, source UTC form in
brackets, per the human's standing ruling (2026-09-25). Calendar-only dates
are unaffected.

R8 (mechanism/STRING overlay) is the DESIGNATED CUT for a shorter slot --
see INCLUDE_R8_MECHANISM below.
"""
from __future__ import annotations

import csv
import io
import json
import subprocess
from datetime import datetime, timezone, timedelta

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas as pdfcanvas
from reportlab.platypus import Paragraph, Table, TableStyle, Frame, Image as RLImage

REPO = "/Users/kaisar/workspace/github/geneformer-lung-tcell"
ASSETS = "/Users/kaisar/workspace/office_hive/hive/reports/lung_tcell_talk_imrad_assets_20260925"
ASSETS_24 = "/Users/kaisar/workspace/office_hive/hive/reports/lung_tcell_talk_assets_20260924"
OUT = "/Users/kaisar/workspace/office_hive/hive/reports/lung_tcell_talk_imrad_20260925.pdf"

# Designated cut for a shorter slot: flip to False to drop R8 (mechanism)
# without touching anything else -- a one-line change, per Michael's
# instruction that a 15-minute version must not require a rewrite.
INCLUDE_R8_MECHANISM = True

JST = timezone(timedelta(hours=9))


def jst_utc(dt_str: str) -> str:
    """Format a git-style timestamp (e.g. '2026-09-24 18:57:39 +0900') as
    'HH:MM JST (HH:MMZ)'. Git commit times from this repo are already
    +0900 (JST); this makes that explicit and derives the UTC form rather
    than asserting it."""
    dt = datetime.strptime(dt_str[:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=JST)
    utc = dt.astimezone(timezone.utc)
    return f"{dt.strftime('%Y-%m-%d %H:%M')} JST ({utc.strftime('%H:%MZ')})"


def git_show(ref: str, path: str) -> str:
    out = subprocess.run(["git", "show", f"{ref}:{path}"], cwd=REPO,
                          capture_output=True, check=True, text=True)
    return out.stdout


def git_log_time(ref: str) -> str:
    out = subprocess.run(["git", "log", "-1", "--format=%ci", ref], cwd=REPO,
                          capture_output=True, check=True, text=True)
    return out.stdout.strip()


def read_csv(ref: str, path: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(git_show(ref, path))))


# --------------------------------------------------------------------------
# Pull every number used on a slide, at build time, from source.
# --------------------------------------------------------------------------

top_candidates = read_csv("b99365b", "sclc_validation/primary_test_perturbation/tables/isp_plausibility_top_candidates.csv")
s100_rows = [r for r in top_candidates if r["Gene_name"] in ("S100A8", "S100A9")]
assert len(s100_rows) == 12
N_SAME_SIGN = sum(1 for r in s100_rows if float(r["delete_shift"]) * float(r["overexpress_shift"]) > 0)
N_CONCORDANT = sum(1 for r in s100_rows if r["concordant"] == "True")
S100A8_luad_sclc = next(r for r in s100_rows if r["Gene_name"] == "S100A8" and r["comparison"] == "luad_to_sclc")
FDR_EXAMPLE = float(S100A8_luad_sclc["delete_fdr"])
N_EXAMPLE = int(float(S100A8_luad_sclc["delete_n"]))
SHIFT_EXAMPLE = float(S100A8_luad_sclc["delete_shift"])

top120 = [r for r in top_candidates if r["selected_by"] == "top_n_by_effect"]
assert len(top120) == 120
N_AMBIENT_FLAG = sum(1 for r in top120 if r["ambient_flag"] == "True")
N_AMBIENT_ANCHOR = sum(1 for r in top120 if r["is_ambient_anchor"] == "True")
N_BOTH = sum(1 for r in top120 if r["is_ambient_anchor"] == "True" and r["ambient_flag"] == "True")
N_FLAG_ONLY = N_AMBIENT_FLAG - N_BOTH
N_ANCHOR_ONLY = N_AMBIENT_ANCHOR - N_BOTH
N_UNION = N_FLAG_ONLY + N_ANCHOR_ONLY + N_BOTH
manifest = json.loads(git_show("b99365b", "sclc_validation/primary_test_perturbation/tables/ambient_risk_manifest.json"))
N_TCELL_ANCHORS = len(manifest["anchors_tcell"])
N_TCELL_PRESENT = sum(1 for r in top120 if r["is_tcell_anchor"] == "True")

nbeal1_rows = [r for r in top_candidates if r["Gene_name"] == "NBEAL1"]
N_NBEAL1 = len(nbeal1_rows)
N_NBEAL1_CONCORDANT = sum(1 for r in nbeal1_rows if r["concordant"] == "True")
NBEAL1_DETECT = float(nbeal1_rows[0]["detect_frac"])

t6_recon = read_csv("9ec518c", "sclc_validation/immune_axis_test/results/t6_weighting_reconciliation.csv")
t6_donors = read_csv("9ec518c", "sclc_validation/immune_axis_test/results/t6_donor_level_scores.csv")
t6_perm = read_csv("9ec518c", "sclc_validation/immune_axis_test/results/t6_permutation_tests.csv")


def recon(pop, state):
    return next(r for r in t6_recon if r["population"] == pop and r["state"] == state)


COMPLETE_SCLC = recon("complete", "sclc")
COMPLETE_LUAD = recon("complete", "luad")
TEST_SCLC = recon("test_only", "sclc")
TEST_LUAD = recon("test_only", "luad")
PE_SHARE = float(TEST_SCLC["max_donor_cell_share"])

COMPLETE_TEST = next(r for r in t6_perm if r["comparison"] == "complete_sclc_vs_luad")
TESTONLY_TEST = next(r for r in t6_perm if r["comparison"] == "test_only_sclc_vs_luad")
P_COMPLETE = float(COMPLETE_TEST["p_two_sided"])
N_DONORS_A = int(COMPLETE_TEST["n_donors_a"])
N_DONORS_B = int(COMPLETE_TEST["n_donors_b"])
FLOOR_P = float(TESTONLY_TEST["min_attainable_two_sided_p"])
N_ASSIGN = int(TESTONLY_TEST["n_assignments"])

METHODS_MD = git_show("origin/main", "sclc_validation/perturbation_workflow/METHODS.md")
assert "RU675" in METHODS_MD and "42 donors" in METHODS_MD
assert "Normal-class\nmetrics are flagged explicitly given 1 test donor." in METHODS_MD or \
    "metrics are flagged explicitly given 1 test donor." in METHODS_MD

BF16_RESULTS = git_show("afa0564", "sclc_validation/bf16_bench/RESULTS_BF16.md")
assert "0.9998" in BF16_RESULTS and "FAIL" in BF16_RESULTS

# Base-branch provenance (freshly fetched at build time, not assumed).
ORIGIN_MAIN_SHA = subprocess.run(["git", "rev-parse", "origin/main"], cwd=REPO,
                                  capture_output=True, check=True, text=True).stdout.strip()
ORIGIN_MAIN_TIME = jst_utc(git_log_time("origin/main"))

BUILD_TIME = jst_utc(datetime.now(JST).strftime("%Y-%m-%d %H:%M:%S %z").replace("+0900", "+0900"))

# --------------------------------------------------------------------------
# House-style helpers (from make_lung_tcell_talk_slides_20260924.py /
# hive/reports/make_s100_slides_pdf.py) -- unchanged, reused verbatim.
# --------------------------------------------------------------------------

W, H = 338 * mm, 190 * mm
INK = colors.HexColor("#12211E")
PAPER = colors.HexColor("#F6F8F6")
WHITE = colors.HexColor("#FFFFFF")
LINE = colors.HexColor("#DCE3DF")
BAND = colors.HexColor("#EAF0ED")
ACC = colors.HexColor("#0E6E5C")
WARM = colors.HexColor("#B4552F")
DARK = colors.HexColor("#12211E")
DARKP = colors.HexColor("#2A201B")
WARMT = colors.HexColor("#E08A5F")
CREAM = colors.HexColor("#F7F2EE")
MUT = colors.HexColor("#5C6B66")
M = 20 * mm
PAGES = 18


def st(size, color=INK, bold=False, leading=None, space=3):
    return ParagraphStyle("s%d%s%s" % (size, bold, id(color)), fontName="Helvetica-Bold" if bold else "Helvetica",
                           fontSize=size, leading=leading or size * 1.32, textColor=color,
                           alignment=TA_LEFT, spaceAfter=space)


def frame(c, x, y, w, h, flows):
    Frame(x, y, w, h, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0,
          showBoundary=0).addFromList(list(flows), c)


def card(c, x, y, w, h, fill=WHITE, stroke=LINE, accent=None):
    c.setFillColor(fill)
    c.setStrokeColor(stroke)
    c.setLineWidth(0.6)
    c.roundRect(x, y, w, h, 4 * mm, stroke=1, fill=1)
    if accent:
        c.setFillColor(accent)
        c.roundRect(x, y, 2.4 * mm, h, 1.2 * mm, stroke=0, fill=1)


def footer(c, n, right=None, dark=False):
    col = colors.HexColor("#9A877C") if dark else colors.HexColor("#7A8A85")
    c.setFillColor(col)
    c.setFont("Courier", 8.5)
    c.drawString(M, 11 * mm, "%d / %d" % (n, PAGES))
    if right:
        c.setFont("Helvetica", 8.5)
        c.drawRightString(W - M, 11 * mm, right)


def bg(c, col=PAPER):
    c.setFillColor(col)
    c.rect(0, 0, W, H, stroke=0, fill=1)


def eyebrow(c, txt, y, col=ACC):
    c.setFillColor(col)
    c.setFont("Courier-Bold", 9)
    c.drawString(M, y, txt.upper())


def title(c, txt, y, size=25, col=INK, width=None):
    h = size * 2.6
    frame(c, M, y - h, width or (W - 2 * M), h, [Paragraph(txt, st(size, col, True, size * 1.16))])


def img(c, path, x, y, w, h):
    frame(c, x, y, w, h, [RLImage(path, width=w, height=h, kind="proportional")])


def three_col(c, n0, items, top=True):
    cw = (W - 2 * M - 12 * mm) / 3
    for i, (hd, txt) in enumerate(items):
        x = M + i * (cw + 6 * mm)
        y0 = H - 108 * mm if top else 96 * mm
        card(c, x, y0, cw, 56 * mm)
        frame(c, x + 7 * mm, y0 + 4 * mm, cw - 14 * mm, 48 * mm, [
            Paragraph(hd, st(12.5, INK, True, 16)),
            Paragraph(txt, st(9.6, colors.HexColor("#4A5A55"), False, 13.4)),
        ])


c = pdfcanvas.Canvas(OUT, pagesize=(W, H))
c.setTitle("A Foundation-Model T-cell Dysfunction Screen: Translational Candidates, Then a Full Audit")
c.setAuthor("Kaisar Dauyey; Shinji Nakaoka")
c.setSubject("IMRaD restructure -- geneformer-lung-tcell, 2026-09-25")

# ============================================================ 1. TITLE ===
bg(c, DARK)
eyebrow(c, "Geneformer · SCLC/LUAD/normal T-cell dysfunction", H - 22 * mm, WARMT)
frame(c, M, H - 100 * mm, W - 2 * M, 70 * mm, [
    Paragraph("A Foundation-Model T-cell Dysfunction Screen:", st(27, CREAM, True, 32)),
    Paragraph("Translational Candidates, Then a Full Audit", st(27, CREAM, True, 32)),
    Paragraph("Four checkpoint/persistence hits, tissue-validated -- and what a full concordance, "
              "ambient, and donor-level audit does and does not change about them.",
              st(12.5, WARMT, False, 17)),
])
frame(c, M, 24 * mm, W - 2 * M, 26 * mm, [
    Paragraph("<b>Kaisar Dauyey</b> &nbsp;&#183;&nbsp; <b>Shinji Nakaoka</b>", st(13, CREAM, False, 17)),
    Paragraph("Laboratory of Mathematical Biology, Faculty of Advanced Life Science, Hokkaido University "
              "&nbsp;&#183;&nbsp; snakaoka@sci.hokudai.ac.jp", st(10.5, WARMT, False, 15)),
    Paragraph(f"geneformer-lung-tcell &nbsp;&#183;&nbsp; restructured from poster_final (2026-08-15) and "
              f"JSDP_P25_talk (2026-08-17) &nbsp;&#183;&nbsp; built {BUILD_TIME}",
              st(9, colors.HexColor("#B8A79A"), False, 13)),
])
img(c, ASSETS + "/qr.png", W - M - 26 * mm, H - 40 * mm, 20 * mm, 20 * mm)
footer(c, 1, "DRAFT -- outline stage, not camera-ready", dark=True)
c.showPage()

# ======================================================= 2. BACKGROUND ===
bg(c)
eyebrow(c, "Background", H - 24 * mm)
title(c, "A perturbation screen can look statistically overwhelming and still be wrong", H - 42 * mm, 21, width=270*mm)
frame(c, M, 30 * mm, W - 2 * M, 92 * mm, [
    Paragraph("Foundation models fine-tuned on single-cell transcriptomes can be probed for candidate "
              "disease-state drivers: delete or overexpress a gene's simulated expression in a "
              "held-out cell, and read the shift in a fine-tuned classifier's predicted state. A "
              "dose-responsive driver should move the model toward its target when overexpressed and "
              "away from it when deleted.", st(12.5, colors.HexColor("#31423D"), False, 18)),
    Paragraph("A classifier trained on real tissue also learns whatever else correlates with its "
              "labels -- ambient RNA from other cell types, detection-rate differences, cohort "
              "composition, tissue of origin. None of that is cell-intrinsic biology, and a screen "
              "does not tell you on its own which one it found.", st(12.5, colors.HexColor("#31423D"), False, 18)),
        Paragraph("This talk presents a set of translational candidates identified this way, then applies "
              "the audit a screen like this needs: a concordance test, a screen-level ambient-contamination "
              "check, and a donor-level reanalysis of a cross-cohort ordering claim.",
              st(12.5, INK, True, 18)),
])
footer(c, 2)
c.showPage()

# =================================================== 3. INTRODUCTION I1 ===
bg(c)
eyebrow(c, "Introduction · the original study", H - 24 * mm)
title(c, "Four replicated checkpoint/persistence hits, tissue-validated", H - 42 * mm, 21, width=270*mm)
img(c, ASSETS + "/screen_b.png", M, 30 * mm, 130 * mm, 100 * mm)
frame(c, M + 138 * mm, 70 * mm, W - M - 138 * mm - M, 60 * mm, [
    Paragraph("A donor-disjoint SCLC/LUAD/normal T-cell classifier supported an applied in-silico "
              "perturbation screen for candidate CAR-T/ICI-relevant drivers. Four hits replicated "
              "across the screen: <b>TIGIT, TIM-3 (HAVCR2), CTLA-4, IL7R.</b>",
              st(11.5, colors.HexColor("#31423D"), False, 16)),
    Paragraph("Independent Visium spatial validation (GSE263196, 5 SCLC specimens, 15,632 in-tissue "
              "spots) supported a tissue-validated antigen-presentation programme.",
              st(11.5, colors.HexColor("#31423D"), False, 16)),
])
card(c, M + 138 * mm, 30 * mm, W - M - 138 * mm - M, 36 * mm, BAND, BAND)
frame(c, M + 146 * mm, 34 * mm, W - M - 138 * mm - M - 16 * mm, 28 * mm, [
    Paragraph("Presented at JSDP, 2026-08-17. No caveat attached here -- none of these four hits is "
              "implicated in anything that follows.", st(10, colors.HexColor("#31423D"), True, 14)),
])
footer(c, 3, "talk/JSDP_P25_talk.pptx; poster/poster_final.html")
c.showPage()

# =================================================== 4. INTRODUCTION I2 ===
bg(c)
eyebrow(c, "Introduction · what was already open, and what this talk adds", H - 24 * mm)
title(c, "One test was already flagged as pending; one new screen needed the same scrutiny", H - 44 * mm, 17, width=290*mm)
three_col(c, 0, [
    ("The poster's own hedge",
     "“The model implies Normal &lt; SCLC &lt; LUAD on the checkpoint axis. This conflicts with "
     "the clinical tumour-level picture of SCLC as cold and ICI-resistant, so it <b>remains an open "
     "question rather than a reconciled conclusion.</b>”"),
    ("The talk's own QC slide, “TESTED SINCE THE POSTER”",
     "“The states are not on a line, so there is no ordering to contradict... Test 1 of 4 &#8212; "
     "sign directions, not centroid geometry. <b>Pending T2 and T3.</b>” The donor-level test was "
     "named as open before this work began."),
    ("What this talk adds",
     "A later, expanded screen surfaced two new high-FDR candidates needing the same concordance "
     "scrutiny the original four received. This talk <b>completes the pending test</b> and applies "
     "that scrutiny &#8212; it does not reopen a settled claim."),
])
footer(c, 4, "poster/poster_final.html; talk/JSDP_P25_talk.pdf (\"TESTED SINCE THE POSTER\", p.7)")
c.showPage()

# ============================================== 5. METHODS M1 · cohort ===
bg(c)
eyebrow(c, "Methods · cohort, classifier", H - 24 * mm)
title(c, "46,140 T cells, 42 donors, one fine-tuned classifier", H - 42 * mm, 22, width=200*mm)
img(c, ASSETS + "/confusion.png", M, 26 * mm, 120 * mm, 96 * mm)
frame(c, M + 128 * mm, 84 * mm, W - M - 128 * mm - M, 40 * mm, [
    Paragraph("Held-out test performance", st(12.5, ACC, True, 16)),
    Paragraph("<font face='Courier-Bold' size=20>91.9%</font> accuracy, macro F1 0.903, "
              "donor-disjoint split. Geneformer V2 104M fine-tune.", st(10.5, INK, False, 15)),
])
card(c, M + 128 * mm, 26 * mm, W - M - 128 * mm - M, 52 * mm, BAND, BAND)
frame(c, M + 136 * mm, 30 * mm, W - M - 128 * mm - M - 16 * mm, 44 * mm, [
    Paragraph("Normal-class test performance rests on 1 held-out donor (566 cells)", st(11, ACC, True, 15)),
    Paragraph("A third of this macro-F1. Already flagged in the repo's own methods "
              "(perturbation_workflow/METHODS.md, lines 57-67 the split table, line 112: “Normal-class "
              "metrics are flagged explicitly given 1 test donor”).",
              st(9.4, colors.HexColor("#31423D"), False, 13)),
    Paragraph("19 SCLC + 22 LUAD + 4 normal = 45, but 42 unique donors: three donor_ids (RU675, RU682, "
              "RU684) each carry two biospecimens from the same clinical patient, one LUAD tumour and "
              "one normal-tissue sample (METHODS.md, lines 33-39).",
              st(9.4, colors.HexColor("#31423D"), False, 13)),
])
footer(c, 5, "README.md; perturbation_workflow/METHODS.md")
c.showPage()

# ======================================== 6. METHODS M2 · perturbation ===
bg(c)
eyebrow(c, "Methods · in-silico perturbation", H - 24 * mm)
title(c, "Delete or overexpress, in held-out cells, across six comparisons", H - 44 * mm, 21, width=270*mm)
frame(c, M, 30 * mm, W - 2 * M, 92 * mm, [
    Paragraph("Each gene is deleted (its token removed from a cell's rank-ordered input) or overexpressed "
              "(promoted to the front of the rank) in held-out cells, and the shift in the classifier's "
              "predicted disease-state probability is read for both operations.",
              st(12.5, colors.HexColor("#31423D"), False, 18)),
    Paragraph("Six source&#8594;goal disease-state comparisons: SCLC&#8594;LUAD, SCLC&#8594;normal, "
              "LUAD&#8594;SCLC, LUAD&#8594;normal, normal&#8594;SCLC, normal&#8594;LUAD.",
              st(12.5, colors.HexColor("#31423D"), False, 18)),
    Paragraph("A dose-responsive driver should move the model toward its target under overexpression and "
              "away from it under deletion &#8212; the two operations are read together, never one alone.",
              st(12.5, INK, True, 18)),
])
footer(c, 6, "README.md")
c.showPage()

# ========================================= 7. METHODS M3 · concordance ===
bg(c)
eyebrow(c, "Methods · how a hit earns the name", H - 24 * mm)
title(c, "A real driver must move opposite ways under deletion and overexpression", H - 44 * mm, 19, width=280*mm)
card(c, M, 90 * mm, W - 2 * M, 46 * mm, WHITE, LINE, ACC)
frame(c, M + 12 * mm, 94 * mm, W - 2 * M - 24 * mm, 38 * mm, [
    Paragraph("The <font face='Courier'>concordant</font> criterion, used throughout this repo's "
              "perturbation work:", st(12, INK, True, 16)),
    Paragraph("<font face='Courier'>delete_fdr &lt; 0.05 AND overexpress_fdr &lt; 0.05 AND "
              "delete_shift &#215; overexpress_shift &lt; 0 AND delete_n &#8805; 25</font>",
              st(11.5, WARM, False, 16)),
])
frame(c, M, 30 * mm, W - 2 * M, 50 * mm, [
    Paragraph("Significance in one arm is not enough. A gene must clear both arms and reverse sign "
              "between them &#8212; the signature a dose-responsive driver actually produces. A gene can "
              "be significant in both arms and still fail this test if deletion and overexpression push "
              "the classifier the <b>same</b> direction &#8212; the model reacting to the gene's presence "
              "at all, not to a dose-responsive relationship with the labelled state.",
              st(12, colors.HexColor("#31423D"), False, 17)),
])
footer(c, 7, "primary_test_perturbation/scripts/build_delete_overexpress_shift_report.py @ b99365b")
c.showPage()

# ============================================= 8. METHODS M4 · audits ===
bg(c)
eyebrow(c, "Methods · screen-level and donor-level audits", H - 24 * mm)
title(c, "Two further checks, applied after the screen, not instead of it", H - 44 * mm, 20, width=280*mm)
three_col(c, 0, [
    ("Ambient-risk audit",
     "Every top-effect gene is checked against a curated ambient-contamination panel and a curated "
     "T-cell-intrinsic anchor panel. The flag threshold is the 25th percentile of the full ambient-anchor "
     "population's own risk scores &#8212; a within-panel, not an absolute, cutoff."),
    ("Donor-weighted reanalysis",
     "Any cross-cohort comparison is recomputed donor-level (mean within donor first, then an unweighted "
     "mean of donor means) alongside the original cell-weighted statistic, to separate a real group "
     "difference from one donor's cell count dominating the mean."),
    ("Why both, not either",
     "A significance test alone cannot distinguish a cell-intrinsic effect from an ambient or "
     "cell-count artifact. Both audits are applied to every candidate discussed in Results, not "
     "selectively to ones that already look questionable."),
])
footer(c, 8, "ambient_risk_diagnostic.py line 198 @ b99365b; immune_axis_test/donor_robustness.py @ 9ec518c")
c.showPage()

# ========================================== 9. RESULTS R1 · original 4 ===
bg(c)
eyebrow(c, "Results · the original screen", H - 24 * mm)
title(c, "Four hits replicate and validate independently in tissue", H - 42 * mm, 21, width=260*mm)
img(c, ASSETS + "/spatial.png", M, 26 * mm, 130 * mm, 100 * mm)
frame(c, M + 138 * mm, 66 * mm, W - M - 138 * mm - M, 60 * mm, [
    Paragraph("TIGIT, TIM-3 (HAVCR2), CTLA-4, IL7R", st(15, ACC, True, 19)),
    Paragraph("Four replicated checkpoint/persistence edits from the applied screen (screen_b.png, "
              "Introduction).", st(10.5, colors.HexColor("#31423D"), False, 15)),
    Paragraph("GSE263196 Visium, 5 SCLC specimens, 15,632 in-tissue spots: antigen-presentation "
              "programme <b>rho = 0.361, 7.95 sigma above the null.</b>",
              st(10.5, colors.HexColor("#31423D"), False, 15)),
])
card(c, M + 138 * mm, 26 * mm, W - M - 138 * mm - M, 32 * mm, BAND, BAND)
frame(c, M + 146 * mm, 29 * mm, W - M - 138 * mm - M - 16 * mm, 26 * mm, [
    Paragraph("No caveat on this slide. These four hits and this spatial-validation number are "
              "unaffected by anything in the rest of Results.", st(9.6, colors.HexColor("#31423D"), True, 13.5)),
])
footer(c, 9, "talk/JSDP_P25_talk.pptx/pdf (rho, sigma: talk only, 0 hits in poster_final)")
c.showPage()

# =========================================== 10. RESULTS R2 · own QC ===
bg(c)
eyebrow(c, "Results · the original talk's own quality-control check", H - 24 * mm)
title(c, "Sparse genes fake big effects, and the poster's own ordering fails its own test", H - 42 * mm, 17, width=290*mm)
img(c, ASSETS + "/detection.png", M, 26 * mm, 150 * mm, 100 * mm)
frame(c, M + 158 * mm, 60 * mm, W - M - 158 * mm - M, 66 * mm, [
    Paragraph("Detection vs. effect size: Spearman rho = -0.60 (p = 2.7e-5, n = 42) across the full "
              "50-gene panel.", st(11, ACC, True, 15)),
    Paragraph("“TESTED SINCE THE POSTER” (talk, slide 7): “the states are not on a line, so "
              "there is no ordering to contradict... Working reading: a triangle &#8212; SCLC &#8804; Normal "
              "&lt; LUAD... Test 1 of 4 &#8212; sign directions, not centroid geometry. <b>Pending T2 and "
              "T3.</b>”", st(10, colors.HexColor("#31423D"), False, 14.5)),
])
card(c, M + 158 * mm, 26 * mm, W - M - 158 * mm - M, 30 * mm, BAND, BAND)
frame(c, M + 166 * mm, 29 * mm, W - M - 158 * mm - M - 16 * mm, 24 * mm, [
    Paragraph("This is the original team's own self-correction, dated before this work began &#8212; "
              "not a finding of ours.", st(9.4, colors.HexColor("#31423D"), True, 13)),
])
footer(c, 10, "talk/JSDP_P25_talk.pdf, slide 7")
c.showPage()

# ================================== 11. RESULTS R3 · completes the test ===
bg(c)
eyebrow(c, "Results · completing the pending test", H - 24 * mm)
title(c, "The open question is now closed as null", H - 44 * mm, 24)
img(c, ASSETS_24 + "/fig_t6_reversal.png", M, 40 * mm, 195 * mm, 80 * mm)
frame(c, M + 203 * mm, H - 60 * mm, W - M - 203 * mm - M, 48 * mm, [
    Paragraph("Donor-level: no reversal", st(13, ACC, True, 17)),
    Paragraph(f"Complete population: SCLC {float(COMPLETE_SCLC['donor_level_mean']):.3f} &gt; LUAD "
              f"{float(COMPLETE_LUAD['donor_level_mean']):.3f}. Test-only: SCLC "
              f"{float(TEST_SCLC['donor_level_mean']):.3f} &gt; LUAD {float(TEST_LUAD['donor_level_mean']):.3f}. "
              f"Two-sided <b>Monte Carlo</b> test, 100,000 replicates: p = {P_COMPLETE:.3f}, "
              f"{N_DONORS_A} vs. {N_DONORS_B} donors. <b>Not significant.</b>",
              st(9.8, colors.HexColor("#4A5A55"), False, 13.6)),
])
card(c, M + 203 * mm, 40 * mm, W - M - 203 * mm - M, 78 * mm, BAND, BAND)
frame(c, M + 211 * mm, 44 * mm, W - M - 203 * mm - M - 16 * mm, 70 * mm, [
    Paragraph("PleuralEffusion: one genuine donor", st(11.5, ACC, True, 15)),
    Paragraph(f"Carries {PE_SHARE*100:.1f}% of test-only SCLC cells. <font face='Courier'>"
              f"HTAN_Participant_ID</font> confirms exactly one participant (<font face='Courier'>"
              f"HTA8_2001</font>), not several merged (METHODS.md, lines 29-32).",
              st(9.4, colors.HexColor("#31423D"), False, 13.2)),
    Paragraph(f"A separate, <b>exact-enumeration</b> test on the 3-vs-4-donor test-only slice: only "
              f"{N_ASSIGN} label assignments exist; the smallest attainable two-sided p is "
              f"{FLOOR_P:.4f} (2/{N_ASSIGN}) &#8212; that design cannot reach significance regardless "
              f"of the true effect.", st(9.2, colors.HexColor("#31423D"), False, 13)),
])
footer(c, 11, "t6_weighting_reconciliation.csv, t6_permutation_tests.csv @ 9ec518c")
c.showPage()

# ================================================ 12. RESULTS R4 · S100 ===
bg(c)
eyebrow(c, "Results · a later, expanded screen", H - 22 * mm)
title(c, "S100A8 and S100A9: an overwhelming FDR, and no concordance at all", H - 32 * mm, 20, width=152 * mm)
img(c, ASSETS_24 + "/fig_concordance.png", M, 40 * mm, 150 * mm, 80 * mm)
frame(c, M + 158 * mm, 92 * mm, W - M - 158 * mm - M, 28 * mm, [
    Paragraph(f"{N_SAME_SIGN} of 12 rows same-sign; {N_CONCORDANT} of 6 comparisons concordant",
              st(13, WARM, True, 16)),
    Paragraph(f"Example: S100A8, LUAD&#8594;SCLC deletion &#8212; FDR = {FDR_EXAMPLE:.2e}, "
              f"N = {N_EXAMPLE} detections, shift = {SHIFT_EXAMPLE:+.4f}. "
              "<b>The FDR is a certainty claim, not a size claim</b> &#8212; it says nothing about "
              "direction. First appearance of these two genes in this deck: caveated from the start, "
              "not a later correction.", st(9.0, colors.HexColor("#4A5A55"), False, 12.6)),
])
card(c, M + 158 * mm, 40 * mm, W - M - 158 * mm - M, 48 * mm, BAND, BAND)
frame(c, M + 166 * mm, 43 * mm, W - M - 158 * mm - M - 16 * mm, 42 * mm, [
    Paragraph("What the LUAD-normal and SCLC-normal rows can and cannot lean on", st(11, ACC, True, 14)),
    Paragraph("Of 12 rows, 8 touch the normal class (4 comparisons). For luad&#8596;normal, 3 of 4 "
              "normal donors are the <i>same patients</i> as 3 LUAD donors &#8212; real protection "
              "against a separate-cohort confound, though only 3 of 4. <b>No such pairing exists for "
              "sclc&#8596;normal.</b>", st(9.2, colors.HexColor("#31423D"), False, 12.6)),
])
footer(c, 12, "isp_plausibility_top_candidates.csv @ b99365b")
c.showPage()

# ============================================== 13. RESULTS R5 · ambient ===
bg(c, DARK)
eyebrow(c, "Results · screen-level ambient audit", H - 22 * mm, WARMT)
title(c, "Half the top 120 is ambient-flagged or a curated ambient anchor", H - 40 * mm, 20, CREAM, width=152 * mm)
img(c, ASSETS_24 + "/fig_ambient_breakdown.png", M, 40 * mm, 150 * mm, 78 * mm)
for i, (big, small) in enumerate([
    (f"{N_UNION}/120", f"rows are ambient-flagged, a curated ambient anchor, or both: {N_FLAG_ONLY} "
     f"flagged-only, {N_ANCHOR_ONLY} anchor-only, {N_BOTH} both &#8212; a real partition of the 120, "
     "not two overlapping counts added together"),
    (f"{N_ANCHOR_ONLY}/120", "rows are anchors below the flag threshold &#8212; mechanical, not a "
     "miss: the threshold is the anchors' own 25th percentile of risk, so some anchors falling under "
     "it is built in."),
    (f"{N_TCELL_PRESENT}/{N_TCELL_ANCHORS}", "of the 36 curated T-cell-intrinsic anchor genes (a "
     "different denominator: the full curated list, not a row count) appear anywhere in the 120."),
]):
    x = M + 158 * mm
    y = H - 46 * mm - i * 36 * mm
    card(c, x, y - 30 * mm, W - x - M, 30 * mm, DARKP, DARKP)
    frame(c, x + 6 * mm, y - 26 * mm, W - x - M - 12 * mm, 22 * mm, [
        Paragraph(f"<font face='Courier-Bold' size=19 color='#E08A5F'>{big}</font>", st(19, WARMT, False, 23)),
        Paragraph(small, st(8.1, colors.HexColor("#D6C7BE"), False, 10.9)),
    ])
footer(c, 13, "top_n_by_effect rows, ambient_risk_manifest.json @ b99365b", dark=True)
c.showPage()

# =============================================== 14. RESULTS R6 · NBEAL1 ===
bg(c)
eyebrow(c, "Results · not a null study", H - 24 * mm)
title(c, "One candidate survives every check", H - 38 * mm, 24)
card(c, M, H - 118 * mm, W - 2 * M, 56 * mm, WHITE, LINE, ACC)
frame(c, M + 12 * mm, H - 110 * mm, W - 2 * M - 24 * mm, 44 * mm, [
    Paragraph("NBEAL1", st(26, ACC, True, 30)),
    Paragraph(f"{N_NBEAL1_CONCORDANT} of {N_NBEAL1} comparisons concordant &#8212; every one it appears "
              f"in &#8212; detected in {NBEAL1_DETECT*100:.1f}% of cells.", st(13, INK, False, 18)),
    Paragraph("Passes the same concordance test that eliminates S100A8/S100A9, at a detection rate well "
              "above the ambient range flagged for those two. Does not make it a validated driver on "
              "its own &#8212; makes it the one candidate this screen has not been able to eliminate.",
              st(10.5, colors.HexColor("#4A5A55"), False, 14.5)),
])
frame(c, M, 22 * mm, W - 2 * M, 22 * mm, [
    Paragraph("Not everything restored from the original screen fails, and the one thing that doesn't "
              "fail is more informative for having survived a real test.",
              st(11.5, colors.HexColor("#4A5A55"), False, 15.5)),
])
footer(c, 14, "isp_plausibility_top_candidates.csv @ b99365b")
c.showPage()

# ================================================= 15. RESULTS R7 · bf16 ===
bg(c, DARK)
eyebrow(c, "Results · control, is any of this a precision artifact?", H - 24 * mm, WARMT)
title(c, "No: precision replicates cleanly; model size is the axis that moves things", H - 40 * mm, 17, CREAM)
img(c, ASSETS_24 + "/fig_bf16_canary.png", M, 40 * mm, 110 * mm, 78 * mm)
frame(c, M + 118 * mm, 76 * mm, W - M - 118 * mm - M, 42 * mm, [
    Paragraph("Same model, fp32 vs. bf16 (half precision)", st(12.5, WARMT, True, 16)),
    Paragraph("104M panel (300 rows): rho = 0.9998, sign agreement 298/300. 316M canary (60 rows): "
              "rho = 0.9987, sign agreement 60/60.", st(9.6, colors.HexColor("#CFE0DA"), False, 13.4)),
    Paragraph("The pre-registered 104M sign-agreement gate is still recorded <b>FAIL</b> &#8212; a "
              "bitwise-deterministic fp32-vs-fp32 rerun made the noise floor zero, not a usable sign "
              "threshold; the domain reviewer's ruling was that the frozen gate stands unamended, not "
              "that bf16 disagrees with fp32. Both are true at once.",
              st(9.2, colors.HexColor("#CFE0DA"), False, 12.8)),
])
frame(c, M + 118 * mm, 40 * mm, W - M - 118 * mm - M, 32 * mm, [
    Paragraph("Same model, changed precision vs. changed size", st(12.5, WARMT, True, 16)),
    Paragraph("Precision (104M fp32&#8594;bf16): rho = 0.9998. Size (104M&#8594;316M, both bf16): "
              "rho = 0.489. Using (1&#8722;rho) as the divergence measure, size moves results ~2,500x "
              "(~3.4 orders of magnitude) more than precision. None of this deck's findings are "
              "precision artifacts.", st(9.4, colors.HexColor("#CFE0DA"), False, 13.2)),
])
footer(c, 15, "committed_run_stats/{bf16,316m_bf16,...} @ origin/main; RESULTS_BF16.md @ afa0564", dark=True)
c.showPage()

# ============================================= 16. RESULTS R8 · mechanism ===
if INCLUDE_R8_MECHANISM:
    bg(c)
    eyebrow(c, "Results · mechanism (designated cut for a shorter slot)", H - 24 * mm)
    title(c, "The four original hits sit in one connected interaction neighbourhood", H - 42 * mm, 18, width=270*mm)
    img(c, ASSETS + "/network_c.png", M, 26 * mm, 150 * mm, 100 * mm)
    frame(c, M + 158 * mm, 60 * mm, W - M - 158 * mm - M, 66 * mm, [
        Paragraph("From prior literature, not inferred from these cells", st(12, ACC, True, 16)),
        Paragraph("TIGIT: CAR-T support. TIM-3: solid-tumour knockout precedent. CTLA-4: real signal, "
                  "less CAR-T-specific. IL7R: persistence, not exhaustion. 54 STRING edges among 16 "
                  "context genes; TOX and LAYN have no edge above threshold; PD-L1 is absent from the "
                  "T-cell atlas.", st(10, colors.HexColor("#31423D"), False, 14.5)),
    ])
    card(c, M + 158 * mm, 26 * mm, W - M - 158 * mm - M, 30 * mm, BAND, BAND)
    frame(c, M + 166 * mm, 29 * mm, W - M - 158 * mm - M - 16 * mm, 24 * mm, [
        Paragraph("No caveat &#8212; these four hits are unaffected by anything in this deck's audit.",
                  st(9.4, colors.HexColor("#31423D"), True, 13)),
    ])
    footer(c, 16, "talk/JSDP_P25_talk.pdf, slide 10 (STRING literature overlay)")
    c.showPage()

# =================================================== 17. DISCUSSION D1 ===
bg(c)
eyebrow(c, "Discussion · synthesis", H - 24 * mm)
title(c, "What this establishes, and what it rules out", H - 42 * mm, 22)
three_col(c, 0, [
    ("Establishes",
     "Four checkpoint/persistence hits (TIGIT, TIM-3, CTLA-4, IL7R) replicate and validate "
     "independently in spatial tissue data. A concordance test and a screen-level ambient audit are "
     "necessary companions to significance for any further perturbation-based driver claim. One "
     "further candidate (NBEAL1) clears every check applied."),
    ("Resolves",
     "The checkpoint-axis ordering the original team flagged as open (“pending T2 and T3”) is "
     "closed: donor-level SCLC-vs-LUAD shows no significant difference in either direction, once a "
     "single donor's cell-count share is accounted for."),
    ("Rules out",
     "S100A8/S100A9 as T-cell-intrinsic dysfunction drivers, at high confidence, independent of model "
     "precision. Model precision (bf16) as an explanation for any finding in this deck."),
])
footer(c, 17)
c.showPage()

# ============================================ 18. DISCUSSION D2+D3 close ===
bg(c)
eyebrow(c, "Discussion · limitations, then where this leaves us", H - 24 * mm)
title(c, "What this study does not establish", H - 40 * mm, 19)
for i, (hd, txt) in enumerate([
    ("Detection-confound axis is an imperfect proxy",
     "Ambient risk is scored partly against the classifier being audited &#8212; evidence of "
     "resemblance to known contaminants, not a direct contamination measurement."),
    ("Single cohort; curated lists",
     "One HTAN/CELLxGENE cohort; the ambient-anchor and T-cell-anchor gene lists are curated, not "
     "ground truth. A different cohort or anchor set could shift specific gene calls."),
    ("Thin normal class, both counts",
     "1 held-out test donor AND 4 cohort donors overall &#8212; already a named limitation in the "
     "repo's own methods (METHODS.md, lines 24-25), carried forward here, not newly flagged."),
    ("Tissue-of-origin skew, per contrast",
     "SCLC: 7/19 (37%) lung, 12/19 metastatic. LUAD: 16/22 (73%) lung. Normal: 4/4 (100%) lung. Skews "
     "sclc_to_normal/normal_to_sclc, not luad_to_normal/normal_to_luad (which benefits from the 3-of-4 "
     "shared-patient pairing above)."),
]):
    cw = (W - 2 * M - 18 * mm) / 4
    x = M + i * (cw + 6 * mm)
    card(c, x, 70 * mm, cw, 60 * mm)
    frame(c, x + 6 * mm, 74 * mm, cw - 12 * mm, 52 * mm, [
        Paragraph(hd, st(10.8, INK, True, 13.6)),
        Paragraph(txt, st(8.6, colors.HexColor("#4A5A55"), False, 11.8)),
    ])
frame(c, M, 30 * mm, W - 2 * M, 32 * mm, [
    Paragraph("Open questions: the original talk's CRISPR TIGIT/HAVCR2 knockout proposal; does NBEAL1 "
              "replicate independently; the talk's own still-open question &#8212; “if SCLC T cells "
              "are not exhausted, what are they, and why does ICI still fail?” A clean negative, "
              "honestly obtained, is a contribution.", st(10.8, colors.HexColor("#31423D"), False, 15)),
])
img(c, ASSETS + "/qr.png", W - M - 22 * mm, 30 * mm, 18 * mm, 18 * mm)
footer(c, 18, "perturbation_workflow/METHODS.md lines 24-25; talk/JSDP_P25_talk.pdf, slide 11")
c.showPage()

c.save()
print(f"wrote {OUT}")
print(f"base branch tip at build time: origin/main = {ORIGIN_MAIN_SHA} ({ORIGIN_MAIN_TIME})")
