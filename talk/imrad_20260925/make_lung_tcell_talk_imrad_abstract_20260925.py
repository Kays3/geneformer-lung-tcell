#!/usr/bin/env python3
"""Generate the IMRaD abstract for the lung T-cell oral presentation.

The word count in the closing note used to be typed by hand and was wrong
three times in a row (197, then 300, then 323 -- the last one truncated by
an off-by-one line range in a throwaway `sed` check, not re-verified before
being written down). Every other number in this deck's slides is pulled
through `git_show()` at build time so it cannot go stale; this script does
the equivalent for the one number in the deck that was still hand-typed --
the abstract computes its own word count from the actual body text at
generation time, so it is correct by construction on every rebuild.
"""
from __future__ import annotations

from pathlib import Path

OUT_PATH = Path(__file__).parent / "lung_tcell_talk_imrad_abstract_20260925.md"

TITLE = (
    "A Foundation-Model T-cell Dysfunction Screen: Translational Candidates, Then a "
    "Full Audit"
)

BACKGROUND = (
    "Foundation models fine-tuned on single-cell transcriptomes can be probed for "
    "candidate disease-state drivers via in-silico perturbation, but a classifier "
    "trained on real tissue also learns whatever else correlates with its labels "
    "— ambient RNA, detection artifacts, cohort composition — and a screen "
    "alone cannot say which one it found."
)

METHODS = (
    "A donor-disjoint SCLC/LUAD/normal T-cell classifier (91.9% accuracy, macro F1 "
    "0.903 -- a third of that macro-F1 from a single held-out normal-class donor, "
    "566 cells; 46,140 cells, 42 donors overall) supported an in-silico "
    "deletion/overexpression screen, scored by significance and by bidirectional "
    "concordance between the two operations. A screen-level ambient-risk audit and "
    "a donor-weighted reanalysis were applied after the screen."
)

RESULTS = (
    "Four originally identified checkpoint/persistence candidates (TIGIT, TIM-3, "
    "CTLA-4, IL7R) replicate and validate independently in spatial tissue data "
    "(antigen-presentation programme, rho=0.361, 7.95 sigma above the null). A "
    "checkpoint-axis ordering the original work had already flagged as an open, "
    "pending question is completed here: donor-level SCLC-vs-LUAD shows no "
    "significant difference (p=.216, two-sided Monte Carlo, 100,000 replicates, 19 "
    "v 22 donors) once a single donor's 74.9% cell share is accounted for. A "
    "subsequent, expanded screen's two most significant hits by FDR (S100A8, "
    "S100A9; FDR~1e-224) fail bidirectional concordance in all six comparisons "
    "tested; screen-wide, half of the top 120 candidates by effect are "
    "ambient-flagged, a curated ambient anchor, or both. One candidate, NBEAL1, "
    "survives every check applied. Model precision (bf16 vs fp32) replicates "
    "cleanly (rho>0.998) and explains none of these results; the pre-registered "
    "104M sign-agreement gate nonetheless still reads FAIL, final, because a "
    "bitwise-deterministic fp32-vs-fp32 rerun made its noise floor zero rather "
    "than a usable sign threshold -- a floor technicality, not evidence that bf16 "
    "disagrees with fp32, and the domain reviewer's ruling left the frozen gate "
    "unamended."
)

DISCUSSION = (
    "A perturbation screen can support real translational candidates and still "
    "require gene-by-gene auditing for confounds; significance and effect size "
    "alone are not sufficient evidence of a cell-intrinsic driver, and an ordering "
    "claim built on a cell-weighted cross-cohort comparison needs a donor-level "
    "check before it is reported either way."
)

SECTIONS = [
    ("Background", BACKGROUND),
    ("Methods", METHODS),
    ("Results", RESULTS),
    ("Discussion", DISCUSSION),
]

DEFAULT_WORD_BUDGET = "250-300"


def word_count(text: str) -> int:
    """Whitespace-token count -- the same method a plain `wc -w` on the
    rendered prose gives, and the one both Stanley and Michael converged on
    independently for the without-section-labels figure."""
    return len(text.split())


def build() -> str:
    body_no_labels = " ".join(text for _, text in SECTIONS)
    n_no_labels = word_count(body_no_labels)
    body_with_labels = " ".join(f"{label}: {text}" for label, text in SECTIONS)
    n_with_labels = word_count(body_with_labels)

    lines = [
        f"# Abstract (structured) — IMRaD restructure, 2026-09-25",
        "",
        f"**Working title:** {TITLE}",
        "",
    ]
    for label, text in SECTIONS:
        lines.append(f"**{label}:** {text}")
        lines.append("")

    note = (
        f"({n_no_labels} words, computed from the section text at generation time "
        f"by this script -- not typed or hand-recounted; {n_with_labels} if the "
        f"four section-label words themselves are also counted. Both figures are "
        f"over the {DEFAULT_WORD_BUDGET}-word default the 2026-09-24 abstract used. "
        f"This count replaces three prior hand-typed figures (197, then 300, then "
        f"323) that were each wrong or stale by the next edit -- the last one "
        f"specifically truncated by an off-by-one line range in a throwaway `sed` "
        f"check that was never re-verified against the file it was checking. Not "
        f"trimmed to force the count under the default -- accuracy set the length; "
        f"that is a consequence for the human's venue/word-limit decision to see, "
        f"not something to compress back out quietly. Author list: Kaisar Dauyey, "
        f"Shinji Nakaoka — matching `poster_final` and `JSDP_P25_talk` exactly; "
        f"no AI-co-authorship line, pending the human's ruling, per the standing "
        f"instruction. Venue, duration and any abstract word limit remain "
        f"unconfirmed with the human.)"
    )
    lines.append(note)
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    content = build()
    OUT_PATH.write_text(content, encoding="utf-8")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
