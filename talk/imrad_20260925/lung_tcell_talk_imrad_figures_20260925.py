"""Figure staging for the IMRaD-restructured oral-presentation deck, 2026-09-25.

Two kinds of figure in this deck:

1. Four figures already built and independently verified (twice) for the
   2026-09-24 deck (fig_concordance, fig_ambient_breakdown, fig_t6_reversal,
   fig_bf16_canary). They are reused UNMODIFIED here -- referenced directly
   from lung_tcell_talk_assets_20260924/ by the slide generator, not copied
   or regenerated. Their own generator (lung_tcell_talk_figures_20260924.py)
   remains the source of truth for how they were built.

2. Six original JSDP figures (confusion, screen_b, detection, spatial,
   network_c, qr) restored as Results/Introduction figures now that the
   original results are back in as results. talk/assets/*.png is NOT
   tracked in git (talk/ is gitignored, same as talk/JSDP_P25_talk.pptx/pdf
   and poster/poster_final.*) -- these are local filesystem artifacts, not
   git-citable objects, so they are staged here by a verified filesystem
   copy (size-checked), the same treatment the sources note already gives
   JSDP_P25_talk and poster_final themselves.
"""
from __future__ import annotations

import hashlib
import os
import shutil

REPO_ASSETS = "/Users/kaisar/workspace/github/geneformer-lung-tcell/talk/assets"
OUT_DIR = "/Users/kaisar/workspace/office_hive/hive/reports/lung_tcell_talk_imrad_assets_20260925"
os.makedirs(OUT_DIR, exist_ok=True)

JSDP_ASSETS = {
    "confusion.png": 87715,
    "screen_b.png": 80663,
    "detection.png": 308598,
    "spatial.png": 2520539,
    "network_c.png": 163632,
    "qr.png": 1598,
}

REUSED_20260924 = "/Users/kaisar/workspace/office_hive/hive/reports/lung_tcell_talk_assets_20260924"
REUSED_FIGS = ["fig_concordance.png", "fig_ambient_breakdown.png", "fig_t6_reversal.png", "fig_bf16_canary.png"]


def stage_jsdp_assets() -> dict[str, str]:
    hashes = {}
    for name, expected_size in JSDP_ASSETS.items():
        src = os.path.join(REPO_ASSETS, name)
        actual_size = os.path.getsize(src)
        assert actual_size == expected_size, (
            f"{name}: expected {expected_size} bytes (recorded 2026-09-24), got {actual_size} -- "
            "the JSDP asset changed since it was last measured; do not silently reuse a different file."
        )
        dst = os.path.join(OUT_DIR, name)
        shutil.copyfile(src, dst)
        hashes[name] = hashlib.sha256(open(dst, "rb").read()).hexdigest()[:16]
    return hashes


def verify_reused_figures_present() -> None:
    for name in REUSED_FIGS:
        path = os.path.join(REUSED_20260924, name)
        assert os.path.exists(path), f"{name} missing from the 2026-09-24 assets dir -- was it moved or deleted?"


if __name__ == "__main__":
    h = stage_jsdp_assets()
    verify_reused_figures_present()
    for name, digest in h.items():
        print(f"staged {name}: sha256[:16]={digest}")
    print(f"confirmed present (referenced, not copied): {REUSED_FIGS}")
