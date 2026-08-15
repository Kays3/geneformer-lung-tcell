#!/usr/bin/env python3
"""Rebuild the talk deck's image assets from the poster's figures.

`talk/build_deck.js` reads these PNGs with addImage(); it does not create them.
Until now they were hand-made crops of poster figures with no generator, which
is why the deck silently kept showing superseded figures after the poster was
rebuilt - the deck and the board are supposed to read as one project, and there
was nothing keeping them in step.

    python talk/make_assets.py

Two of the assets are whole poster figures, copied and downscaled. Two are
single panels of multi-panel poster figures: rather than pixel-cropping the
composite - which silently breaks whenever the source layout shifts - they are
re-rendered from the same panel functions the poster figure uses, so they are
the same code producing the same marks at a size that suits a 16:9 slide.

The deck reads each PNG's real pixel dimensions and derives its own aspect
ratio (see `pngSize` in build_deck.js), so changing a size here does not
distort a slide.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import qrcode
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
ASSETS = HERE / "assets"
ASSETS.mkdir(parents=True, exist_ok=True)

CKPT = REPO / "sclc_validation" / "checkpoint_cart_perturbation"
PERTURB = REPO / "sclc_validation" / "perturbation_workflow" / "figures"
SPATIAL = REPO / "sclc_validation" / "spatial_validation" / "figures"

GITHUB_URL = "https://github.com/Kays3/geneformer-lung-tcell"

# Slide images are shown a few inches wide on a 13.3 x 7.5 in stage, so there is
# no point carrying poster-resolution pixels into a 3.7 MB pptx.
MAX_WIDTH = 2000


def _load(script: Path):
    """Import a figure generator by path so its panel functions can be reused."""
    spec = importlib.util.spec_from_file_location(script.stem, script)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _copy(source: Path, name: str) -> Path:
    out = ASSETS / name
    image = Image.open(source)
    if image.width > MAX_WIDTH:
        height = round(image.height * MAX_WIDTH / image.width)
        image = image.resize((MAX_WIDTH, height), Image.LANCZOS)
    image.convert("RGB").save(out, optimize=True)
    return out


def _panel(name: str, figsize, draw) -> Path:
    """Render one panel standalone, using the poster figure's own code."""
    fig, ax = plt.subplots(figsize=figsize)
    draw(ax)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.84, bottom=0.20)
    out = ASSETS / name
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> list[Path]:
    written = []

    # ---- whole figures, straight from the poster -------------------------
    written.append(_copy(PERTURB / "sclc_confusion_matrix.png", "confusion.png"))
    written.append(_copy(CKPT / "figures" / "cart_overexpression.png", "detection.png"))
    written.append(_copy(SPATIAL / "spatial_tissue_validation_panel.png", "spatial.png"))

    # ---- single panels, re-rendered rather than cropped ------------------
    ici_module = _load(CKPT / "scripts" / "make_ici_cart_figure.py")
    cart = pd.read_csv(CKPT / "tables" / "cart_engineering_perturbation.csv")
    nodes = pd.read_csv(CKPT / "tables" / "network_node_perturbation.csv")
    edges = pd.read_csv(CKPT / "tables" / "string_network_edges.csv")
    replicated = ici_module._replicated(cart)
    candidates = set(nodes[nodes.is_candidate].gene)

    written.append(_panel(
        "screen_b.png", (5.6, 5.0),
        lambda ax: ici_module._panel_rank(ax, cart, replicated)))
    written.append(_panel(
        "network_c.png", (7.0, 4.8),
        lambda ax: ici_module._panel_network(ax, nodes, edges, candidates)))

    # ---- QR, same target as the poster's ---------------------------------
    qr = qrcode.QRCode(version=2, box_size=8, border=2,
                       error_correction=qrcode.constants.ERROR_CORRECT_M)
    qr.add_data(GITHUB_URL)
    qr.make(fit=True)
    out = ASSETS / "qr.png"
    qr.make_image(fill_color="#10243c", back_color="white").save(out)
    written.append(out)

    return written


if __name__ == "__main__":
    for path in main():
        image = Image.open(path)
        size = path.stat().st_size / 1024
        print(f"  wrote {path.relative_to(REPO)}  {image.width}x{image.height}  {size:.0f} KB")
