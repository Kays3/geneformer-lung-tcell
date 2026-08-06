#!/usr/bin/env python3
"""Generate committable PNG snapshots of large binary deliverables.

Large generated files (poster PDFs, exported HTML) are git-ignored because every
revision would store a fresh multi-megabyte blob. This script renders a
lightweight PNG preview of each so the repository still shows what was produced.

Run from anywhere:

    python tools/make_snapshots.py            # refresh snapshots
    python tools/make_snapshots.py --check    # report staleness, write nothing

Snapshots land in snapshots/ and ARE tracked. The sources they preview are not.

Note: source figures that notebooks read as inputs are never snapshotted or
ignored, even when they exceed the size threshold - degrading or dropping an
input would break reproduction. See .gitignore for the explicit exception.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SNAP_DIR = REPO / "snapshots"

# Long edge of the rendered preview, in pixels. Large enough to read panel
# structure and section headings; small enough to stay well under the limit.
TARGET_LONG_EDGE = 2000

# Files to preview. Globs are relative to the repository root.
# Only GENERATED deliverables belong here - never tracked source such as
# poster/poster_template*.html.
#
# poster/poster_*.pdf is the current A0 poster written by
# poster/generate_poster.ipynb. output/pdf/*.pdf is the earlier 900x2100 mm
# full-board export; both are kept because they are different physical
# deliverables, and the flattened snapshot filenames keep them apart.
SOURCES = [
    "poster/poster_*.pdf",
    "output/pdf/*.pdf",
]

SIZE_LIMIT_MB = 5.0

# Anything smaller than this is committed directly; a preview would be pointless.
MIN_SOURCE_MB = 5.0


def snapshot_name(src: Path) -> str:
    """Flatten a repo-relative path into a single snapshot filename."""
    rel = src.relative_to(REPO)
    return rel.with_suffix(".png").as_posix().replace("/", "__")


def render_pdf(src: Path):
    import pypdfium2 as pdfium

    doc = pdfium.PdfDocument(src)
    try:
        page = doc[0]
        scale = TARGET_LONG_EDGE / max(page.get_width(), page.get_height())
        return page.render(scale=scale).to_pil().convert("RGB")
    finally:
        doc.close()


def render_html(src: Path):
    """HTML previews need a browser engine, which is not always available."""
    raise RuntimeError(
        "HTML rendering needs a headless browser. Export the poster to PDF "
        "first (Chrome > Print > Save as PDF), then re-run this script."
    )


RENDERERS = {".pdf": render_pdf, ".html": render_html}


def discover() -> list[Path]:
    """Generated deliverables over the size threshold, deduplicated."""
    found: list[Path] = []
    for pattern in SOURCES:
        for p in sorted(REPO.glob(pattern)):
            if p.is_file() and p.stat().st_size / 1e6 >= MIN_SOURCE_MB:
                found.append(p)
    return sorted(set(found))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="report missing or stale snapshots, write nothing")
    args = ap.parse_args()

    sources = discover()
    if not sources:
        print("No large deliverables found. Nothing to snapshot.")
        return 0

    SNAP_DIR.mkdir(exist_ok=True)
    stale: list[str] = []
    failures: list[str] = []

    for src in sources:
        rel = src.relative_to(REPO)
        size_mb = src.stat().st_size / 1e6
        out = SNAP_DIR / snapshot_name(src)

        fresh = out.exists() and out.stat().st_mtime >= src.stat().st_mtime
        if args.check:
            state = "ok" if fresh else ("missing" if not out.exists() else "stale")
            if state != "ok":
                stale.append(str(rel))
            print(f"  {state:8s} {rel}  ({size_mb:.1f} MB)")
            continue

        if fresh:
            print(f"  up to date  {rel}")
            continue

        renderer = RENDERERS.get(src.suffix.lower())
        if renderer is None:
            print(f"  SKIP        {rel}  (no renderer for {src.suffix})")
            continue

        try:
            img = renderer(src)
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures.append(f"{rel}: {exc}")
            print(f"  FAILED      {rel}\n              {exc}")
            continue

        buf = io.BytesIO()
        img.save(buf, "PNG", optimize=True)
        data = buf.getvalue()
        out.write_bytes(data)
        print(f"  wrote       {out.relative_to(REPO)}  "
              f"({size_mb:.1f} MB -> {len(data)/1e6:.2f} MB, {img.size[0]}x{img.size[1]})")

        if len(data) / 1e6 > SIZE_LIMIT_MB:
            failures.append(
                f"{out.name} is {len(data)/1e6:.1f} MB, over the {SIZE_LIMIT_MB} MB "
                f"limit - lower TARGET_LONG_EDGE")

    if args.check and stale:
        print(f"\n{len(stale)} snapshot(s) missing or stale. Run without --check.")
        return 1
    if failures:
        print("\nProblems:")
        for f in failures:
            print(f"  - {f}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
