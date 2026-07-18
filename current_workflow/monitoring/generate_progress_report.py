#!/usr/bin/env python3
"""Generate the live perturbation report from monitoring snapshots.

The job is intentionally deterministic:
- read the current snapshot from latest_status.json;
- read the monitor history from hourly_history.csv;
- derive a short delta summary from the latest two samples;
- regenerate the small SVG assets; and
- rewrite GPU_PROGRESS_REPORT.md.

The history filename is legacy. The job treats it as the source of truth for
the recent monitor timeline, regardless of cadence.
"""

from __future__ import annotations

import csv
import html
import json
import math
import os
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from PIL import Image, ImageDraw, ImageFilter, ImageFont


HERE = Path(__file__).resolve().parent
STATUS_PATH = HERE / "latest_status.json"
HISTORY_PATH = HERE / "hourly_history.csv"
REPORT_PATH = HERE / "GPU_PROGRESS_REPORT.md"
SNAPSHOT_GALLERY_DIR = HERE / "snapshot_gallery"
PROGRESS_GIF_PATH = HERE / "progress_animation.gif"
PROGRESS_SVG_PATH = HERE / "progress_animation.svg"
DIAGRAM_SVG_PATH = HERE / "cell_interaction_diagram.svg"
STATS_DIR = Path(os.environ.get("PERTURBATION_STATS_DIR", HERE.parent / "stats"))

SOURCE_ORDER = ("lusc", "luad", "normal")
SOURCE_LABELS = {
    "lusc": "LUSC",
    "luad": "LUAD",
    "normal": "NORMAL",
}
STATISTICAL_COMPARISONS = (
    ("lusc", "normal"),
    ("lusc", "luad"),
    ("luad", "normal"),
    ("luad", "lusc"),
    ("normal", "luad"),
    ("normal", "lusc"),
)



@dataclass(frozen=True)
class SourceProgress:
    name: str
    completed_cells: int
    total_cells: int
    completed_shards: int
    total_shards: int
    completed_deletions: int
    total_deletions: int
    raw_files: int

    @property
    def percent(self) -> float:
        return 100.0 * self.completed_cells / self.total_cells if self.total_cells else 0.0


@dataclass(frozen=True)
class StatisticalComparison:
    source: str
    target: str
    state: str
    result_rows: int | None
    updated_at: str | None
    filename: str

    @property
    def label(self) -> str:
        return f"{SOURCE_LABELS[self.source]} → {SOURCE_LABELS[self.target]}"

    @property
    def complete(self) -> bool:
        return self.state == "COMPLETE"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_history(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    parsed: list[dict[str, Any]] = []
    for row in rows:
        parsed.append(
            {
                "timestamp": row["timestamp"],
                "run_active": row["run_active"].lower() == "true",
                "gpu_utilization_percent": float(row["gpu_utilization_percent"]),
                "gpu_temperature_c": float(row["gpu_temperature_c"]),
                "gpu_power_w": float(row["gpu_power_w"]),
                "gpu_process_memory_mib": float(row["gpu_process_memory_mib"]),
                "system_used_gib": float(row["system_used_gib"]),
                "system_available_gib": float(row["system_available_gib"]),
                "lusc_cells": int(row["lusc_cells"]),
                "luad_cells": int(row["luad_cells"]),
                "normal_cells": int(row["normal_cells"]),
                "completed_shards": int(row["completed_shards"]),
            }
        )
    return parsed


def overall_progress_from_status(status: dict[str, Any]) -> tuple[int, int, float]:
    completed = 0
    total = 0
    for key in SOURCE_ORDER:
        source = status["progress"][key]
        completed += int(source["completed_cells"])
        total += int(source["total_cells"])
    percent = 100.0 * completed / total if total else 0.0
    return completed, total, percent


def overall_progress_from_history_row(row: dict[str, Any], total: int) -> tuple[int, float]:
    completed = int(row["lusc_cells"]) + int(row["luad_cells"]) + int(row["normal_cells"])
    percent = 100.0 * completed / total if total else 0.0
    return completed, percent


def load_sources(status: dict[str, Any]) -> dict[str, SourceProgress]:
    result: dict[str, SourceProgress] = {}
    for name in SOURCE_ORDER:
        source = status["progress"][name]
        result[name] = SourceProgress(
            name=name,
            completed_cells=int(source["completed_cells"]),
            total_cells=int(source["total_cells"]),
            completed_shards=int(source["completed_shards"]),
            total_shards=int(source["total_shards"]),
            completed_deletions=int(source["completed_marker_deletions"]),
            total_deletions=int(source["total_deletions"]),
            raw_files=int(source["raw_files"]),
        )
    return result


def csv_result_rows(path: Path) -> int:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        return sum(1 for row in reader if row)


def load_statistical_comparisons(
    status: dict[str, Any], sources: dict[str, SourceProgress]
) -> list[StatisticalComparison]:
    """Inspect final comparison artifacts without interpreting their biology."""
    perturbation_complete = all(
        source.completed_shards == source.total_shards for source in sources.values()
    )
    comparisons: list[StatisticalComparison] = []
    for source, target in STATISTICAL_COMPARISONS:
        filename = f"heldout_allgene_{source}_to_{target}.csv"
        path = STATS_DIR / filename
        if path.is_file():
            try:
                result_rows = csv_result_rows(path)
                state = "COMPLETE" if result_rows > 0 else "EMPTY OUTPUT"
            except (OSError, csv.Error):
                result_rows = None
                state = "UNREADABLE OUTPUT"
            updated_at = datetime.fromtimestamp(
                path.stat().st_mtime
            ).astimezone().isoformat(timespec="seconds")
        else:
            result_rows = None
            updated_at = None
            if perturbation_complete:
                state = "QUEUED / RUNNING" if status["run_active"] else "PENDING"
            else:
                state = "WAITING FOR PERTURBATION"
        comparisons.append(
            StatisticalComparison(
                source=source,
                target=target,
                state=state,
                result_rows=result_rows,
                updated_at=updated_at,
                filename=filename,
            )
        )
    return comparisons


def statistical_comparison_summary(
    comparisons: list[StatisticalComparison], sources: dict[str, SourceProgress]
) -> str:
    complete = sum(comparison.complete for comparison in comparisons)
    total = len(comparisons)
    perturbation_complete = all(
        source.completed_shards == source.total_shards for source in sources.values()
    )
    if complete == total:
        return (
            f"All {total} directional comparisons have generated non-empty result "
            "tables. Statistical outputs are available for review."
        )
    if not perturbation_complete:
        completed_shards = sum(source.completed_shards for source in sources.values())
        total_shards = sum(source.total_shards for source in sources.values())
        return (
            f"{complete} / {total} comparisons are complete. Final aggregation waits "
            f"for the deletion screen, currently {completed_shards:,} / "
            f"{total_shards:,} shards."
        )
    return (
        f"{complete} / {total} comparisons are complete. The remaining comparison "
        "jobs are queued or running."
    )


def delta_summary(
    status: dict[str, Any],
    history: list[dict[str, Any]],
    sources: dict[str, SourceProgress],
) -> str:
    completed, total, percent = overall_progress_from_status(status)
    if len(history) < 2:
        return (
            f"The run is at {completed:,} / {total:,} cells ({percent:.2f}%). "
            "This is the first generated snapshot."
        )

    previous = history[-2]
    previous_completed, previous_percent = overall_progress_from_history_row(previous, total)
    shard_delta = sources["lusc"].completed_shards + sources["luad"].completed_shards + sources["normal"].completed_shards - previous["completed_shards"]
    cell_delta = completed - previous_completed

    clauses: list[str] = [
        (
            f"The run advanced by {cell_delta:,} cells and {shard_delta:,} shards, "
            f"lifting completion from {previous_percent:.2f}% to {percent:.2f}%."
        )
    ]

    previous_source_cells = {
        "lusc": int(previous["lusc_cells"]),
        "luad": int(previous["luad_cells"]),
        "normal": int(previous["normal_cells"]),
    }

    source_clauses: list[str] = []
    for name in SOURCE_ORDER:
        current = sources[name]
        previous_cells = previous_source_cells[name]
        label = SOURCE_LABELS[name]
        if current.completed_cells == current.total_cells and previous_cells == current.total_cells:
            source_clauses.append(f"{label} remained complete")
        elif current.completed_cells == current.total_cells and previous_cells < current.total_cells:
            source_clauses.append(f"{label} finished")
        elif current.completed_cells == previous_cells:
            if current.completed_cells == 0:
                source_clauses.append(f"{label} has not started")
            else:
                source_clauses.append(f"{label} held at {current.completed_cells:,} cells")
        else:
            source_clauses.append(
                f"{label} moved from {previous_cells:,} to {current.completed_cells:,} cells"
            )

    gpu_now = float(status["gpu"]["utilization_percent"])
    gpu_prev = float(previous["gpu_utilization_percent"])
    if gpu_now != gpu_prev:
        direction = "rose" if gpu_now > gpu_prev else "fell"
        source_clauses.append(
            f"GPU utilization {direction} from {gpu_prev:.0f}% to {gpu_now:.0f}%"
        )
    else:
        source_clauses.append(f"GPU utilization stayed at {gpu_now:.0f}%")

    clauses.append("; ".join(source_clauses) + ".")
    return " ".join(clauses)


def fmt_int(value: int | float) -> str:
    return f"{int(round(value)):,}"


def fmt_float(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}f}"


def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/dejavu/DejaVuSans.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def _round_rect(draw: ImageDraw.ImageDraw, box: tuple[float, float, float, float], radius: float, fill: tuple[int, int, int, int], outline: tuple[int, int, int, int] | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def render_progress_gif(status: dict[str, Any], output: Path) -> None:
    completed, total, percent = overall_progress_from_status(status)
    base_width = 812 * percent / 100.0
    handle_left = 47
    handle_width = 34
    blade_left = handle_left + handle_width + 8
    blade_track = 812 - handle_width - 8
    frames: list[Image.Image] = []
    title_font = _font(22, bold=True)
    body_font = _font(16)

    for frame_index in range(18):
        pulse = 8.0 * (1.0 + math.sin((frame_index / 18.0) * 2.0 * math.pi))
        width = min(float(blade_track), base_width + pulse)
        tip_x = blade_left + width

        canvas = Image.new("RGBA", (900, 150), (255, 255, 255, 255))
        draw = ImageDraw.Draw(canvas)

        draw.text((44, 14), "Held-out all-gene deletion progress", font=title_font, fill=(39, 55, 70, 255))
        _round_rect(draw, (44, 52, 856, 86), 17, (232, 238, 242, 255), (182, 200, 210, 255))

        glow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        glow_draw = ImageDraw.Draw(glow)
        glow_draw.rounded_rectangle((blade_left, 55, tip_x, 83), 14, fill=(120, 255, 180, 120))
        glow_draw.rounded_rectangle((blade_left, 58, tip_x, 80), 12, fill=(32, 216, 120, 140))
        glow = glow.filter(ImageFilter.GaussianBlur(5))
        canvas = Image.alpha_composite(canvas, glow)
        draw = ImageDraw.Draw(canvas)

        if width > 0:
            _round_rect(draw, (blade_left, 55, tip_x, 83), 14, (32, 216, 120, 255))
            if width > 4:
                _round_rect(draw, (blade_left + 2, 57, tip_x - 2, 81), 12, (185, 255, 213, 190))
            if width > 14:
                _round_rect(draw, (blade_left + 5, 59, tip_x - 8, 79), 10, (22, 184, 102, 255))

        tip_glow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        tip_draw = ImageDraw.Draw(tip_glow)
        tip_draw.ellipse((tip_x - 8, 61, tip_x + 8, 77), fill=(242, 255, 247, 210))
        tip_glow = tip_glow.filter(ImageFilter.GaussianBlur(4))
        canvas = Image.alpha_composite(canvas, tip_glow)
        draw = ImageDraw.Draw(canvas)

        _round_rect(draw, (handle_left, 58, handle_left + handle_width, 80), 5, (37, 44, 52, 255), (118, 129, 137, 255), 1)
        _round_rect(draw, (handle_left + 4, 60, handle_left + 9, 78), 2, (219, 227, 233, 255))
        _round_rect(draw, (handle_left + 11, 59, handle_left + 24, 79), 3, (57, 68, 77, 255))
        _round_rect(draw, (handle_left + 20, 61, handle_left + 27, 77), 2, (156, 167, 176, 255))
        draw.line((handle_left + handle_width + 4, 69, blade_left + 4, 69), fill=(255, 255, 255, 150), width=2)

        draw.text(
            (44, 111),
            f"{completed:,} / {total:,} cells  |  {percent:.2f}%  |  {'RUNNING' if status['run_active'] else 'IDLE'}",
            font=body_font,
            fill=(66, 91, 112, 255),
        )

        frames.append(canvas.convert("P", palette=Image.ADAPTIVE))

    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=120,
        loop=0,
        disposal=2,
        optimize=False,
    )


def render_progress_svg(status: dict[str, Any], output: Path) -> None:
    completed, total, percent = overall_progress_from_status(status)
    fill_width = 812 * percent / 100.0
    animated_width = min(812, fill_width + 18)
    handle_left = 47
    handle_width = 34
    blade_left = handle_left + handle_width + 8
    blade_track = 812 - handle_width - 8
    blade_width = max(0.0, min(blade_track, fill_width - handle_width - 4))
    animated_blade_width = max(0.0, min(blade_track, animated_width - handle_width - 4))
    label = (
        f"{completed:,} / {total:,} cells  |  {percent:.2f}%  |  "
        f"{'RUNNING' if status['run_active'] else 'IDLE'}"
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="900" height="150" viewBox="0 0 900 150" role="img" aria-labelledby="title desc">
  <title id="title">Held-out all-gene deletion progress</title>
  <desc id="desc">A green lightsaber-like progress bar showing {percent:.2f} percent completion.</desc>
  <defs>
    <linearGradient id="track" x1="0" x2="1"><stop stop-color="#e8eef2"/><stop offset="1" stop-color="#dce6ea"/></linearGradient>
    <linearGradient id="blade" x1="0" x2="1"><stop stop-color="#20d878"/><stop offset=".5" stop-color="#b9ffd5"/><stop offset="1" stop-color="#16b866"/></linearGradient>
    <filter id="glow" x="-20%" y="-100%" width="140%" height="300%"><feGaussianBlur stdDeviation="5" result="blur"/><feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge></filter>
  </defs>
  <rect width="900" height="150" fill="#fff"/>
  <text x="44" y="31" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="22" font-weight="700" fill="#273746">Held-out all-gene deletion progress</text>
  <rect x="44" y="52" width="812" height="34" rx="17" fill="url(#track)" stroke="#b6c8d2"/>
  <rect x="{blade_left}" y="55" width="{blade_width:.1f}" height="28" rx="14" fill="url(#blade)" filter="url(#glow)">
    <animate attributeName="width" values="{blade_width:.1f};{animated_blade_width:.1f};{blade_width:.1f}" dur="2.4s" repeatCount="indefinite"/>
  </rect>
  <circle cx="{blade_left + blade_width:.1f}" cy="69" r="8" fill="#f2fff7" filter="url(#glow)"><animate attributeName="cx" values="{blade_left + blade_width:.1f};{blade_left + animated_blade_width:.1f};{blade_left + blade_width:.1f}" dur="2.4s" repeatCount="indefinite"/></circle>
  <rect x="{handle_left}" y="58" width="{handle_width}" height="22" rx="5" fill="#252c34" stroke="#7d8790" stroke-width="1"/>
  <rect x="{handle_left + 4}" y="60" width="5" height="18" rx="2" fill="#dfe6ea"/>
  <rect x="{handle_left + 11}" y="59" width="13" height="20" rx="3" fill="#3f4a55"/>
  <rect x="{handle_left + 20}" y="61" width="7" height="16" rx="2" fill="#a7b0b8"/>
  <rect x="{handle_left + handle_width + 4}" y="67" width="4" height="4" rx="2" fill="#ffffff" opacity=".75"/>
  <text x="44" y="119" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="16" fill="#425b70">{label}</text>
</svg>
"""
    output.write_text(svg)


def render_diagram_svg(output: Path) -> None:
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="760" height="210" viewBox="0 0 760 210" role="img" aria-labelledby="title desc">
  <title id="title">Single-cell-inspired interaction sketch</title>
  <desc id="desc">A shotgun-pathology-inspired compact diagram connecting a T cell, tumor state, and reference centroids.</desc>
  <defs><marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#25b86b"/></marker></defs>
  <rect width="760" height="210" rx="16" fill="#f7fbf9" stroke="#c7dfd1"/>
  <g font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif">
    <text x="22" y="28" font-size="13" font-weight="700" fill="#35634a" letter-spacing="1">SINGLE-CELL INTERACTION SNAPSHOT</text>
    <g transform="translate(82 106)"><circle r="47" fill="#dff7e8" stroke="#29b86a" stroke-width="3"/><circle cx="-13" cy="-5" r="15" fill="#9edbb6"/><circle cx="18" cy="10" r="10" fill="#72c999"/><circle cx="1" cy="-18" r="5" fill="#fff"/><text x="0" y="70" text-anchor="middle" font-size="13" fill="#24563c">T cell</text></g>
    <path d="M135 106 C190 68 210 68 263 106" fill="none" stroke="#25b86b" stroke-width="3" marker-end="url(#arrow)"/>
    <g transform="translate(320 106)"><circle r="54" fill="#fff1e7" stroke="#ef9961" stroke-width="3"/><path d="M-34-7 Q-10-40 16-18 T34 14 Q2 43-31 20Z" fill="#f6c39d"/><circle cx="-12" cy="0" r="7" fill="#d77b4e"/><circle cx="16" cy="-7" r="5" fill="#d77b4e"/><text x="0" y="77" text-anchor="middle" font-size="13" fill="#70462f">Shotgun tissue field</text></g>
    <path d="M378 106 C431 68 451 68 503 106" fill="none" stroke="#25b86b" stroke-width="3" marker-end="url(#arrow)"/>
    <g transform="translate(580 75)"><circle r="31" fill="#e4f0ff" stroke="#5f98d5" stroke-width="2"/><text text-anchor="middle" y="5" font-size="12" fill="#2c5684">LUAD</text></g>
    <g transform="translate(580 143)"><circle r="31" fill="#f0e8ff" stroke="#9d78d2" stroke-width="2"/><text text-anchor="middle" y="5" font-size="12" fill="#5a4084">LUSC</text></g>
    <path d="M520 93 L547 81 M520 119 L547 137" stroke="#7b91a4" stroke-width="2" stroke-dasharray="4 4" marker-end="url(#arrow)"/>
    <text x="680" y="102" text-anchor="middle" font-size="11" fill="#577064">model state</text><text x="680" y="118" text-anchor="middle" font-size="11" fill="#577064">reference shift</text>
  </g>
</svg>
"""
    output.write_text(svg)


def markdown_table(headers: Iterable[str], rows: Iterable[Iterable[str]]) -> str:
    header_row = "| " + " | ".join(headers) + " |"
    separator = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([header_row, separator, *body])


def progress_blocks(percent: float, width: int = 12) -> str:
    filled = max(0, min(width, int(round((percent / 100.0) * width))))
    return "█" * filled + "░" * (width - filled)


def snapshot_slug(timestamp: str) -> str:
    return datetime.fromisoformat(timestamp).strftime("%Y%m%dT%H%M%S%z")


def snapshot_path(row: dict[str, Any]) -> Path:
    return SNAPSHOT_GALLERY_DIR / f"snapshot_{snapshot_slug(row['timestamp'])}.svg"


def render_snapshot_thumbnail_svg(row: dict[str, Any], total_cells: int, output: Path) -> None:
    cells, percent = overall_progress_from_history_row(row, total_cells)
    lusc = int(row["lusc_cells"])
    luad = int(row["luad_cells"])
    normal = int(row["normal_cells"])
    source_total = max(1, lusc + luad + normal)
    total_shards = max(1, int(row["completed_shards"]))
    active = "RUNNING" if row["run_active"] else "IDLE"
    stamp = html.escape(datetime.fromisoformat(row["timestamp"]).strftime("%Y-%m-%d %H:%M"))
    timestamp = html.escape(row["timestamp"])
    title = (
        f"Snapshot {timestamp} | {percent:.2f}% complete | "
        f"{cells:,} cells | {active}"
    )
    r = 36
    circumference = 2 * math.pi * r
    dash = max(0.0, min(circumference, circumference * percent / 100.0))

    def source_badge(x: int, y: int, label: str, value: int, total: int, fill: str, pale: str) -> str:
        pct = 100.0 * value / total if total else 0.0
        return f"""<g transform="translate({x} {y})">
      <rect x="0" y="0" width="122" height="30" rx="10" fill="{pale}" stroke="#d8e3ea"/>
      <circle cx="17" cy="15" r="8" fill="{fill}"/>
      <text x="32" y="13" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="11" font-weight="700" fill="#24313d">{label}</text>
      <text x="32" y="24" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="10" fill="#51606d">{value:,} / {total:,} ({pct:.1f}%)</text>
    </g>"""

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="260" height="180" viewBox="0 0 260 180" role="img" aria-labelledby="title desc">
  <title id="title">{title}</title>
  <desc id="desc">A compact single-cell diagram for the {timestamp} snapshot, showing overall job progress and NSCLC pathology state.</desc>
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#ffffff"/>
      <stop offset="1" stop-color="#f7fbf9"/>
    </linearGradient>
    <linearGradient id="cell" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="#e7f9ef"/>
      <stop offset="1" stop-color="#d7f2e3"/>
    </linearGradient>
    <linearGradient id="progress" x1="0" x2="1">
      <stop offset="0" stop-color="#20d878"/>
      <stop offset="1" stop-color="#0ea55b"/>
    </linearGradient>
  </defs>
  <rect x="0" y="0" width="260" height="180" rx="16" fill="url(#bg)" stroke="#ceddde"/>
  <text x="14" y="20" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="11" font-weight="700" fill="#365368">{stamp}</text>
  <text x="14" y="34" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="10" fill="#607381">{cells:,} cells · {percent:.2f}% · {active}</text>
  <g transform="translate(68 88)">
    <circle r="41" fill="none" stroke="#dce7eb" stroke-width="10"/>
    <circle r="41" fill="none" stroke="url(#progress)" stroke-width="10" stroke-linecap="round" transform="rotate(-90)" stroke-dasharray="{dash:.2f} {circumference - dash:.2f}"/>
    <circle r="31" fill="url(#cell)" stroke="#2bb56f" stroke-width="2"/>
    <circle cx="-9" cy="-6" r="10" fill="#9edcb6"/>
    <circle cx="11" cy="8" r="7" fill="#73c898"/>
    <circle cx="0" cy="-15" r="4" fill="#ffffff"/>
    <circle cx="17" cy="-11" r="2.5" fill="#ffffff" opacity=".85"/>
    <circle cx="-16" cy="13" r="2.5" fill="#ffffff" opacity=".85"/>
    <text x="0" y="50" text-anchor="middle" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="10" font-weight="700" fill="#2c5a3f">job cell</text>
  </g>
  <rect x="14" y="138" width="94" height="28" rx="10" fill="#eef7ff" stroke="#d7e8fb"/>
  <text x="20" y="150" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="9" font-weight="700" fill="#2c5684">Progress</text>
  <text x="20" y="160" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="10" fill="#2c5684">{progress_blocks(percent, 14)}</text>
  <text x="93" y="157" text-anchor="end" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="9" fill="#2c5684">{percent:.2f}%</text>
  {source_badge(112, 56, "LUAD", luad, source_total, "#4f83c6", "#ebf3ff")}
  {source_badge(112, 88, "LUSC", lusc, source_total, "#a36ad9", "#f4ecff")}
  {source_badge(112, 120, "NORMAL", normal, source_total, "#43b87a", "#ecfbf3")}
  <text x="200" y="156" text-anchor="middle" font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif" font-size="9" fill="#5f7180">shards {total_shards:,}</text>
</svg>
"""
    output.write_text(svg)


def render_snapshot_gallery(history: list[dict[str, Any]], total_cells: int) -> list[dict[str, str]]:
    SNAPSHOT_GALLERY_DIR.mkdir(parents=True, exist_ok=True)
    gallery: list[dict[str, str]] = []
    for row in history:
        path = snapshot_path(row)
        render_snapshot_thumbnail_svg(row, total_cells, path)
        cells, percent = overall_progress_from_history_row(row, total_cells)
        gallery.append(
            {
                "path": f"snapshot_gallery/{path.name}",
                "timestamp": row["timestamp"],
                "cells": f"{cells:,}",
                "percent": f"{percent:.2f}",
                "gpu": f"{float(row['gpu_utilization_percent']):.0f}%",
                "temp": f"{float(row['gpu_temperature_c']):.0f} C",
                "power": f"{float(row['gpu_power_w']):.1f} W",
                "shards": f"{int(row['completed_shards']):,}",
                "run_state": "RUNNING" if row["run_active"] else "IDLE",
                "slug": snapshot_slug(row["timestamp"]),
            }
        )
    return gallery


def snapshot_gallery(gallery: list[dict[str, str]], columns: int = 3, limit: int = 6) -> str:
    recent = list(reversed(gallery[-limit:]))
    if not recent:
        return "_No snapshot history available._"

    rows: list[str] = []
    for start in range(0, len(recent), columns):
        chunk = recent[start : start + columns]
        cells = []
        for item in chunk:
            cells.append(
                "<td align=\"center\" valign=\"top\">"
                f"<img src=\"{item['path']}\" alt=\"Snapshot {html.escape(item['timestamp'])}\" width=\"260\"/>"
                "<br/>"
                f"<sub>{html.escape(item['timestamp'])} · {item['percent']}% · {item['cells']} cells</sub>"
                "</td>"
            )
        while len(cells) < columns:
            cells.append("<td></td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")

    header = "".join("<th align=\"center\">Single-cell snapshot</th>" for _ in range(columns))
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def build_report(status: dict[str, Any], history: list[dict[str, Any]], generated_at: str) -> str:
    sources = load_sources(status)
    comparisons = load_statistical_comparisons(status, sources)
    completed, total, percent = overall_progress_from_status(status)
    summary = delta_summary(status, history, sources)
    gallery_items = render_snapshot_gallery(history, total)
    gallery = snapshot_gallery(gallery_items)

    rows = [
        ("Generated", generated_at),
        ("Run status", "RUNNING" if status["run_active"] else "IDLE"),
        ("Overall cell progress", f"{completed:,} / {total:,} ({percent:.2f}%)"),
        ("GPU", status["gpu"]["name"]),
        ("GPU utilization", f"{fmt_float(status['gpu']['utilization_percent'], 0)}%"),
        ("GPU temperature", f"{fmt_float(status['gpu']['temperature_c'], 0)} C"),
        ("GPU power", f"{fmt_float(status['gpu']['power_w'], 1)} W"),
        ("Perturbation GPU memory", f"{fmt_int(status['gpu']['perturbation_process']['memory_mib'])} MiB"),
        ("System memory used", f"{fmt_float(status['memory']['used_gib'], 1)} GiB"),
        ("System memory available", f"{fmt_float(status['memory']['available_gib'], 1)} GiB"),
        ("Swap used", f"{fmt_float(status['memory']['swap_used_gib'], 1)} GiB"),
    ]
    current_table = markdown_table(("Metric", "Value"), rows)

    source_rows = []
    for name in SOURCE_ORDER:
        source = sources[name]
        source_rows.append(
            (
                SOURCE_LABELS[name],
                f"{source.completed_cells:,} / {source.total_cells:,} ({source.percent:.2f}%)",
                f"{source.completed_shards:,} / {source.total_shards:,}",
                f"{source.raw_files:,}",
                f"{source.completed_deletions:,} / {source.total_deletions:,}",
            )
        )
    source_table = markdown_table(
        ("Source", "Cells", "Shards", "Raw files", "Marker deletions"),
        source_rows,
    )

    comparison_rows = []
    for comparison in comparisons:
        comparison_rows.append(
            (
                comparison.label,
                comparison.state,
                f"{comparison.result_rows:,}" if comparison.result_rows is not None else "—",
                comparison.updated_at or "—",
                f"`{comparison.filename}`",
            )
        )
    comparison_table = markdown_table(
        ("Comparison", "State", "Result rows", "Updated", "Output"),
        comparison_rows,
    )
    comparison_summary = statistical_comparison_summary(comparisons, sources)

    history_rows = []
    for row in reversed(history[-8:]):
        overall_cells, overall_pct = overall_progress_from_history_row(row, total)
        history_rows.append(
            (
                row["timestamp"],
                f"{overall_cells:,}",
                f"{overall_pct:.2f}%",
                f"{row['gpu_utilization_percent']:.0f}%",
                f"{row['gpu_temperature_c']:.0f} C",
                f"{row['gpu_power_w']:.1f} W",
                f"{row['completed_shards']:,}",
            )
        )
    history_table = markdown_table(
        ("Timestamp", "Cells", "Progress", "GPU util", "Temp", "Power", "Shards"),
        history_rows,
    )

    return f"""# Held-out all-gene deletion status

![Animated overall cell progress](progress_animation.gif)

> **15-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary, and the gallery below shows the recent single-cell snapshots. The
> thumbnails are orientation aids only; they are not measured ligand-receptor
> or pathology results.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Current snapshot

**What changed since the prior report:** {summary}

{current_table}

### Progress by source

{source_table}

## Final statistical comparisons

**{comparison_summary}**

{comparison_table}

Result-row counts confirm artifact generation only; they do not establish
biological significance. Gene rankings should be interpreted only after all
six comparisons complete and coverage, FDR, and donor-consistency checks pass.

## Snapshot gallery

Each thumbnail is a single-cell diagram that encodes overall progress and the
LUAD, LUSC, and NORMAL pathology balance at that snapshot.

{gallery}

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

{history_table}

## Job notes

- Job entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Statistics source: `{STATS_DIR}` (override with `PERTURBATION_STATS_DIR`)
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, `cell_interaction_diagram.svg`, and `snapshot_gallery/*.svg`
- Cadence: 15 minutes
"""


def write_atomic(path: Path, content: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(content)
    tmp.replace(path)


def main() -> None:
    if not STATUS_PATH.exists():
        raise SystemExit(f"Missing status file: {STATUS_PATH}")
    if not HISTORY_PATH.exists():
        raise SystemExit(f"Missing history file: {HISTORY_PATH}")

    status = read_json(STATUS_PATH)
    history = read_history(HISTORY_PATH)
    generated_at = status["timestamp"]

    render_progress_gif(status, PROGRESS_GIF_PATH)
    render_progress_svg(status, PROGRESS_SVG_PATH)
    render_diagram_svg(DIAGRAM_SVG_PATH)
    report = build_report(status, history, generated_at)
    write_atomic(REPORT_PATH, report)


if __name__ == "__main__":
    main()
