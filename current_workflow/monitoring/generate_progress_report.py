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
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from PIL import Image, ImageDraw, ImageFilter, ImageFont


HERE = Path(__file__).resolve().parent
STATUS_PATH = HERE / "latest_status.json"
HISTORY_PATH = HERE / "hourly_history.csv"
REPORT_PATH = HERE / "GPU_PROGRESS_REPORT.md"
PROGRESS_GIF_PATH = HERE / "progress_animation.gif"
PROGRESS_SVG_PATH = HERE / "progress_animation.svg"
DIAGRAM_SVG_PATH = HERE / "cell_interaction_diagram.svg"

SOURCE_ORDER = ("lusc", "luad", "normal")
SOURCE_LABELS = {
    "lusc": "LUSC",
    "luad": "LUAD",
    "normal": "NORMAL",
}


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
    frames: list[Image.Image] = []
    title_font = _font(22, bold=True)
    body_font = _font(16)

    for frame_index in range(18):
        pulse = 10.0 * (1.0 + math.sin((frame_index / 18.0) * 2.0 * math.pi))
        width = min(812.0, base_width + pulse)
        tip_x = 47.0 + width

        canvas = Image.new("RGBA", (900, 150), (255, 255, 255, 255))
        draw = ImageDraw.Draw(canvas)

        draw.text((44, 14), "Held-out all-gene deletion progress", font=title_font, fill=(39, 55, 70, 255))
        _round_rect(draw, (44, 52, 856, 86), 17, (232, 238, 242, 255), (182, 200, 210, 255))

        glow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        glow_draw = ImageDraw.Draw(glow)
        glow_draw.rounded_rectangle((47, 55, tip_x, 83), 14, fill=(120, 255, 180, 120))
        glow_draw.rounded_rectangle((47, 58, tip_x, 80), 12, fill=(32, 216, 120, 140))
        glow = glow.filter(ImageFilter.GaussianBlur(5))
        canvas = Image.alpha_composite(canvas, glow)
        draw = ImageDraw.Draw(canvas)

        core_width = max(0.0, width - 18.0)
        _round_rect(draw, (47, 55, 47 + core_width, 83), 14, (32, 216, 120, 255))
        if core_width > 4:
            _round_rect(draw, (49, 57, 47 + core_width - 2, 81), 12, (185, 255, 213, 190))
        if core_width > 14:
            _round_rect(draw, (52, 59, 47 + core_width - 8, 79), 10, (22, 184, 102, 255))

        tip_glow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
        tip_draw = ImageDraw.Draw(tip_glow)
        tip_draw.ellipse((tip_x - 8, 61, tip_x + 8, 77), fill=(242, 255, 247, 210))
        tip_glow = tip_glow.filter(ImageFilter.GaussianBlur(4))
        canvas = Image.alpha_composite(canvas, tip_glow)
        draw = ImageDraw.Draw(canvas)

        handle_x = max(47.0, tip_x - 16.0)
        _round_rect(draw, (handle_x - 10, 60, handle_x + 6, 78), 4, (12, 143, 78, 230))

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
  <rect x="47" y="55" width="{fill_width:.1f}" height="28" rx="14" fill="url(#blade)" filter="url(#glow)">
    <animate attributeName="width" values="{fill_width:.1f};{animated_width:.1f};{fill_width:.1f}" dur="2.4s" repeatCount="indefinite"/>
  </rect>
  <circle cx="{47 + fill_width:.1f}" cy="69" r="8" fill="#f2fff7" filter="url(#glow)"><animate attributeName="cx" values="{47 + fill_width:.1f};{47 + animated_width:.1f};{47 + fill_width:.1f}" dur="2.4s" repeatCount="indefinite"/></circle>
  <rect x="{37 + fill_width:.1f}" y="60" width="16" height="18" rx="4" fill="#0c8f4e" opacity=".9"><animate attributeName="x" values="{37 + fill_width:.1f};{37 + animated_width:.1f};{37 + fill_width:.1f}" dur="2.4s" repeatCount="indefinite"/></rect>
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


def build_report(status: dict[str, Any], history: list[dict[str, Any]], generated_at: str) -> str:
    sources = load_sources(status)
    completed, total, percent = overall_progress_from_status(status)
    summary = delta_summary(status, history, sources)

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
> summary plus the compact sketch below. The sketch is an orientation aid only;
> it is not a measured ligand-receptor or pathology result.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Current snapshot

**What changed since the prior report:** {summary}

{current_table}

### Progress by source

{source_table}

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

{history_table}

## Job notes

- Job entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, and `cell_interaction_diagram.svg`
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
