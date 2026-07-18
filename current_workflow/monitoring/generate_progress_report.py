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
DISEASE_COMPLETION_DIR = HERE / "disease_completion"
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

SUMMARY_LOG_START = "<!-- JOB_RUN_SUMMARIES_START -->"
SUMMARY_LOG_END = "<!-- JOB_RUN_SUMMARIES_END -->"
SUMMARY_LOG_LIMIT = 48



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


def render_disease_completion_svg(source: SourceProgress, output: Path) -> None:
    colors = {
        "lusc": ("#7c3aed", "#f3e8ff"),
        "luad": ("#2563eb", "#eaf2ff"),
        "normal": ("#059669", "#e8f8f2"),
    }
    accent, pale = colors[source.name]
    label = SOURCE_LABELS[source.name]
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="360" height="190" viewBox="0 0 360 190" role="img" aria-labelledby="title desc">
  <title id="title">{label} perturbation complete: {source.completed_cells:,} of {source.total_cells:,} cells</title>
  <desc id="desc">Final disease completion diagram with exact cell, shard, and deletion counts.</desc>
  <rect width="360" height="190" rx="18" fill="{pale}" stroke="{accent}" stroke-width="2"/>
  <circle cx="58" cy="65" r="34" fill="#fff" stroke="{accent}" stroke-width="4"/>
  <circle cx="48" cy="59" r="10" fill="{accent}" opacity=".55"/>
  <circle cx="69" cy="72" r="8" fill="{accent}" opacity=".75"/>
  <circle cx="64" cy="49" r="5" fill="{accent}" opacity=".35"/>
  <g font-family="-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif">
    <text x="108" y="39" font-size="16" font-weight="700" fill="{accent}">{label} COMPLETE</text>
    <text x="108" y="75" font-size="28" font-weight="800" fill="#17212b">{source.completed_cells:,} / {source.total_cells:,}</text>
    <text x="108" y="98" font-size="15" font-weight="600" fill="#465563">cells · {source.percent:.2f}%</text>
    <rect x="24" y="122" width="312" height="45" rx="10" fill="#fff" opacity=".92"/>
    <text x="42" y="142" font-size="12" fill="#52616d">Completed shards</text>
    <text x="42" y="158" font-size="15" font-weight="700" fill="#17212b">{source.completed_shards:,} / {source.total_shards:,}</text>
    <text x="190" y="142" font-size="12" fill="#52616d">Marker deletions</text>
    <text x="190" y="158" font-size="15" font-weight="700" fill="#17212b">{source.completed_deletions:,} / {source.total_deletions:,}</text>
  </g>
</svg>
"""
    output.write_text(svg)


def disease_completion_gallery(sources: dict[str, SourceProgress]) -> str:
    completed_sources = [
        sources[name]
        for name in SOURCE_ORDER
        if sources[name].completed_cells == sources[name].total_cells
        and sources[name].completed_shards == sources[name].total_shards
    ]
    if not completed_sources:
        return ""

    DISEASE_COMPLETION_DIR.mkdir(parents=True, exist_ok=True)
    cells: list[str] = []
    for source in completed_sources:
        filename = f"{source.name}_complete.svg"
        render_disease_completion_svg(source, DISEASE_COMPLETION_DIR / filename)
        label = SOURCE_LABELS[source.name]
        cells.append(
            '<td align="center" valign="top">'
            f'<img src="disease_completion/{filename}" alt="{label} perturbation complete" width="360"/>'
            "</td>"
        )
    return (
        "\n\n## Disease completion diagrams\n\n"
        "One final diagram is appended for each source after its perturbation "
        "screen completes. Cell totals are shown explicitly.\n\n"
        "<table><tbody><tr>" + "".join(cells) + "</tr></tbody></table>"
    )


def existing_summary_entries(path: Path) -> list[str]:
    if not path.exists():
        return []
    content = path.read_text()
    if SUMMARY_LOG_START not in content or SUMMARY_LOG_END not in content:
        return []
    block = content.split(SUMMARY_LOG_START, 1)[1].split(SUMMARY_LOG_END, 1)[0]
    return [line for line in block.splitlines() if line.startswith("- **")]


def update_summary_entries(
    status: dict[str, Any],
    history: list[dict[str, Any]],
    sources: dict[str, SourceProgress],
) -> list[str]:
    completed, total, percent = overall_progress_from_status(status)
    timestamp = status["timestamp"]
    state = "RUNNING" if status["run_active"] else "IDLE"
    summary = delta_summary(status, history, sources)
    new_entry = (
        f"- **{timestamp}** — {state}; {completed:,} / {total:,} cells "
        f"({percent:.2f}%). {summary}"
    )
    retained = [
        entry
        for entry in existing_summary_entries(REPORT_PATH)
        if not entry.startswith(f"- **{timestamp}**")
    ]
    return [new_entry, *retained][:SUMMARY_LOG_LIMIT]


def build_report(status: dict[str, Any], history: list[dict[str, Any]], generated_at: str) -> str:
    sources = load_sources(status)
    comparisons = load_statistical_comparisons(status, sources)
    completed, total, percent = overall_progress_from_status(status)
    summary = delta_summary(status, history, sources)
    summary_entries = update_summary_entries(status, history, sources)
    summary_log = "\n".join(summary_entries)
    completion_gallery = disease_completion_gallery(sources)

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

> **30-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary. A single final diagram is appended for each disease source when its
> perturbation screen completes.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Job run summaries

Newest refreshes are appended at the top and retained for the most recent
{SUMMARY_LOG_LIMIT} runs.

{SUMMARY_LOG_START}
{summary_log}
{SUMMARY_LOG_END}

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

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

{history_table}

## Job notes

- Scheduled entrypoint: `current_workflow/monitoring/refresh_live_report.sh`
- Render entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Statistics source: `{STATS_DIR}` (override with `PERTURBATION_STATS_DIR`)
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, `cell_interaction_diagram.svg`, and `disease_completion/*.svg`
- Cadence: 30 minutes
{completion_gallery}
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
