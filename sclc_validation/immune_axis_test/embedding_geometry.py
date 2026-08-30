#!/usr/bin/env python3
"""T3: measure centroid geometry and reconstruct in-plane OE responses.

Geneformer reports the change in cosine similarity to each state centroid. The
committed tables retain the two non-source-state shifts. Given the centroids,
this script reconstructs the unique vector in the affine centroid direction
plane that reproduces those scores.

That vector is a model-based in-plane representation, not an observed embedding
displacement: two scores cannot identify an out-of-plane component. Centroid
geometry itself is measured without that reconstruction assumption.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
DEFAULT_CENTROIDS = RESULTS / "training_donor_disease_centroids.npz"
DEFAULT_SHIFTS = (
    HERE.parent
    / "perturbation_workflow/targeted_panel/results/targeted_panel_delete_overexpress_merged.csv"
)

STATES = ("normal", "sclc", "luad")
STATE_LABELS = {"normal": "Normal", "sclc": "SCLC", "luad": "LUAD"}
STATE_COLORS = {"normal": "#6f7f8f", "sclc": "#b34f6a", "luad": "#157a79"}
PROGRAMS = {
    "exhaustion": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"],
    "cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG"],
    "progenitor": ["TCF7", "SLAMF6", "IL7R", "CCR7"],
    "sclc_subtype_tf": ["ASCL1", "NEUROD1", "POU2F3", "YAP1"],
}
COMPLEMENT = (
    ("sclc", "luad", "normal"),
    ("normal", "luad", "sclc"),
    ("luad", "normal", "sclc"),
)


class GeometryError(ValueError):
    """Raised when the centroids cannot support a 2-D reconstruction."""


def load_centroids(path: Path) -> dict[str, np.ndarray]:
    if path.suffix != ".npz":
        raise ValueError("T3 accepts .npz centroids; use export_state_centroids.py first")
    with np.load(path, allow_pickle=False) as payload:
        missing = set(STATES) - set(payload.files)
        if missing:
            raise ValueError(f"Centroid file is missing {sorted(missing)}")
        centroids = {state: np.asarray(payload[state], dtype=float).squeeze() for state in STATES}
    shapes = {value.shape for value in centroids.values()}
    if len(shapes) != 1 or next(iter(shapes)) == ():
        raise ValueError(f"Centroids must be equal-length vectors; got {sorted(shapes)}")
    for state, vector in centroids.items():
        if vector.ndim != 1 or not np.isfinite(vector).all() or np.linalg.norm(vector) == 0:
            raise ValueError(f"Invalid centroid for {state}: shape={vector.shape}")
    return centroids


def normalized(centroids: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {state: vector / np.linalg.norm(vector) for state, vector in centroids.items()}


def centroid_plane(
    unit: dict[str, np.ndarray], tolerance: float = 1e-10
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Return origin, oriented basis, and coordinates for the affine plane."""
    origin = unit["normal"]
    normal_to_luad = unit["luad"] - origin
    e1_norm = np.linalg.norm(normal_to_luad)
    if e1_norm <= tolerance:
        raise GeometryError("Normal and LUAD centroids coincide after cosine normalization")
    e1 = normal_to_luad / e1_norm
    normal_to_sclc = unit["sclc"] - origin
    residual = normal_to_sclc - np.dot(normal_to_sclc, e1) * e1
    e2_norm = np.linalg.norm(residual)
    if e2_norm <= tolerance:
        raise GeometryError("The three normalized centroids are collinear; the response plane is undefined")
    e2 = residual / e2_norm
    basis = np.column_stack((e1, e2))
    coords = {state: (vector - origin) @ basis for state, vector in unit.items()}
    return origin, basis, coords


def geometry_tables(
    raw: dict[str, np.ndarray],
    unit: dict[str, np.ndarray],
    basis: np.ndarray,
    coords: dict[str, np.ndarray],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    state_rows = []
    for state in STATES:
        perpendicular = unit[state] - basis @ (basis.T @ unit[state])
        state_rows.append(
            {
                "state": state,
                "label": STATE_LABELS[state],
                "embedding_dimension": raw[state].size,
                "raw_norm": np.linalg.norm(raw[state]),
                "plane_x": coords[state][0],
                "plane_y": coords[state][1],
                "common_plane_offset_norm": np.linalg.norm(perpendicular),
            }
        )

    pair_rows = []
    for index, left in enumerate(STATES):
        for right in STATES[index + 1 :]:
            similarity = float(np.dot(unit[left], unit[right]))
            pair_rows.append(
                {
                    "state_a": left,
                    "state_b": right,
                    "cosine_similarity": similarity,
                    "cosine_distance": 1.0 - similarity,
                    "unit_chord_distance": float(np.linalg.norm(unit[left] - unit[right])),
                }
            )

    angle_rows = []
    for vertex in STATES:
        others = [state for state in STATES if state != vertex]
        a = unit[others[0]] - unit[vertex]
        b = unit[others[1]] - unit[vertex]
        cosine = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
        angle_rows.append(
            {
                "vertex_state": vertex,
                "toward_state_a": others[0],
                "toward_state_b": others[1],
                "interior_angle_degrees": float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))),
            }
        )
    return pd.DataFrame(state_rows), pd.DataFrame(pair_rows), pd.DataFrame(angle_rows)


def load_overexpression_shifts(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    columns = ["Gene_name", "source_state", "goal_state", "overexpress_shift", "gene_set"]
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"Shift table is missing columns: {sorted(missing)}")
    frame = frame[columns].copy()
    frame["source_state"] = frame["source_state"].str.lower()
    frame["goal_state"] = frame["goal_state"].str.lower()
    frame["overexpress_shift"] = pd.to_numeric(frame["overexpress_shift"], errors="raise")
    frame = frame[frame["source_state"].isin(STATES) & frame["goal_state"].isin(STATES)]
    counts = frame.groupby(["Gene_name", "source_state"])["goal_state"].nunique()
    if len(counts) != 50 * 3 or not (counts == 2).all():
        raise ValueError("Expected 50 genes x 3 sources with exactly two goal shifts each")
    return frame


def program_for_gene(gene: str) -> str:
    for program, genes in PROGRAMS.items():
        if gene in genes:
            return program
    return "none"


def reconstruct_in_plane(
    shifts: pd.DataFrame,
    unit: dict[str, np.ndarray],
    basis: np.ndarray,
    coords: dict[str, np.ndarray],
    max_condition: float = 1e8,
) -> pd.DataFrame:
    rows = []
    for (gene, source), group in shifts.groupby(["Gene_name", "source_state"], sort=True):
        goals = [state for state in STATES if state != source]
        values = group.set_index("goal_state")["overexpress_shift"]
        if set(values.index) != set(goals):
            raise ValueError(f"{gene}/{source} does not have shifts for {goals}")
        score_matrix = np.vstack([unit[goal] @ basis for goal in goals])
        condition = float(np.linalg.cond(score_matrix))
        if not np.isfinite(condition) or condition > max_condition:
            raise GeometryError(f"Ill-conditioned reconstruction for {source}: condition={condition:.3g}")
        observed = np.array([values[goal] for goal in goals], dtype=float)
        vector = np.linalg.solve(score_matrix, observed)
        predicted = score_matrix @ vector
        row = {
            "gene": gene,
            "gene_set": group["gene_set"].iloc[0],
            "program": program_for_gene(gene),
            "source_state": source,
            "start_x": coords[source][0],
            "start_y": coords[source][1],
            "shift_x": vector[0],
            "shift_y": vector[1],
            "end_x": coords[source][0] + vector[0],
            "end_y": coords[source][1] + vector[1],
            "vector_norm": float(np.linalg.norm(vector)),
            "score_matrix_condition": condition,
            "max_score_reconstruction_error": float(np.max(np.abs(predicted - observed))),
            "assumption": "unique displacement in affine centroid direction plane; out-of-plane component unidentifiable",
        }
        for state in STATES:
            row[f"similarity_shift_{state}"] = float(values[state]) if state in values else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_programs(vectors: pd.DataFrame) -> pd.DataFrame:
    rows = []
    named = vectors[vectors["program"] != "none"]
    groups = [(program, source, group) for (source, program), group in named.groupby(["source_state", "program"])]
    groups += [("all_50_genes", source, group) for source, group in vectors.groupby("source_state")]
    for program, source, group in groups:
        mean_vector = group[["shift_x", "shift_y"]].mean().to_numpy()
        mean_gene_norm = group["vector_norm"].mean()
        rows.append(
            {
                "program": program,
                "source_state": source,
                "n_genes": len(group),
                "mean_shift_x": mean_vector[0],
                "mean_shift_y": mean_vector[1],
                "mean_vector_norm": float(np.linalg.norm(mean_vector)),
                "mean_gene_vector_norm": float(mean_gene_norm),
                "directional_concentration": float(np.linalg.norm(mean_vector) / mean_gene_norm) if mean_gene_norm else np.nan,
            }
        )
    return pd.DataFrame(rows)


def complement_geometry(
    unit: dict[str, np.ndarray], basis: np.ndarray, observed_path: Path
) -> pd.DataFrame:
    observed = pd.read_csv(observed_path).set_index("source_state") if observed_path.exists() else None
    rows = []
    for source, predictor, response in COMPLEMENT:
        predictor_plane = basis.T @ unit[predictor]
        response_plane = basis.T @ unit[response]
        rows.append(
            {
                "source_state": source,
                "predictor_goal": predictor,
                "response_goal": response,
                "observed_slope": float(observed.loc[source, "slope"]) if observed is not None else np.nan,
                "full_space_isotropic_slope": float(np.dot(unit[response], unit[predictor])),
                "centroid_plane_isotropic_slope": float(
                    np.dot(response_plane, predictor_plane) / np.dot(predictor_plane, predictor_plane)
                ),
                "interpretation": "isotropic reference only; actual slope depends on displacement covariance",
            }
        )
    return pd.DataFrame(rows)


def plot_geometry(
    coords: dict[str, np.ndarray], vectors: pd.DataFrame, programs: pd.DataFrame, output: Path
) -> float:
    import matplotlib.pyplot as plt

    exhaustion = vectors[vectors["program"] == "exhaustion"]
    means = programs[programs["program"] == "exhaustion"].set_index("source_state")
    triangle_span = max(np.linalg.norm(coords[a] - coords[b]) for a in STATES for b in STATES)
    max_mean = means["mean_vector_norm"].max()
    scale = 0.22 * triangle_span / max_mean if max_mean > 0 else 1.0

    fig, ax = plt.subplots(figsize=(7.2, 6.2), dpi=300)
    polygon = np.vstack([coords[state] for state in STATES] + [coords[STATES[0]]])
    ax.plot(polygon[:, 0], polygon[:, 1], color="#8793a1", linewidth=1.2, zorder=1)
    for state in STATES:
        start = coords[state]
        for row in exhaustion[exhaustion["source_state"] == state].itertuples():
            ax.arrow(
                start[0], start[1], row.shift_x * scale, row.shift_y * scale,
                width=0, head_width=0, color=STATE_COLORS[state], alpha=0.18,
                length_includes_head=True, zorder=2,
            )
        mean = means.loc[state]
        ax.arrow(
            start[0], start[1], mean.mean_shift_x * scale, mean.mean_shift_y * scale,
            width=0.006 * triangle_span, head_width=0.045 * triangle_span,
            color=STATE_COLORS[state], length_includes_head=True, zorder=4,
        )
        ax.scatter(*start, s=130, color=STATE_COLORS[state], edgecolor="white", linewidth=1.2, zorder=5)
        ax.annotate(STATE_LABELS[state], start, xytext=(7, 7), textcoords="offset points", weight="bold")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("centroid-plane coordinate 1 (Normal → LUAD)")
    ax.set_ylabel("centroid-plane coordinate 2 (SCLC-positive)")
    ax.set_title("T3 — exhaustion-program overexpression in the centroid plane", loc="left", weight="bold")
    ax.text(
        0.01, 0.01,
        f"Thin arrows: 7 genes; bold: mean; arrow scale = {scale:.2g}×\n"
        "Vectors assume displacement lies in the centroid direction plane.",
        transform=ax.transAxes, fontsize=8, color="#4a5768", va="bottom",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return float(scale)


def run_analysis(centroid_path: Path, shift_path: Path, results_dir: Path, figures_dir: Path) -> dict:
    raw = load_centroids(centroid_path)
    unit = normalized(raw)
    _, basis, coords = centroid_plane(unit)
    states, pairs, angles = geometry_tables(raw, unit, basis, coords)
    shifts = load_overexpression_shifts(shift_path)
    vectors = reconstruct_in_plane(shifts, unit, basis, coords)
    programs = summarize_programs(vectors)
    complement = complement_geometry(unit, basis, results_dir / "axis_goal_complement_regression.csv")

    results_dir.mkdir(parents=True, exist_ok=True)
    states.to_csv(results_dir / "t3_centroid_coordinates.csv", index=False)
    pairs.to_csv(results_dir / "t3_centroid_pairwise_geometry.csv", index=False)
    angles.to_csv(results_dir / "t3_centroid_interior_angles.csv", index=False)
    vectors.to_csv(results_dir / "t3_overexpression_vectors.csv", index=False)
    programs.to_csv(results_dir / "t3_program_vectors.csv", index=False)
    complement.to_csv(results_dir / "t3_complement_geometry.csv", index=False)
    arrow_scale = plot_geometry(coords, vectors, programs, figures_dir / "t3_embedding_geometry.png")

    sclc_angle = float(angles.set_index("vertex_state").loc["sclc", "interior_angle_degrees"])
    base = float(np.linalg.norm(coords["luad"] - coords["normal"]))
    height = float(abs(coords["sclc"][1]))
    offsets = states["common_plane_offset_norm"].to_numpy()
    summary = {
        "status": "complete",
        "embedding_dimension": int(next(iter(raw.values())).size),
        "sclc_interior_angle_degrees": sclc_angle,
        "sclc_height_over_normal_luad_base": height / base,
        "triangle_area_unit_centroid_space": 0.5 * base * height,
        "common_plane_offset_norm_range": [float(offsets.min()), float(offsets.max())],
        "n_genes": int(vectors["gene"].nunique()),
        "n_sources": int(vectors["source_state"].nunique()),
        "arrow_display_scale": arrow_scale,
        "reconstruction_assumption": "displacement lies in affine centroid direction plane",
        "identifiability_warning": (
            "Two retained goal shifts do not identify out-of-plane displacement. "
            "Complement slope is not fixed by centroids alone; isotropic reference slopes are reported."
        ),
    }
    (results_dir / "t3_geometry_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centroids", type=Path, default=DEFAULT_CENTROIDS)
    parser.add_argument("--shifts", type=Path, default=DEFAULT_SHIFTS)
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--figures-dir", type=Path, default=FIGURES)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(json.dumps(run_analysis(args.centroids, args.shifts, args.results_dir, args.figures_dir), indent=2))


if __name__ == "__main__":
    main()
