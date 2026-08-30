# Results — T3, embedding geometry

**Status: run.** The trusted three-state training-donor centroid artifact was
exported to a portable `.npz` and analyzed with
[`embedding_geometry.py`](embedding_geometry.py). This is a geometry and
identifiability result; it does not authorize changing poster or talk wording
without the review required in [PLAN.md §10](PLAN.md#10-if-the-resolution-holds).

## What ran

The input contains three aggregate vectors (`normal`, `sclc`, `luad`), each of
dimension 768, with no cell- or donor-level records. The source pickle's SHA256
and the export metadata are recorded in
[`results/training_donor_disease_centroids.json`](results/training_donor_disease_centroids.json);
the portable vectors are in
[`results/training_donor_disease_centroids.npz`](results/training_donor_disease_centroids.npz).
The analysis normalizes each centroid before measuring pairwise geometry. It
then reconstructs a unique displacement in the affine centroid-direction
plane from each gene's two retained non-source overexpression shifts.

Outputs:

- [`results/t3_centroid_pairwise_geometry.csv`](results/t3_centroid_pairwise_geometry.csv)
- [`results/t3_centroid_interior_angles.csv`](results/t3_centroid_interior_angles.csv)
- [`results/t3_centroid_coordinates.csv`](results/t3_centroid_coordinates.csv)
- [`results/t3_overexpression_vectors.csv`](results/t3_overexpression_vectors.csv)
- [`results/t3_program_vectors.csv`](results/t3_program_vectors.csv)
- [`results/t3_complement_geometry.csv`](results/t3_complement_geometry.csv)
- [`results/t3_geometry_summary.json`](results/t3_geometry_summary.json)
- [`figures/t3_embedding_geometry.png`](figures/t3_embedding_geometry.png)

## Finding 1 — the three states are not collinear

The normalized centroids form a genuine triangle, not a one-dimensional line:

| Pair or vertex | Result |
|---|---:|
| Normal–SCLC cosine similarity | −0.3712 |
| Normal–LUAD cosine similarity | −0.3937 |
| SCLC–LUAD cosine similarity | −0.2578 |
| Normal interior angle | 56.97° |
| **SCLC interior angle** | **61.94°** |
| LUAD interior angle | 61.09° |
| Triangle area in unit-centroid space | 1.1590 |
| SCLC height / Normal–LUAD base | 0.8316 |

The centroid direction plane has a common offset of 0.3244 from the origin in
the normalized 768-dimensional space. The substantial SCLC height and
non-zero triangle area reject the geometric premise needed to read two shifts
from one source state as a position on a Normal → SCLC → LUAD line. This
confirms T1's qualitative warning, now measured directly from the centroids.

## Finding 2 — in-plane reconstructions are numerically exact, but conditional

For each of 50 genes and each of the three source states, the two retained
goal shifts solve a 2×2 system in the centroid plane. The maximum absolute
reconstruction error across all 150 vectors is below `7e-17`; this is an
algebraic check, not evidence that the true 768-dimensional displacement has no
component outside the plane. Two scalar shifts cannot identify that component.

The curated exhaustion program's mean in-plane vectors are:

| Source | Mean Δcoordinate 1 (Normal → LUAD) | Mean Δcoordinate 2 | Mean vector norm | Directional concentration |
|---|---:|---:|---:|---:|
| Normal | 0.0153 | −0.0100 | 0.01825 | 0.992 |
| SCLC | 0.0096 | −0.0091 | 0.01319 | 0.792 |
| LUAD | 0.0101 | −0.0055 | 0.01153 | 0.909 |

The full program table also reports cytotoxicity, progenitor/memory, and
SCLC-subtype-TF contrasts. Those means are descriptive summaries of the
existing single-gene OE shifts; they are not program-level perturbation results
and have no uncertainty estimate.

## Finding 3 — T1c's complement slope has no unique centroid-only prediction

The observed complement regressions and two isotropic reference calculations
are:

| Source | Observed slope | Full-space isotropic reference | Centroid-plane isotropic reference |
|---|---:|---:|---:|
| SCLC | −0.202 | −0.394 | −0.558 |
| Normal | −0.397 | −0.258 | −0.406 |
| LUAD | −0.563 | −0.371 | −0.533 |

The references differ because the slope depends on displacement covariance, not
centroid geometry alone. They are therefore benchmarks only; neither is a
unique prediction of T1c. The observed slopes remain the measured ISP result,
while the in-plane vectors are explicitly model-based reconstructions.

## Conclusion and limits

T3 resolves the geometric question: the three normalized state centroids are
non-collinear, so the original one-dimensional axis reading is not supported by
the embedding geometry. It does **not** establish that the model's direction is
biologically meaningless, nor does it identify a unique out-of-plane response or
complement slope. The next unresolved experiment is T4's program-level,
nested-set, matched-null overexpression run on the configured GPU host; its
runner is prepared but has not been launched from this laptop.

No poster or talk wording was changed.
