#!/usr/bin/env python3
"""Small dependency-light regression tests for T3 geometry."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import embedding_geometry as eg


class EmbeddingGeometryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.raw = {
            "normal": np.array([1.0, 0.0, 0.2]),
            "sclc": np.array([0.7, 0.7, 0.2]),
            "luad": np.array([0.0, 1.0, 0.2]),
        }
        self.unit = eg.normalized(self.raw)
        _, self.basis, self.coords = eg.centroid_plane(self.unit)

    def test_npz_loader_rejects_missing_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "centroids.npz"
            np.savez(path, normal=self.raw["normal"], sclc=self.raw["sclc"])
            with self.assertRaisesRegex(ValueError, "missing"):
                eg.load_centroids(path)

    def test_triangle_geometry_is_finite(self) -> None:
        states, pairs, angles = eg.geometry_tables(self.raw, self.unit, self.basis, self.coords)
        self.assertEqual(len(states), 3)
        self.assertEqual(len(pairs), 3)
        self.assertEqual(len(angles), 3)
        self.assertTrue(np.isfinite(angles["interior_angle_degrees"]).all())
        self.assertAlmostEqual(float(angles["interior_angle_degrees"].sum()), 180.0, places=8)

    def test_in_plane_reconstruction_reproduces_both_scores(self) -> None:
        rows = []
        expected = {
            "normal": np.array([0.02, -0.01]),
            "sclc": np.array([-0.01, 0.03]),
            "luad": np.array([0.01, 0.02]),
        }
        for source, vector in expected.items():
            for goal in [state for state in eg.STATES if state != source]:
                rows.append(
                    {
                        "Gene_name": "PDCD1",
                        "source_state": source,
                        "goal_state": goal,
                        "overexpress_shift": float(np.dot(self.unit[goal] @ self.basis, vector)),
                        "gene_set": "panel",
                    }
                )
        reconstructed = eg.reconstruct_in_plane(pd.DataFrame(rows), self.unit, self.basis, self.coords)
        for row in reconstructed.itertuples():
            np.testing.assert_allclose([row.shift_x, row.shift_y], expected[row.source_state], atol=1e-12)
            self.assertLess(row.max_score_reconstruction_error, 1e-12)

    def test_collinear_centroids_are_rejected(self) -> None:
        collinear = {
            "normal": np.array([1.0, 0.0]),
            "sclc": np.array([0.0, 1.0]),
            "luad": np.array([-1.0, 2.0]),
        }
        with self.assertRaisesRegex(eg.GeometryError, "collinear"):
            eg.centroid_plane(collinear)


if __name__ == "__main__":
    unittest.main()
