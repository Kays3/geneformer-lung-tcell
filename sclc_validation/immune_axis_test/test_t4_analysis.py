#!/usr/bin/env python3
"""CPU-only tests for T4 raw-shift summaries and decision tables."""
from __future__ import annotations

import pickle
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_t4_overexpression as analysis


class T4AnalysisTests(unittest.TestCase):
    def test_load_raw_shifts(self) -> None:
        payload = {}
        for state in analysis.runner.STATE_NAMES.values():
            payload[state] = defaultdict(list, {((1, 2), "cell_emb"): [0.1, -0.2]})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.pickle"
            with path.open("wb") as handle:
                pickle.dump(payload, handle)
            shifts = analysis.load_raw_shifts(path)
        self.assertEqual(set(shifts), set(analysis.runner.STATE_NAMES))
        np.testing.assert_allclose(shifts["sclc"], [0.1, -0.2])

    def test_run_summary_reports_donor_and_subtype(self) -> None:
        run = {"id": "program_exhaustion", "phase": "program", "program": "exhaustion"}
        metadata = pd.DataFrame(
            {
                "individual": ["D1", "D1", "D2", "D2"],
                "celltype": ["x"] * 4,
                "length": [4, 3, 2, 1],
                "cd4cd8": ["CD4", "CD4", "CD8", "CD8"],
            }
        )
        shifts = {state: np.array([0.1, 0.3, -0.1, 0.1]) for state in analysis.runner.STATE_NAMES}
        summary, donors = analysis.summarize_run(run, "sclc", shifts, metadata)
        overall = next(row for row in summary if row["target_state"] == "luad" and row["scope"] == "all")
        self.assertAlmostEqual(overall["mean_shift"], 0.1)
        self.assertEqual(overall["n_donors"], 2)
        self.assertEqual({row["scope"] for row in summary}, {"all", "CD4", "CD8"})
        self.assertEqual(len(donors), 12)

    def test_matched_null_empirical_p(self) -> None:
        rows = [
            {
                "item_id": "program_exhaustion",
                "phase": "program",
                "program": "exhaustion",
                "source_state": "sclc",
                "target_state": "luad",
                "scope": "all",
                "mean_shift": 1.0,
            }
        ]
        rows += [
            {
                "item_id": f"null_{index}",
                "phase": "null",
                "program": "exhaustion",
                "source_state": "sclc",
                "target_state": "luad",
                "scope": "all",
                "mean_shift": value,
            }
            for index, value in enumerate(np.linspace(-0.2, 0.2, 20))
        ]
        result = analysis.matched_null_tests(pd.DataFrame(rows), expected_null=20).iloc[0]
        self.assertTrue(result["null_complete"])
        self.assertAlmostEqual(result["empirical_p_directional"], 1 / 21)

    def test_nested_monotonicity_identifies_primary_path(self) -> None:
        frame = pd.DataFrame(
            [
                {
                    "phase": "nested", "scope": "all", "source_state": "sclc",
                    "target_state": "luad", "nested_size": size, "mean_shift": size * 0.01,
                }
                for size in range(1, 8)
            ]
        )
        result = analysis.nested_monotonicity(frame).iloc[0]
        self.assertTrue(result["fully_non_decreasing"])
        self.assertTrue(result["pre_registered_primary"])
        self.assertAlmostEqual(result["spearman_size_vs_shift"], 1.0)


if __name__ == "__main__":
    unittest.main()
