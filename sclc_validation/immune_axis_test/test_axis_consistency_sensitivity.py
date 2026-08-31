#!/usr/bin/env python3
"""Tests for the T1 ambient/stress/ribosomal sensitivity arm."""
from __future__ import annotations

import unittest

import axis_consistency as axis
import axis_consistency_sensitivity as sensitivity


class AxisSensitivityTests(unittest.TestCase):
    def test_flagged_classes_include_expected_genes(self) -> None:
        excluded, _ = sensitivity.flagged_genes()
        self.assertTrue({"HBA1", "HBA2", "HBB", "HSPA1B", "RPL21", "RPS26", "S100A8"} <= excluded)

    def test_clean_scenario_remains_analyzable(self) -> None:
        wide = axis.load_shifts()
        excluded, _ = sensitivity.flagged_genes()
        clean = wide.drop(index=sorted(excluded & set(wide.index)))
        row, tables = sensitivity.scenario_summary("clean", clean)
        self.assertGreater(row["n_genes"], 20)
        self.assertEqual(len(tables["reciprocity"]), 3)
        self.assertEqual(len(tables["complement"]), 3)
        self.assertAlmostEqual(tables["variance"]["variance_explained"].sum(), 1.0)


if __name__ == "__main__":
    unittest.main()
