#!/usr/bin/env python3
"""Regression tests for the deterministic T4 design manifest."""
from __future__ import annotations

import unittest

import build_t4_manifest as t4


class T4ManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = t4.build_manifest()

    def test_expected_run_counts(self) -> None:
        self.assertEqual(
            self.manifest["run_counts"],
            {
                "program_definitions": 4,
                "nested_definitions": 7,
                "null_definitions": 80,
                "gpu_runs_all_sources": 273,
            },
        )

    def test_nested_sets_are_strict_prefixes(self) -> None:
        nested = [run for run in self.manifest["runs"] if run["phase"] == "nested"]
        ordered = self.manifest["design"]["exhaustion_nested_order"]
        self.assertEqual(len(nested), len(ordered))
        for run in nested:
            symbols = [gene["gene"] for gene in run["genes"]]
            self.assertEqual(symbols, ordered[: run["nested_size"]])

    def test_null_sets_match_size_and_exclude_targets(self) -> None:
        targets = {gene for genes in t4.PROGRAMS.values() for gene in genes}
        for run in self.manifest["runs"]:
            if run["phase"] != "null":
                continue
            symbols = [gene["gene"] for gene in run["genes"]]
            self.assertEqual(len(symbols), len(t4.PROGRAMS[run["program"]]))
            self.assertEqual(len(symbols), len(set(symbols)))
            self.assertFalse(set(symbols) & targets)

    def test_same_seed_is_fully_reproducible(self) -> None:
        again = t4.build_manifest()
        self.assertEqual(self.manifest, again)


if __name__ == "__main__":
    unittest.main()
