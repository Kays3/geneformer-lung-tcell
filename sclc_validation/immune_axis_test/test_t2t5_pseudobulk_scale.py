#!/usr/bin/env python3
"""Unit tests for the denominator audit; no H5AD or compute host required."""
from __future__ import annotations

import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

import audit_t2t5_pseudobulk_scale as m


class NormalizationTests(unittest.TestCase):
    def test_equal_denominators_make_t2_and_t5_identical(self) -> None:
        counts = sp.csr_matrix([[2.0, 8.0], [1.0, 9.0]])
        denominators = np.array([10.0, 10.0])
        np.testing.assert_allclose(
            m.log1p_cp10k(counts, denominators).toarray(),
            m.log1p_cp10k(counts, np.asarray(counts.sum(axis=1)).ravel()).toarray(),
        )

    def test_stored_counts_larger_than_current_matrix_lower_t2_scores(self) -> None:
        counts = sp.csr_matrix([[2.0, 8.0]])
        t2 = m.log1p_cp10k(counts, np.array([20.0])).toarray()
        t5 = m.log1p_cp10k(counts, np.array([10.0])).toarray()
        self.assertTrue(np.all(t5 > t2))

    def test_ratio_summary_detects_nonmatching_denominators(self) -> None:
        summary = m.summarize_ratio(np.array([20.0, 10.0]), np.array([10.0, 10.0]))
        self.assertEqual(summary["median"], 0.75)
        self.assertEqual(summary["fraction_exactly_one"], 0.5)

    def test_state_means_are_cell_weighted(self) -> None:
        norm = sp.csr_matrix([[1.0], [3.0], [5.0], [7.0]])
        state = pd.Series(["sclc", "sclc", "luad", "normal"])
        values = m.state_program_means(norm, state)
        self.assertEqual(values["sclc"], 2.0)
        self.assertEqual(values["luad"], 5.0)
        self.assertEqual(values["normal"], 7.0)


if __name__ == "__main__":
    unittest.main()
