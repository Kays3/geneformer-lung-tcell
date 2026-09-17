#!/usr/bin/env python3
"""Regression tests for T2's current-matrix library-size denominator."""
from __future__ import annotations

import unittest

import numpy as np
import scipy.sparse as sp

import measure_baseline_expression as m


class CurrentLibrarySizeTests(unittest.TestCase):
    def test_uses_the_current_matrix_row_sum(self) -> None:
        matrix = sp.csr_matrix([[2.0, 8.0], [1.0, 9.0]])
        np.testing.assert_allclose(m.current_library_sizes(matrix), [10.0, 10.0])

    def test_rejects_zero_library_cells(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-positive"):
            m.current_library_sizes(sp.csr_matrix([[0.0, 0.0], [1.0, 2.0]]))

    def test_strict_gene_subset_keeps_full_matrix_denominator(self) -> None:
        full = sp.csr_matrix([[2.0, 8.0, 100.0], [1.0, 9.0, 200.0]])
        selected = full[:, :2]
        counts, sizes = m.selected_counts_with_library_sizes(full, selected)
        np.testing.assert_allclose(counts.toarray(), [[2.0, 8.0], [1.0, 9.0]])
        np.testing.assert_allclose(sizes, [110.0, 210.0])
        self.assertFalse(np.allclose(sizes, np.asarray(selected.sum(axis=1)).ravel()))

    def test_rejects_selected_matrix_with_different_cells(self) -> None:
        with self.assertRaisesRegex(ValueError, "same cells"):
            m.selected_counts_with_library_sizes(
                sp.csr_matrix([[1.0], [2.0]]), sp.csr_matrix([[1.0]])
            )


if __name__ == "__main__":
    unittest.main()
