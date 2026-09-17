#!/usr/bin/env python3
"""Regression tests for T6 donor-balanced robustness.

unittest, not pytest (pytest is not installed on this machine; `python3 -m pytest`
returning "No module named pytest" is a missing tool, not a failing test).

Two of these encode mistakes that were actually made and caught while writing T6:
  * `LodoSignMetricTests` - the LODO sign-stability metric was first written comparing
    each drop against the FIRST drop rather than against the full-data value. That is
    silently wrong whenever the first drop is itself the sign-flipping one.
  * `ScaleSeparationTests` - two committed tables covering the same cells differ by
    ~0.79x, so anything that pools or differences across them is invalid.
"""
from __future__ import annotations

import unittest
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

import donor_robustness as m

RESULTS = Path(__file__).resolve().parent / "results"


class PermutationFloorTests(unittest.TestCase):
    """The floor must be computed from the design, before any p-value is read."""

    def test_balanced_split_matches_two_over_central_binomial(self) -> None:
        for k in (2, 3, 4, 5):
            floor = m.permutation_floor(k, k)
            self.assertEqual(floor["n_assignments"], comb(2 * k, k))
            self.assertAlmostEqual(floor["min_two_sided_p"], 2 / comb(2 * k, k))

    def test_three_versus_four_floor_exceeds_alpha_05(self) -> None:
        """The test-only design: 2/35 = 0.0571, so a two-sided 0.05 gate cannot pass."""
        floor = m.permutation_floor(3, 4)
        self.assertEqual(floor["n_assignments"], 35)
        self.assertAlmostEqual(floor["min_two_sided_p"], 2 / 35)
        self.assertGreater(floor["min_two_sided_p"], 0.05)
        self.assertFalse(floor["two_sided_alpha_05_attainable"])

    def test_nineteen_versus_twentytwo_floor_permits_a_real_test(self) -> None:
        floor = m.permutation_floor(19, 22)
        self.assertEqual(floor["n_assignments"], comb(41, 19))
        self.assertTrue(floor["two_sided_alpha_05_attainable"])
        self.assertLess(floor["min_two_sided_p"], 1e-10)

    def test_floor_is_symmetric_in_its_arguments(self) -> None:
        """The floor does not depend on which group is named first. The record of
        which group was which does, so compare the floor, not the whole dict."""
        a, b = m.permutation_floor(3, 4), m.permutation_floor(4, 3)
        self.assertEqual(a["n_assignments"], b["n_assignments"])
        self.assertAlmostEqual(a["min_two_sided_p"], b["min_two_sided_p"])
        self.assertEqual(a["two_sided_alpha_05_attainable"], b["two_sided_alpha_05_attainable"])


class WeightingTests(unittest.TestCase):
    """Cell-weighting vs donor-level, the mechanism T6 is about."""

    def _frame(self, scores, cells, states, population="p"):
        return pd.DataFrame({
            "population": population, "state": states,
            "donor": [f"d{i}" for i in range(len(scores))],
            "score": scores, "n_cells": cells,
        })

    def test_one_dominant_donor_reverses_the_ordering(self) -> None:
        """A synthetic replica of the real failure: group A's mean is higher per donor,
        but a single low donor carrying most cells reverses the cell-weighted ordering."""
        frame = pd.concat([
            self._frame([0.9, 0.9, 0.1], [10, 10, 1000], ["a"] * 3),
            self._frame([0.5, 0.5, 0.5], [100, 100, 100], ["b"] * 3),
        ])
        table = m.weighting_table(frame)
        a = table[table.state == "a"].iloc[0]
        b = table[table.state == "b"].iloc[0]
        self.assertGreater(a.donor_level_mean, b.donor_level_mean)      # donor-level: a > b
        self.assertLess(a.cell_weighted_mean, b.cell_weighted_mean)     # cell-weighted: reversed
        self.assertGreater(a.max_donor_cell_share, 0.9)

    def test_equal_cell_counts_make_the_two_weightings_identical(self) -> None:
        frame = self._frame([0.2, 0.4, 0.9], [50, 50, 50], ["a"] * 3)
        table = m.weighting_table(frame).iloc[0]
        self.assertAlmostEqual(table.cell_weighted_mean, table.donor_level_mean)


class PermutationTestTests(unittest.TestCase):
    def test_exact_enumeration_on_a_small_separable_case(self) -> None:
        """Perfect separation at 3v3 cannot beat the floor: p must equal 2/C(6,3)=0.10."""
        frame = pd.DataFrame({
            "population": "p", "state": ["a"] * 3 + ["b"] * 3,
            "donor": list("uvwxyz"), "score": [1.0, 1.1, 1.2, 5.0, 5.1, 5.2],
            "n_cells": [10] * 6,
        })
        result = m.permutation_test(frame, "a", "b", "t", np.random.default_rng(0))
        self.assertEqual(result["method"], "exact_enumeration")
        self.assertAlmostEqual(result["p_two_sided"], 2 / 20)
        self.assertAlmostEqual(result["p_two_sided"], result["min_attainable_two_sided_p"])

    def test_identical_groups_give_p_of_one(self) -> None:
        frame = pd.DataFrame({
            "population": "p", "state": ["a"] * 3 + ["b"] * 3,
            "donor": list("uvwxyz"), "score": [1.0] * 6, "n_cells": [10] * 6,
        })
        result = m.permutation_test(frame, "a", "b", "t", np.random.default_rng(0))
        self.assertAlmostEqual(result["p_two_sided"], 1.0)

    def test_result_always_carries_its_floor_and_its_null(self) -> None:
        frame = pd.DataFrame({
            "population": "p", "state": ["a"] * 3 + ["b"] * 4,
            "donor": list("uvwxyzq"), "score": [1.0, 2.0, 3.0, 1.5, 2.5, 3.5, 4.0],
            "n_cells": [10] * 7,
        })
        result = m.permutation_test(frame, "a", "b", "t", np.random.default_rng(0))
        self.assertIn("null", result)
        self.assertIn("exchangeable", result["null"])
        self.assertGreaterEqual(result["p_two_sided"], result["min_attainable_two_sided_p"])


class LodoSignMetricTests(unittest.TestCase):
    """Sign stability must be judged against the FULL-data value, never against the
    first drop. Written because the first version of this metric did the latter."""

    def test_sign_stability_is_measured_against_full_data(self) -> None:
        # Full data: a - b is positive. Dropping a's high donor makes it negative.
        frame = pd.DataFrame({
            "population": "p", "state": ["a", "a", "b", "b"],
            "donor": ["a1", "a2", "b1", "b2"], "score": [10.0, 1.0, 2.0, 2.0],
            "n_cells": [10, 10, 10, 10],
        })
        lodo = m.leave_one_donor_out(frame, "a", "b", "p")
        dropped_high = lodo[lodo.dropped_donor == "a1"].iloc[0]
        self.assertLess(dropped_high.donor_level_difference, 0)
        self.assertFalse(bool(dropped_high.donor_level_sign_matches_full))
        self.assertTrue(bool(lodo[lodo.dropped_donor == "a2"].iloc[0].donor_level_sign_matches_full))

    def test_every_donor_is_dropped_exactly_once(self) -> None:
        frame = pd.DataFrame({
            "population": "p", "state": ["a"] * 3 + ["b"] * 4,
            "donor": list("uvwxyzq"), "score": np.linspace(1, 2, 7), "n_cells": [10] * 7,
        })
        lodo = m.leave_one_donor_out(frame, "a", "b", "p")
        self.assertEqual(len(lodo), 7)
        self.assertEqual(sorted(lodo.dropped_donor), sorted("uvwxyzq"))
        self.assertTrue((lodo.n_donors_remaining == 6).all())


class RealDataTests(unittest.TestCase):
    """Against the committed tables. Skip rather than fail if a table is absent."""

    def test_baseline_table_reproduces_t2_program_summary(self) -> None:
        """T2's published exhaustion figures must come back out of its own input."""
        summary_path = RESULTS / "t2_program_summary.csv"
        if not summary_path.exists():
            self.skipTest("t2_program_summary.csv absent")
        donors = m.test_only_donor_scores()
        summary = pd.read_csv(summary_path)
        summary = summary[summary.program == "exhaustion"].set_index("state")
        table = m.weighting_table(donors).set_index("state")
        for state in ("luad", "sclc", "normal"):
            self.assertAlmostEqual(table.loc[state, "cell_weighted_mean"],
                                   summary.loc[state, "mean_log1p_cp10k"], places=6,
                                   msg=f"{state} does not reproduce T2")

    def test_donor_level_reverses_the_test_only_ordering(self) -> None:
        """The headline reconciliation, asserted on the real table."""
        table = m.weighting_table(m.test_only_donor_scores()).set_index("state")
        self.assertGreater(table.loc["luad", "cell_weighted_mean"],
                           table.loc["sclc", "cell_weighted_mean"])
        self.assertGreater(table.loc["sclc", "donor_level_mean"],
                           table.loc["luad", "donor_level_mean"])

    def test_both_populations_agree_at_donor_level(self) -> None:
        if not (RESULTS / "t5a_donor_composition_check.csv").exists():
            self.skipTest("t5a_donor_composition_check.csv absent")
        donors = pd.concat([m.test_only_donor_scores(), m.complete_donor_scores()])
        table = m.weighting_table(donors).set_index(["population", "state"])
        for population in ("test_only", "complete"):
            self.assertGreater(table.loc[(population, "sclc"), "donor_level_mean"],
                               table.loc[(population, "luad"), "donor_level_mean"],
                               msg=f"{population}: donor-level SCLC should exceed LUAD")

    def test_one_sclc_test_donor_holds_most_of_the_cell_mass(self) -> None:
        table = m.weighting_table(m.test_only_donor_scores()).set_index("state")
        self.assertGreater(table.loc["sclc", "max_donor_cell_share"], 0.7)

    def test_normal_test_split_has_no_donor_replication(self) -> None:
        """Guards the report's statement that no Normal contrast is testable."""
        table = m.weighting_table(m.test_only_donor_scores()).set_index("state")
        self.assertEqual(int(table.loc["normal", "n_donors"]), 1)


class ScaleConsistencyTests(unittest.TestCase):
    """Matched T2/T5 inputs must use the same canonical CP10k scale."""

    def test_t2_and_t5_pipelines_agree_on_the_same_test_cells(self) -> None:
        other = RESULTS / "pseudobulk_per_donor_test_only.csv"
        if not other.exists():
            self.skipTest("pseudobulk_per_donor_test_only.csv absent")
        t2 = m.weighting_table(m.test_only_donor_scores()).set_index("state")
        frame = pd.read_csv(other)
        frame = frame[frame.gene_symbol.isin(m.EXHAUSTION)]
        donor = frame.groupby(["state", "donor"]).agg(
            score=("mean_log1p_cp10k", "mean"), n_cells=("n_cells", "first")).reset_index()
        for state in ("normal", "sclc", "luad"):
            sub = donor[donor.state == state]
            t5_value = np.average(sub.score, weights=sub.n_cells)
            self.assertAlmostEqual(
                t5_value, t2.loc[state, "cell_weighted_mean"], places=6,
                msg=f"{state}: matched T2/T5 cells and genes must share the CP10k scale",
            )

    def test_manifest_records_the_normalization_invariant(self) -> None:
        path = RESULTS / "t6_manifest.json"
        if not path.exists():
            self.skipTest("t6_manifest.json absent; run donor_robustness.py")
        import json
        manifest = json.loads(path.read_text())
        self.assertIn("must agree", manifest["normalization_note"])
        self.assertEqual(manifest["replication_unit"], "donor")


class Cd4Cd8Tests(unittest.TestCase):
    def test_all_four_strata_present_and_direction_recorded(self) -> None:
        _, stratified = m.cd4cd8_tables()
        summary = m.stratified_summary(stratified, "sclc", "luad")
        self.assertEqual(len(summary), 4)
        self.assertEqual(set(summary.cd4cd8), {"CD4", "CD4 (Treg)", "CD8", "other"})
        self.assertTrue(summary.direction_matches_unstratified.all(),
                        "report states SCLC >= LUAD in every stratum")

    def test_composition_shares_sum_to_one_per_donor(self) -> None:
        composition, _ = m.cd4cd8_tables()
        totals = composition.groupby(["state", "donor"])["share_of_donor_cells"].sum()
        np.testing.assert_allclose(totals.to_numpy(), 1.0, atol=1e-9)

    def test_pleural_effusion_is_a_composition_outlier(self) -> None:
        """The report calls it a compositional outlier, not just a low scorer."""
        composition, _ = m.cd4cd8_tables()
        sclc = composition[composition.state == "sclc"]
        cd8 = sclc[sclc.cd4cd8 == "CD8"].set_index("donor")["share_of_donor_cells"]
        self.assertLess(cd8["PleuralEffusion"], 0.2)
        self.assertTrue((cd8.drop("PleuralEffusion") > 0.4).all())


if __name__ == "__main__":
    unittest.main()
