#!/usr/bin/env python3
"""Synthetic-fixture test for matched_controls.py, now that
build_matched_control_table() is unblocked (s100-isp-execution-20260922,
human ruling 2026-09-23). No GPU, no real dataset, no `datasets` import at
module import time (load_and_freeze_luad_gene_stats is the only function
that touches the real corpus, and it is not exercised here).

Run: python3 test_matched_controls.py
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, ".")
import matched_controls as mc


def test_median_token_rank_and_detection() -> None:
    cells = [
        [10, 20, 30],
        [20, 30],
        [10, 30, 40],
        [20],
    ]
    gene_token_dict = {"GA": 10, "GB": 20, "GC": 30, "GD": 40, "GE": 50}

    df0 = mc.median_token_rank_and_detection(cells, gene_token_dict, rank_convention="0_based").set_index("ensembl_id")
    assert abs(df0.loc["GA", "detect_frac"] - 0.5) < 1e-12
    assert abs(df0.loc["GA", "median_token_rank"] - 0.0) < 1e-12
    assert abs(df0.loc["GB", "detect_frac"] - 0.75) < 1e-12
    assert abs(df0.loc["GB", "median_token_rank"] - 0.0) < 1e-12
    assert abs(df0.loc["GC", "detect_frac"] - 0.75) < 1e-12
    assert abs(df0.loc["GC", "median_token_rank"] - 1.0) < 1e-12
    assert abs(df0.loc["GD", "detect_frac"] - 0.25) < 1e-12
    assert abs(df0.loc["GD", "median_token_rank"] - 2.0) < 1e-12
    assert df0.loc["GE", "detect_frac"] == 0.0
    assert np.isnan(df0.loc["GE", "median_token_rank"])
    print("median_token_rank_and_detection (0-based): detect_frac + median rank, incl. never-detected gene -- OK")

    df1 = mc.median_token_rank_and_detection(cells, gene_token_dict, rank_convention="1_based").set_index("ensembl_id")
    assert abs(df1.loc["GA", "median_token_rank"] - 1.0) < 1e-12
    assert abs(df1.loc["GC", "median_token_rank"] - 2.0) < 1e-12
    assert abs(df1.loc["GD", "median_token_rank"] - 3.0) < 1e-12
    # 1-based is exactly 0-based + 1 -- convention shifts every value uniformly,
    # never changes a relative comparison (the ruling's own claim, checked directly).
    for gene in ("GA", "GB", "GC", "GD"):
        assert abs(df1.loc[gene, "median_token_rank"] - (df0.loc[gene, "median_token_rank"] + 1)) < 1e-12
    print("median_token_rank_and_detection (1-based == 0-based + 1, uniformly) -- OK")

    try:
        mc.median_token_rank_and_detection(cells, gene_token_dict, rank_convention="bogus")
        raise AssertionError("expected ValueError for unknown rank_convention")
    except ValueError:
        print("median_token_rank_and_detection rejects unknown rank_convention -- OK")


def test_rank_percentile_transform() -> None:
    values = pd.Series([10.0, 20.0, 20.0, 30.0], index=["a", "b", "c", "d"])
    pct = mc.rank_percentile_transform(values)
    assert abs(pct["a"] - 12.5) < 1e-9
    assert abs(pct["b"] - 50.0) < 1e-9
    assert abs(pct["c"] - 50.0) < 1e-9
    assert abs(pct["d"] - 87.5) < 1e-9
    # cross-check against the exact same formula ambient_stats.midrank_percentile
    # uses (rankdata average, (r - 0.5) / N * 100) computed independently here.
    ranks = rankdata(values.values, method="average")
    expected = 100.0 * (ranks - 0.5) / len(values)
    assert np.allclose(pct.values, expected)
    print("rank_percentile_transform matches independent midrank formula, incl. ties -- OK")


def _build_synthetic_universe():
    """2 panel genes (stratA), 1 anchor gene (stratB, its own stratum), and
    three groups of candidates: 25 that qualify for stratA on both axes, 5
    that fail on detection only, 5 that fail on rank only."""
    panel_genes = [
        {"gene": "PANELX", "ensembl_id": "ENSGPANELX", "stratum": "stratA", "role": "flagged"},
        {"gene": "PANELY", "ensembl_id": "ENSGPANELY", "stratum": "stratA", "role": "intermediate"},
        {"gene": "ANCHORZ", "ensembl_id": "ENSGANCHORZ", "stratum": "stratB", "role": "ambient_anchor"},
    ]
    strata = {"stratA": ("PANELX", "PANELY"), "stratB": ("ANCHORZ",)}

    # median_token_rank is tied at 1000.0 across panel members, "good"
    # candidates, and "bad detection" candidates deliberately -- the
    # rank-percentile tolerance (5 points) is tight relative to a large
    # population's rank granularity (100 / N per step), so this fixture
    # isolates each axis: the tied group differs only on detect_frac, and
    # the "bad rank" group differs only on median_token_rank.
    rows = [
        {"ensembl_id": "ENSGPANELX", "detect_frac": 0.10, "median_token_rank": 1000.0},
        {"ensembl_id": "ENSGPANELY", "detect_frac": 0.10, "median_token_rank": 1000.0},
        # Distinct detect_frac from everything else in this fixture, so
        # zero candidates qualify for stratB -- genuinely not_estimable,
        # not tied into stratA's qualifying pool by fixture accident.
        {"ensembl_id": "ENSGANCHORZ", "detect_frac": 0.50, "median_token_rank": 1000.0},
    ]
    for i in range(25):
        rows.append({"ensembl_id": f"ENSGGOOD{i:03d}", "detect_frac": 0.10, "median_token_rank": 1000.0})
    for i in range(5):
        rows.append({"ensembl_id": f"ENSGBADDET{i:03d}", "detect_frac": 0.90, "median_token_rank": 1000.0})
    for i in range(5):
        rows.append({"ensembl_id": f"ENSGBADRANK{i:03d}", "detect_frac": 0.10, "median_token_rank": 9000.0 + i})
    luad_gene_stats = pd.DataFrame(rows)
    return panel_genes, strata, luad_gene_stats


def test_build_matched_control_table_qualifying_and_deterministic() -> None:
    panel_genes, strata, luad_gene_stats = _build_synthetic_universe()

    result1 = mc.build_matched_control_table(panel_genes, luad_gene_stats, strata=strata)
    result2 = mc.build_matched_control_table(panel_genes, luad_gene_stats, strata=strata)

    a = result1["stratA"]
    assert a["status"] == "eligible", a
    assert a["n_candidates"] == 25, a["n_candidates"]
    assert len(a["controls"]) == mc.MIN_COMMON_CONTROLS == 20
    assert all(cid.startswith("ENSGGOOD") for cid in a["controls"]), a["controls"]
    assert len(set(a["controls"])) == 20, "must not draw duplicates"
    # exclusions: panel genes and anchor genes never appear as controls, even
    # though ANCHORZ's own stats would otherwise be close enough to qualify.
    assert "ENSGPANELX" not in a["controls"] and "ENSGPANELY" not in a["controls"]
    assert "ENSGANCHORZ" not in a["controls"]
    # bad-detection and bad-rank candidates never qualify.
    assert not any(cid.startswith("ENSGBADDET") for cid in a["controls"])
    assert not any(cid.startswith("ENSGBADRANK") for cid in a["controls"])
    print(f"stratA: eligible, 20/{a['n_candidates']} candidates drawn, all from the qualifying pool -- OK")

    assert result1["stratA"]["controls"] == result2["stratA"]["controls"], \
        "same seed must reproduce the exact same draw"
    print("build_matched_control_table is deterministic across repeated calls (seed 20260922) -- OK")

    # stratB has only ANCHORZ as a member and zero non-panel/non-anchor
    # candidates were engineered close to it -- not_estimable, not "rescued".
    b = result1["stratB"]
    assert b["status"] == mc.STATUS_NOT_ESTIMABLE_CONTROL_STRATUM, b
    assert b["controls"] is None
    print(f"stratB: not_estimable_control_stratum ({b['n_candidates']} candidates, < {mc.MIN_COMMON_CONTROLS}) -- OK")


def test_build_matched_control_table_missing_member() -> None:
    """PANELX stays a registered panel gene (it always is -- it's pre-
    registered), but its stats row is absent from luad_gene_stats (e.g. zero
    token-positive LUAD cells, so median_token_rank_and_detection() never
    produced a row for it). This must be reported as not_estimable_control_
    stratum for its stratum, not raise a KeyError."""
    panel_genes, strata, luad_gene_stats = _build_synthetic_universe()
    stats_missing = luad_gene_stats[luad_gene_stats["ensembl_id"] != "ENSGPANELX"]
    assert "ENSGPANELX" not in stats_missing["ensembl_id"].values  # sanity on the fixture itself

    result = mc.build_matched_control_table(panel_genes, stats_missing, strata=strata)
    a = result["stratA"]
    assert a["status"] == mc.STATUS_NOT_ESTIMABLE_CONTROL_STRATUM
    assert a["controls"] is None
    assert "ENSGPANELX" in a["reason"]
    print("build_matched_control_table: stratum member absent from luad_gene_stats -> not_estimable, not a KeyError -- OK")


def test_synthetic_control_table_still_fake() -> None:
    table = mc.synthetic_control_table()
    for stratum, ids in table.items():
        assert len(ids) == mc.MIN_COMMON_CONTROLS
        assert all(cid.startswith("ENSG9SYNTH") for cid in ids)
    print("synthetic_control_table() unaffected by unblocking build_matched_control_table() -- OK")


def main() -> None:
    test_median_token_rank_and_detection()
    test_rank_percentile_transform()
    test_build_matched_control_table_qualifying_and_deterministic()
    test_build_matched_control_table_missing_member()
    test_synthetic_control_table_still_fake()
    print("\nALL MATCHED-CONTROLS CHECKS PASSED")


if __name__ == "__main__":
    main()
