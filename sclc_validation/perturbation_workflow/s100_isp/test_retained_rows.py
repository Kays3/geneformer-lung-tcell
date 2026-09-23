#!/usr/bin/env python3
"""Synthetic-fixture test for retained_rows.py (s100-isp-execution-20260922,
"build the analysis layer" card). No GPU, no real dataset, no geneformer
import -- builds a small fake raw_root/paired_eligible_dir/stats_root tree
covering every status branch, then checks build_retained_rows() against it.

Run: python3 test_retained_rows.py
"""
from __future__ import annotations

import json
import pickle
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, ".")
import retained_rows as rr

import os
SCRATCH = Path(os.environ.get("TMPDIR", "/tmp")) / "kevin_s100_analysis_layer_test"


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _write_manifest(path: Path, cell_ids, donors, lengths) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": cell_ids, "donor": donors, "length": lengths}).to_csv(path, index=False)


def _write_pickle(path: Path, cos_sims_dict: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(cos_sims_dict, f)


def build_fixture() -> dict:
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    raw_root = SCRATCH / "raw"
    paired_eligible_dir = SCRATCH / "paired_eligible"
    manifest_dir = paired_eligible_dir / "manifests"
    stats_root = SCRATCH / "stats"

    gene_token_dict = {"ENSGA": 111, "ENSGB": 222}
    panel_genes = [
        {"gene": "GENEA", "ensembl_id": "ENSGA", "role": "clean_control", "stratum": "s1"},
        {"gene": "GENEB", "ensembl_id": "ENSGB", "role": "flagged", "stratum": "s2"},
    ]

    # --- GENEA / sclc: fully eligible + completed for both delete and overexpress ---
    manifest_a_sclc = manifest_dir / "sclc_ENSGA.csv"
    _write_manifest(manifest_a_sclc, [f"c{i}" for i in range(6)],
                     ["D1", "D1", "D1", "D1", "D2", "D3"], [500, 480, 460, 440, 420, 400])
    _write_json(paired_eligible_dir / "sclc_ENSGA.eligibility.json", {
        "eligible": True, "n_cells": 6, "n_donors": 3,
        "donor_cell_counts": {"D1": 4, "D2": 1, "D3": 1},
        "seed": 20260922, "donor_cell_cap": 100, "manifest": str(manifest_a_sclc),
    })
    # Donor cell counts (4/1/1) deliberately UNEQUAL and mirror the real
    # LUAD S100A2 imbalance (40/4/4/25) this is actually protecting -- with
    # equal per-donor counts, a cell-weighted mean and a donor-balanced
    # mean are ARITHMETICALLY IDENTICAL and the assertion below would pass
    # just as happily against the cell-weighted bug this test exists to
    # catch (caught live by Michael, 2026-09-23: the previous 2/2/2 fixture
    # gave both computations 0.3, so it asserted a fact without ever being
    # able to observe its negation -- "a test whose pass carries no
    # information is not protection, it is a claim of protection"). With
    # this fixture: donor-weighted = (0.9+0.1-0.1)/3 = 0.3 (unchanged
    # expected value); cell-weighted = (0.9*4+0.1-0.1)/6 = 0.6 (would now
    # FAIL the assertion below if the code ever regressed to it).
    shifts_luad = [0.9, 0.9, 0.9, 0.9, 0.1, -0.1]  # donor means: D1=.9 D2=.1 D3=-.1 -> balanced mean = .3
    shifts_normal = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    for ptype in ("delete", "overexpress"):
        raw_dir = raw_root / ptype / "sclc"
        _write_json(raw_dir / "targeted_sclc_GENEA.complete.json", {
            "completed_utc": "2026-09-23T00:00:00Z", "perturb_type": ptype, "source": "sclc",
            "gene": "GENEA", "ensembl_id": "ENSGA", "elapsed_seconds": 1.0, "n_raw_files": 1,
        })
        _write_pickle(raw_dir / f"in_silico_{ptype}_targeted_sclc_GENEA_batch0_raw.pickle", {
            "lung adenocarcinoma": {(111, "cell_emb"): shifts_luad},
            "normal": {(111, "cell_emb"): shifts_normal},
        })
    # noop for GENEA/sclc
    raw_noop_dir = raw_root / "noop" / "sclc"
    _write_json(raw_noop_dir / "targeted_sclc_GENEA.complete.json", {
        "completed_utc": "2026-09-23T00:00:00Z", "perturb_type": "noop", "source": "sclc",
        "gene": "GENEA", "ensembl_id": "ENSGA", "elapsed_seconds": 1.0, "n_cells": 6, "n_raw_files": 1,
    })
    _write_pickle(raw_noop_dir / "in_silico_noop_targeted_sclc_GENEA_batch0_raw.pickle", {
        "lung adenocarcinoma": {(111, "cell_emb"): [0.001, -0.001, 0.0005, -0.0005, 0.0002, -0.0002]},
        "normal": {(111, "cell_emb"): [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]},
    })
    # ISP stats CSVs (delete/overexpress x goal luad/normal) -- provenance only
    for ptype in ("delete", "overexpress"):
        for goal_slug, goal_name in (("luad", "lung adenocarcinoma"), ("normal", "normal")):
            stats_path = stats_root / ptype / f"targeted_{ptype}_sclc_to_{goal_slug}.csv"
            stats_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{"Gene_name": "GENEA", "Ensembl_ID": "ENSGA", "Shift_to_goal_end": 0.3}]).to_csv(stats_path, index=False)

    # --- GENEA / luad: not estimable (donor count) ---
    _write_json(paired_eligible_dir / "luad_ENSGA.eligibility.json", {
        "eligible": False, "n_cells": 40, "n_donors": 2,
        "reason": "40 token-positive cells from 2 donors (need >= 50 cells from >= 3 donors)",
    })

    # --- GENEB / sclc: eligible but not yet run (no completion marker) ---
    manifest_b_sclc = manifest_dir / "sclc_ENSGB.csv"
    _write_manifest(manifest_b_sclc, [f"c{i}" for i in range(4)], ["D1", "D1", "D2", "D2"], [300, 290, 280, 270])
    _write_json(paired_eligible_dir / "sclc_ENSGB.eligibility.json", {
        "eligible": True, "n_cells": 4, "n_donors": 2,  # NOTE: 2 donors would actually fail the real gate;
        # kept anyway to exercise the "eligible but not run" branch independent of the gate check itself.
        "donor_cell_counts": {"D1": 2, "D2": 2}, "seed": 20260922, "donor_cell_cap": 100,
        "manifest": str(manifest_b_sclc),
    })

    # --- GENEB / luad: eligible, completed, but the library reported zero cells (run_failed) ---
    manifest_b_luad = manifest_dir / "luad_ENSGB.csv"
    _write_manifest(manifest_b_luad, [f"c{i}" for i in range(4)], ["D1", "D2", "D3", "D4"], [300, 290, 280, 270])
    _write_json(paired_eligible_dir / "luad_ENSGB.eligibility.json", {
        "eligible": True, "n_cells": 4, "n_donors": 4,
        "donor_cell_counts": {"D1": 1, "D2": 1, "D3": 1, "D4": 1}, "seed": 20260922,
        "donor_cell_cap": 100, "manifest": str(manifest_b_luad),
    })
    for ptype in ("delete", "overexpress"):
        _write_json(raw_root / ptype / "luad" / "targeted_luad_GENEB.complete.json", {
            "completed_utc": "2026-09-23T00:00:00Z", "perturb_type": ptype, "source": "luad",
            "gene": "GENEB", "ensembl_id": "ENSGB", "elapsed_seconds": 1.0, "n_raw_files": 0,
            "skipped_zero_cells_detected": True,
        })

    return {
        "raw_root": raw_root, "paired_eligible_dir": paired_eligible_dir,
        "stats_root": stats_root, "gene_token_dict": gene_token_dict, "panel_genes": panel_genes,
    }


def main() -> None:
    fx = build_fixture()
    # Explicit sources=("sclc", "luad") -- deliberately not relying on
    # MODULE_A_SOURCES' default. This test exercises build_retained_rows()'s
    # generic per-source Cartesian/status logic, which is unchanged by the
    # 2026-09-23 ruling that dropped SCLC as a *production* source; the
    # fixture below still needs both source branches to cover every status
    # path (see the module docstring's "genuinely tested" standard).
    result = rr.build_retained_rows(
        panel_id="test-panel", panel_genes=fx["panel_genes"],
        raw_root=fx["raw_root"], paired_eligible_dir=fx["paired_eligible_dir"],
        stats_root=fx["stats_root"], gene_token_dict=fx["gene_token_dict"],
        run_id="test-run", runner_sha256="deadbeef", panel_sha256="cafef00d",
        sources=("sclc", "luad"),
    )

    expected_n = 2 * 2 * 2 * 2  # 2 genes x 2 sources x 2 goals x 2 ops
    assert len(result) == expected_n, f"expected {expected_n} rows, got {len(result)}"
    print(f"Row count: {len(result)} (== {expected_n}) -- OK")

    # every planned combination present exactly once
    dupes = result.duplicated(subset=["gene", "source_state", "goal_state", "perturbation_type"])
    assert not dupes.any(), "duplicate planned rows"
    print("No duplicate planned rows -- OK")

    def row(gene, source, goal, ptype):
        m = result[(result.gene == gene) & (result.source_state == source) &
                    (result.goal_state == goal) & (result.perturbation_type == ptype)]
        assert len(m) == 1
        return m.iloc[0]

    # GENEA/sclc/delete/->luad: eligible_completed, donor-balanced (not cell-weighted) shift
    r = row("GENEA", rr.SCLC, rr.LUAD, "delete")
    assert r["status"] == rr.STATUS_ELIGIBLE_COMPLETED, r["status"]
    assert abs(r["donor_balanced_shift"] - 0.3) < 1e-9, r["donor_balanced_shift"]
    # Donor means: D1=.9, D2=.1, D3=-.1 -> balanced mean = .3 (positive).
    # 2 of 3 donors (D1, D2) share that sign; D3 does not -- proves this is
    # a real per-donor computation, not a cell-weighted mean read back.
    assert abs(r["donor_sign_fraction"] - (2 / 3)) < 1e-9, r["donor_sign_fraction"]
    print("GENEA/sclc/delete/->luad: eligible_completed, donor_balanced_shift=0.3, donor_sign_fraction=2/3 -- OK")

    r2 = row("GENEA", rr.SCLC, rr.NORMAL, "delete")
    assert r2["status"] == rr.STATUS_ELIGIBLE_COMPLETED
    assert abs(r2["donor_balanced_shift"] - 0.05) < 1e-9
    assert r2["paired_cell_manifest_sha256"] is not None and len(r2["paired_cell_manifest_sha256"]) == 64
    assert r2["paired_cell_count"] == 6
    assert r2["stats_row_present"] is True or bool(r2["stats_row_present"]) is True
    print("GENEA/sclc/delete/->normal: manifest hash + paired_cell_count + stats join -- OK")

    r3 = row("GENEA", rr.SCLC, rr.LUAD, "delete")
    assert r3["no_op_status"] == rr.STATUS_ELIGIBLE_COMPLETED
    assert r3["no_op_score"] is not None and abs(r3["no_op_score"]) < 0.01
    print("GENEA/sclc no-op join: eligible_completed, small independent-pass score -- OK")

    # GENEA/luad/*: not_estimable_donor_count (2 donors < 3)
    r4 = row("GENEA", rr.LUAD, rr.SCLC, "delete")
    assert r4["status"] == rr.STATUS_NOT_ESTIMABLE_DONOR_COUNT, r4["status"]
    assert r4["donor_balanced_shift"] is None or pd.isna(r4["donor_balanced_shift"])
    print("GENEA/luad: not_estimable_donor_count, donor_balanced_shift null -- OK")

    # GENEB/sclc/*: eligible but no completion marker -> not_run
    r5 = row("GENEB", rr.SCLC, rr.LUAD, "delete")
    assert r5["status"] == rr.STATUS_NOT_RUN, r5["status"]
    print("GENEB/sclc: not_run (eligible, no marker yet) -- OK")

    # GENEB/luad/*: skipped_zero_cells_detected -> run_failed
    r6 = row("GENEB", rr.LUAD, rr.SCLC, "overexpress")
    assert r6["status"] == rr.STATUS_RUN_FAILED, r6["status"]
    assert "zero cells" in r6["status_reason"]
    print("GENEB/luad: run_failed (skipped_zero_cells_detected surfaced, not silently accepted) -- OK")

    counts = rr.status_counts(result)
    print("\nStatus counts:")
    print(counts.to_string())
    assert counts.sum() == expected_n

    print("\nALL RETAINED-ROWS CHECKS PASSED")


if __name__ == "__main__":
    main()
