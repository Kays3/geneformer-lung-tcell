#!/usr/bin/env python3
"""CPU-only verification for run_targeted_panel.py's paired_eligible_dataset()
(card isp-runner-paired-arms-20260922). Not a synthetic-fixture unit test --
it runs directly against the real held-out dataset and the real 50-gene
panel on the compute host, since the thing being checked (token-positive
counts, donor coverage, the seeded per-donor cap) only means something
against real data. No GPU, no model load, no torch forward pass.

DO NOT run this as-is for any panel other than the default 50-gene
target_gene_panel.json it was written against (2026-09-22, flagged by
Stanley during s100-isp-execution-20260922's gate review, before he ran
it): it hardcodes that file below rather than honoring
TARGET_GENES_FILE_OVERRIDE, so running it unmodified for e.g. the S100
override panel silently verifies the WRONG genes. It also actively
mutates shared state -- recompute_matches() deletes and rebuilds cache
entries under PAIRED_ELIGIBLE_DIR/MANIFEST_DIR as part of verifying
determinism -- which is fine when you are its author verifying your own
just-written code, but is a real hazard for anyone else running it
against a shared checkout expecting a read-only check. Point it at a
scratch/override panel and expect it to touch the cache before running
it against anything other than the default panel.

Checks, for every (gene, source) in target_gene_panel.json:
  - the saved dataset is sorted by "length" descending (the property that
    makes InSilicoPerturber's own internal re-sort a no-op -- see
    run_targeted_panel.py's module docstring);
  - the manifest CSV's (cell_id, donor) rows line up positionally,
    row-for-row, with the saved dataset;
  - the per-donor cap (100 cells, seed 20260922) is applied exactly,
    donor-by-donor;
  - a second call is a pure cache hit (identical eligibility info);
  - recomputing from scratch (cache cleared) reproduces the identical
    manifest, both for a donor under the cap and for one over it --
    the actual claim under test, since a cap only exercises the seeded
    subsample path when some donor in the real panel actually exceeds it.

Run from sclc_validation/perturbation_workflow/targeted_panel/:
    python3 verify_paired_eligible_dataset.py
"""
from __future__ import annotations

import json
import shutil
import sys

sys.argv = ["run_targeted_panel.py", "--run-tag", "isp_paired_verify"]
sys.path.insert(0, ".")
import run_targeted_panel as rtp  # noqa: E402


def recompute_matches(source: str, symbol: str, ensembl: str) -> bool:
    """Clear the cache for one (gene, source) and confirm rebuilding it from
    scratch reproduces the identical manifest -- the real test of seed
    determinism, since a cache hit alone would prove nothing."""
    key = f"{source}_{ensembl}"
    manifest_before = rtp.pd.read_csv(rtp.MANIFEST_DIR / f"{key}.csv")
    shutil.rmtree(rtp.PAIRED_ELIGIBLE_DIR / f"{key}.dataset", ignore_errors=True)
    (rtp.PAIRED_ELIGIBLE_DIR / f"{key}.eligibility.json").unlink()
    rtp.paired_eligible_dataset(source, {"gene": symbol, "ensembl_id": ensembl})
    manifest_after = rtp.pd.read_csv(rtp.MANIFEST_DIR / f"{key}.csv")
    return manifest_before.equals(manifest_after)


def main() -> None:
    genes = json.load(open("target_gene_panel.json"))["genes"]
    cap_exercised_example = None

    for gene in genes:
        symbol, ensembl = gene["gene"], gene["ensembl_id"]
        for source in ("sclc", "luad", "normal"):
            path, info = rtp.paired_eligible_dataset(source, {"gene": symbol, "ensembl_id": ensembl})
            print(f"{symbol}/{source}: eligible={info['eligible']} n_cells={info.get('n_cells')} "
                  f"n_donors={info.get('n_donors')} reason={info.get('reason')}")
            if not info["eligible"]:
                continue

            ds = rtp.load_from_disk(str(path))
            assert len(ds) == info["n_cells"]
            lengths = ds["length"]
            assert lengths == sorted(lengths, reverse=True), f"{symbol}/{source}: not length-sorted"

            manifest = rtp.pd.read_csv(info["manifest"])
            assert len(manifest) == len(ds)
            assert list(manifest["donor"]) == ds["individual"], f"{symbol}/{source}: donor order mismatch"
            assert list(manifest["cell_id"]) == ds["cell_id"], f"{symbol}/{source}: cell_id order mismatch"

            for donor, count in info["donor_cell_counts"].items():
                expected = min(count, rtp.DONOR_CELL_CAP)
                actual = int(manifest["donor"].eq(donor).sum())
                assert actual == expected, f"{symbol}/{source}/{donor}: expected {expected}, got {actual}"
                if count > rtp.DONOR_CELL_CAP and cap_exercised_example is None:
                    cap_exercised_example = (source, symbol, ensembl)
            print("  -> sorted, manifest/dataset aligned, donor caps correct")

    print(f"\nDonor cap (>{rtp.DONOR_CELL_CAP} cells) exercised by this panel: {cap_exercised_example is not None}")

    _, info_a = rtp.paired_eligible_dataset("luad", {"gene": "S100A2", "ensembl_id": "ENSG00000196754"})
    _, info_b = rtp.paired_eligible_dataset("luad", {"gene": "S100A2", "ensembl_id": "ENSG00000196754"})
    assert info_a == info_b, "cache hit returned different eligibility info"
    print("Cache hit: identical eligibility info on repeat call -- OK")

    assert recompute_matches("luad", "S100A2", "ENSG00000196754")
    print("Recompute-from-scratch (no cap needed): identical manifest -- OK")

    if cap_exercised_example:
        source, symbol, ensembl = cap_exercised_example
        assert recompute_matches(source, symbol, ensembl)
        print(f"Recompute-from-scratch ({symbol}/{source}, cap WAS exercised): "
              f"identical manifest -- seeded per-donor subsample is deterministic")
    else:
        print("WARNING: no (gene, source, donor) in this panel exceeded the cap -- "
              "the seeded-subsample branch was not exercised by this run of this script.")

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
