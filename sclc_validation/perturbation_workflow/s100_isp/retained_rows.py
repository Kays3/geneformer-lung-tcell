"""Build s100_isp_retained_rows.csv (retained_rows_spec_20260922.md).

The whole point of this contract is that a missing stats row must remain a
row: build_retained_rows() always returns exactly the planned Cartesian row
count (12 genes x 2 sources x 2 goals x 2 operations = 96 for the core
panel), asserted at the end, regardless of how much of the real run has
completed. A row that can't be computed yet gets a status and a
status_reason explaining why -- it is never dropped, and a not-estimable
or not-yet-run row's numeric columns are null, never zero.

Deliberately has NO dependency on geneformer/torch: the caller passes in
`gene_token_dict` (loaded however the real runner loads it) so this module
stays testable with plain Python + pandas + numpy against synthetic
fixtures, without needing a GPU host or the vendored Geneformer checkout.

KNOWN GAP, surfaced rather than resolved quietly (2026-09-23,
s100-isp-execution-20260922): retained_rows_spec_20260922.md's `status`
enum (eligible_completed / not_estimable_cell_count / not_estimable_
donor_count / not_estimable_control_stratum / run_failed / no_op_failed)
has no value for "eligible, but the run has not happened yet". The spec's
own framing ("Build ... before any result filtering") reads as though it
is meant to run only after a real run completes, when every planned row is
guaranteed to land in one of those six buckets. This module adds a
STATUS_NOT_RUN = "not_run" value outside that enum to describe the
pre-run and partially-run states this table is built in *before* Module A
executes (i.e., right now). Flagged to Michael/Pam rather than assumed;
if the intent really is "this table only exists post-run", STATUS_NOT_RUN
rows simply won't occur in that case and this is dead code, not a wrong
answer.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from provenance_utils import sha256_file

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
STATES = (SCLC, LUAD, NORMAL)
SLUGS = {SCLC: "sclc", LUAD: "luad", NORMAL: "normal"}
STATE_BY_SLUG = {v: k for k, v in SLUGS.items()}

# Module A only ever uses SCLC/LUAD as sources -- Normal has a single donor
# in the held-out split, so the >=3-donor eligibility gate can never pass
# for any gene, by construction (confirmed live, s100-isp-execution-20260922
# gate 3(a)). Including it here would only ever produce not_estimable_donor_
# count rows, which is correct but adds nothing -- the design doc's own
# Module A section already excludes it.
MODULE_A_SOURCES = ("sclc", "luad")
PERTURBATION_TYPES = ("delete", "overexpress")

STATUS_ELIGIBLE_COMPLETED = "eligible_completed"
STATUS_NOT_ESTIMABLE_CELL_COUNT = "not_estimable_cell_count"
STATUS_NOT_ESTIMABLE_DONOR_COUNT = "not_estimable_donor_count"
STATUS_NOT_ESTIMABLE_CONTROL_STRATUM = "not_estimable_control_stratum"
STATUS_RUN_FAILED = "run_failed"
STATUS_NO_OP_FAILED = "no_op_failed"
# Not in retained_rows_spec_20260922.md's enum -- see module docstring.
STATUS_NOT_RUN = "not_run"


def planned_cartesian_rows(
    panel_genes: list[dict],
    sources: tuple[str, ...] = MODULE_A_SOURCES,
    perturbation_types: tuple[str, ...] = PERTURBATION_TYPES,
) -> pd.DataFrame:
    """The fixed Cartesian design: gene x source x goal x perturbation, two
    goals per source (the other two disease states). Generic over
    `panel_genes` so the identical function builds the later control-gene
    rows at the same grain once matched_controls.py is unblocked --
    retained_rows_spec_20260922.md: "Add control rows at the same grain
    after the six fixed 20-control strata have been locked."
    """
    rows = []
    for gene in panel_genes:
        for source_slug in sources:
            source_state = STATE_BY_SLUG[source_slug]
            goals = [s for s in STATES if s != source_state]
            for goal_state in goals:
                alt_state = next(s for s in STATES if s not in (source_state, goal_state))
                for ptype in perturbation_types:
                    rows.append({
                        "gene": gene["gene"],
                        "ensembl_id": gene["ensembl_id"],
                        "role": gene.get("role"),
                        "stratum": gene.get("stratum"),
                        "source_state": source_state,
                        "goal_state": goal_state,
                        "alt_state": alt_state,
                        "perturbation_type": ptype,
                    })
    return pd.DataFrame(rows)


def load_eligibility(paired_eligible_dir, source_slug: str, ensembl_id: str) -> dict | None:
    info_path = Path(paired_eligible_dir) / f"{source_slug}_{ensembl_id}.eligibility.json"
    if not info_path.exists():
        return None
    return json.loads(info_path.read_text())


def load_completion_marker(raw_root, perturb_type: str, source_slug: str, symbol: str) -> dict | None:
    marker_path = Path(raw_root) / perturb_type / source_slug / f"targeted_{source_slug}_{symbol}.complete.json"
    if not marker_path.exists():
        return None
    return json.loads(marker_path.read_text())


def completion_marker_path(raw_root, perturb_type: str, source_slug: str, symbol: str) -> Path:
    return Path(raw_root) / perturb_type / source_slug / f"targeted_{source_slug}_{symbol}.complete.json"


def find_raw_pickle(raw_root, perturb_type: str, source_slug: str, symbol: str):
    """Exactly the same glob pattern run_gene()/run_noop_gene() use for
    their own cleanup/output-counting. NPROC=1 in the current runner
    configuration guarantees at most one matching file per (gene, source,
    perturb_type) in practice -- but this is asserted, not assumed: more
    than one match is reported as an ambiguous run_failed rather than
    silently taking matches[0], since concatenation order across multiple
    worker-pickle files has not been verified against the real Geneformer
    source and a wrong guess there would silently corrupt donor alignment.
    """
    raw_dir = Path(raw_root) / perturb_type / source_slug
    prefix = f"targeted_{source_slug}_{symbol}"
    matches = sorted(raw_dir.glob(f"in_silico_{perturb_type}_{prefix}_*_raw.pickle"))
    if len(matches) > 1:
        return "ambiguous", matches
    return (matches[0] if matches else None), matches


def donor_balanced_shift(
    raw_pickle_path,
    goal_state_name: str,
    ensembl_id: str,
    gene_token_dict: dict,
    manifest_df: pd.DataFrame,
) -> dict:
    """The design's actual required statistic (hive/reports/
    s100-isp-design-20260922.md: "compute a mean shift within each donor
    first, then give each donor equal weight") -- deliberately NOT
    InSilicoPerturberStats' own Shift_to_goal_end column, which is a plain
    mean over all cells (cell-weighted, not donor-weighted). Reads the raw
    per-cell pickle directly and aligns it against the paired-eligible
    manifest's donor column, which is written in the exact same cell order
    by construction (run_targeted_panel.py's paired_eligible_dataset(),
    PR #23) -- so this alignment is sound by construction, not
    reconstruction, the same property that replaced the old
    donor_consistency.py script.
    """
    import pickle

    token = gene_token_dict.get(ensembl_id)
    if token is None:
        return {"donor_balanced_shift": None, "donor_sign_fraction": None,
                "n_donors_scored": None, "status": STATUS_RUN_FAILED,
                "reason": f"{ensembl_id} not present in the token dictionary"}
    with open(raw_pickle_path, "rb") as f:
        cos_sims_dict = pickle.load(f)
    inner = cos_sims_dict.get(goal_state_name)
    if not inner:
        return {"donor_balanced_shift": None, "donor_sign_fraction": None,
                "n_donors_scored": None, "status": STATUS_RUN_FAILED,
                "reason": f"goal state {goal_state_name!r} missing from raw pickle {Path(raw_pickle_path).name}"}
    key = (token, "cell_emb")
    if key not in inner:
        return {"donor_balanced_shift": None, "donor_sign_fraction": None,
                "n_donors_scored": None, "status": STATUS_RUN_FAILED,
                "reason": f"key {key!r} missing from raw pickle {Path(raw_pickle_path).name}"}
    values = inner[key]
    if len(values) != len(manifest_df):
        return {"donor_balanced_shift": None, "donor_sign_fraction": None,
                "n_donors_scored": None, "status": STATUS_RUN_FAILED,
                "reason": f"length mismatch: {len(values)} shift values vs {len(manifest_df)} manifest rows"}
    df = manifest_df.copy()
    df["shift"] = values
    donor_means = df.groupby("donor")["shift"].mean()
    donor_balanced = float(donor_means.mean())
    sign = np.sign(donor_balanced)
    same_sign = (np.sign(donor_means) == sign) | (donor_means == 0)
    return {
        "donor_balanced_shift": donor_balanced,
        "donor_sign_fraction": float(same_sign.mean()),
        "n_donors_scored": int(len(donor_means)),
        "status": STATUS_ELIGIBLE_COMPLETED,
        "reason": None,
    }


def build_retained_rows(
    panel_id: str,
    panel_genes: list[dict],
    *,
    raw_root,
    paired_eligible_dir,
    stats_root,
    gene_token_dict: dict,
    run_id: str,
    runner_sha256: str,
    panel_sha256: str,
    sources: tuple[str, ...] = MODULE_A_SOURCES,
    perturbation_types: tuple[str, ...] = PERTURBATION_TYPES,
    min_cells_eligible: int = 50,
    min_donors_eligible: int = 3,
) -> pd.DataFrame:
    planned = planned_cartesian_rows(panel_genes, sources, perturbation_types)
    manifest_cache: dict[str, pd.DataFrame] = {}
    stats_cache: dict[str, pd.DataFrame] = {}
    out_rows = []

    for _, planned_row in planned.iterrows():
        row = planned_row.to_dict()
        source_slug = SLUGS[row["source_state"]]
        goal_slug = SLUGS[row["goal_state"]]
        symbol, ensembl_id, ptype = row["gene"], row["ensembl_id"], row["perturbation_type"]

        row["panel_id"] = panel_id
        row["run_id"] = run_id
        row["runner_sha256"] = runner_sha256
        row["panel_sha256"] = panel_sha256
        row.update({
            "status": None, "status_reason": None,
            "n_token_positive_cells": None, "n_eligible_donors": None, "donor_cell_counts": None,
            "paired_cell_manifest_sha256": None, "paired_cell_count": None,
            "no_op_status": None, "no_op_score": None,
            "raw_completion_marker": str(completion_marker_path(raw_root, ptype, source_slug, symbol)),
            "raw_stats_path": None, "raw_stats_sha256": None, "stats_row_present": False,
            "donor_balanced_shift": None, "donor_sign_fraction": None, "matched_control_percentile_q": None,
        })

        elig = load_eligibility(paired_eligible_dir, source_slug, ensembl_id)
        if elig is None:
            row["status"] = STATUS_NOT_RUN
            row["status_reason"] = "no eligibility record -- paired_eligible_dataset() has not been run yet for this (gene, source)"
            out_rows.append(row)
            continue

        row["n_token_positive_cells"] = elig.get("n_cells")
        row["n_eligible_donors"] = elig.get("n_donors")
        if elig.get("donor_cell_counts") is not None:
            row["donor_cell_counts"] = json.dumps(elig["donor_cell_counts"], sort_keys=True)

        if not elig.get("eligible"):
            n_donors = elig.get("n_donors") or 0
            row["status"] = (
                STATUS_NOT_ESTIMABLE_DONOR_COUNT if n_donors < min_donors_eligible
                else STATUS_NOT_ESTIMABLE_CELL_COUNT
            )
            row["status_reason"] = elig.get("reason")
            out_rows.append(row)
            continue

        manifest_path = elig.get("manifest")
        if manifest_path not in manifest_cache:
            manifest_cache[manifest_path] = pd.read_csv(manifest_path)
        manifest_df = manifest_cache[manifest_path]
        row["paired_cell_manifest_sha256"] = sha256_file(manifest_path)
        row["paired_cell_count"] = int(len(manifest_df))

        marker = load_completion_marker(raw_root, ptype, source_slug, symbol)
        if marker is None:
            row["status"] = STATUS_NOT_RUN
            row["status_reason"] = "eligible but no completion marker -- this (gene, source, operation) has not been run yet"
            out_rows.append(row)
            continue

        if marker.get("skipped_zero_cells_detected"):
            row["status"] = STATUS_RUN_FAILED
            row["status_reason"] = "library reported zero cells detected despite passing the paired-eligibility gate -- unexpected; flagged rather than accepted silently"
            out_rows.append(row)
            continue

        raw_pickle, matches = find_raw_pickle(raw_root, ptype, source_slug, symbol)
        if raw_pickle is None:
            row["status"] = STATUS_RUN_FAILED
            row["status_reason"] = "completion marker present but no matching raw pickle found on disk"
            out_rows.append(row)
            continue
        if raw_pickle == "ambiguous":
            row["status"] = STATUS_RUN_FAILED
            row["status_reason"] = f"{len(matches)} raw pickle files matched -- refusing to guess concatenation order (NPROC=1 should make this impossible)"
            out_rows.append(row)
            continue

        shift = donor_balanced_shift(raw_pickle, row["goal_state"], ensembl_id, gene_token_dict, manifest_df)
        row["donor_balanced_shift"] = shift["donor_balanced_shift"]
        row["donor_sign_fraction"] = shift["donor_sign_fraction"]
        if shift["status"] != STATUS_ELIGIBLE_COMPLETED:
            row["status"] = shift["status"]
            row["status_reason"] = shift["reason"]
            out_rows.append(row)
            continue

        row["status"] = STATUS_ELIGIBLE_COMPLETED
        row["status_reason"] = None

        stats_path = Path(stats_root) / ptype / f"targeted_{ptype}_{source_slug}_to_{goal_slug}.csv"
        row["raw_stats_path"] = str(stats_path)
        if stats_path.exists():
            if str(stats_path) not in stats_cache:
                stats_cache[str(stats_path)] = pd.read_csv(stats_path)
            stats_df = stats_cache[str(stats_path)]
            row["raw_stats_sha256"] = sha256_file(stats_path)
            row["stats_row_present"] = bool((stats_df["Ensembl_ID"] == ensembl_id).any())

        noop_marker = load_completion_marker(raw_root, "noop", source_slug, symbol)
        if noop_marker is None:
            row["no_op_status"] = STATUS_NOT_RUN
        elif noop_marker.get("skipped_zero_cells_detected"):
            row["no_op_status"] = STATUS_NO_OP_FAILED
        else:
            noop_pickle, noop_matches = find_raw_pickle(raw_root, "noop", source_slug, symbol)
            if noop_pickle in (None, "ambiguous"):
                row["no_op_status"] = STATUS_NO_OP_FAILED
            else:
                noop_shift = donor_balanced_shift(noop_pickle, row["goal_state"], ensembl_id, gene_token_dict, manifest_df)
                if noop_shift["status"] == STATUS_ELIGIBLE_COMPLETED:
                    row["no_op_status"] = STATUS_ELIGIBLE_COMPLETED
                    row["no_op_score"] = noop_shift["donor_balanced_shift"]
                else:
                    row["no_op_status"] = STATUS_NO_OP_FAILED

        out_rows.append(row)

    result = pd.DataFrame(out_rows)
    expected_n = len(panel_genes) * len(sources) * 2 * len(perturbation_types)
    assert len(result) == expected_n, (
        f"retained-row count drifted from the planned Cartesian design: "
        f"got {len(result)}, expected {expected_n} -- a row was dropped, which "
        f"this table must never do"
    )
    return result


def status_counts(retained_rows: pd.DataFrame) -> pd.Series:
    """Per retained_rows_spec_20260922.md: 'The final report must give
    counts for every status.'"""
    return retained_rows["status"].value_counts(dropna=False).sort_index()
