"""Matched-control table construction for the S100 ISP design
(hive/reports/s100-isp-design-20260922.md, "Controls and analysis", item 3).

STATUS: BLOCKED. build_matched_control_table() is a stub that raises. Do
not implement it against a guessed definition -- an analysis that silently
runs on a guessed control table produces a number that looks exactly like
the real one, and there is no downstream check that can tell the
difference. Two definitional questions are open with the human/Pam
(card s100-isp-execution-20260922, gate 3(b)):

1. Which file and column define "median token rank" for a gene. The
   existing ambient-risk table (sclc_validation/primary_test_perturbation/
   tables/ambient_risk_all_genes.csv) has a detection fraction (detect_frac)
   for ~11,000 genes but no token-rank column of any kind. Detection
   fraction must not be substituted for it.
2. Whether the 0.5-log2-detection-unit / 5-median-rank-percentile-point
   matching tolerance is computed per source-state (SCLC vs LUAD have
   different expression profiles) or globally across the whole panel
   (consistent with the six strata being fixed once, not per source).

Everything downstream of this module (ambient_stats.py's E/Q/stratum-
rollup/primary-test/stability code) is fully built and tested against
synthetic_control_table() below, which is NOT a guess at the real answer --
it is a structurally-valid fake control assignment (fabricated gene ids
that cannot collide with anything real) used only to prove the downstream
math is correct before the real controls exist. It must never be used for
a real analysis.
"""
from __future__ import annotations

MIN_COMMON_CONTROLS = 20

# The six fixed strata, exactly as registered in the design doc and as
# present in the "stratum" column of s100_gene_panel_20260922.json.
STRATA = {
    "clean_high": ("S100A4", "S100A6", "S100A10", "S100A11"),
    "low_p_a16": ("S100P", "S100A16"),
    "low_a2_b": ("S100A2", "S100B"),
    "singleton_a13": ("S100A13",),
    "singleton_a8": ("S100A8",),
    "low_a9_pbp": ("S100A9", "S100PBP"),
}

_BLOCKED_MESSAGE = (
    "Matched-control construction is blocked (s100-isp-execution-20260922, "
    "gate 3(b)) on two open definitional questions with the human/Pam, "
    "not a missing implementation: "
    "(1) which file/column defines 'median token rank' for a gene -- "
    "ambient_risk_all_genes.csv has detect_frac but no rank column, and "
    "detection fraction must not be substituted for it; "
    "(2) whether the 0.5-log2-detection / 5-median-rank-percentile matching "
    "tolerance is computed per source-state or globally across the panel. "
    "See hive/reports/s100-isp-design-20260922.md, 'Controls and analysis', "
    "item 3, and matched_controls.py's module docstring. Do not resolve "
    "this by guessing a stand-in -- a guessed control table produces a "
    "number indistinguishable from the real one."
)


def build_matched_control_table(
    strata: dict[str, tuple[str, ...]],
    ambient_risk_table_path,
    token_rank_source_path,
    *,
    per_source_state: bool,
    detection_tolerance_log2: float = 0.5,
    rank_tolerance_percentile: float = 5.0,
    min_common_controls: int = MIN_COMMON_CONTROLS,
) -> dict[str, tuple[str, ...]]:
    """Build the six fixed 20-control matched strata tables.

    Intended return shape once unblocked: {stratum_name: (ensembl_id, ...)}
    with exactly `min_common_controls` non-S100, non-ambient-anchor genes
    per stratum, each within `detection_tolerance_log2` log2 detection units
    and `rank_tolerance_percentile` median-rank-percentile points of EVERY
    member of that stratum (not a stratum centroid), sampled once with seed
    20260922 and never re-sampled. A stratum with fewer than
    `min_common_controls` common eligible controls is not_estimable and
    must be reported as such, never widened, split, or given private
    controls.

    `token_rank_source_path` has no confirmed value yet -- that is exactly
    open question (1) above. Its presence in this signature documents what
    the interface needs, not what currently exists.

    Raises NotImplementedError unconditionally. See module docstring.
    """
    raise NotImplementedError(_BLOCKED_MESSAGE)


def synthetic_control_table(
    strata: dict[str, tuple[str, ...]] = STRATA,
    n_controls: int = MIN_COMMON_CONTROLS,
) -> dict[str, tuple[str, ...]]:
    """A structurally-valid FAKE control table for testing ambient_stats.py
    only. Fabricated ensembl-shaped ids ("ENSG9SYNTH...") that cannot
    collide with any real gene, so a caller can never mistake this for real
    controls by accident. Deterministic (no randomness needed -- the ids
    are just enumerated), since determinism of the *test* isn't the thing
    seed 20260922 is protecting.

    NEVER use this for a real analysis. It exists only so the E/Q/
    stratum-rollup/primary-test/stability code can be built and verified
    before build_matched_control_table() is unblocked.
    """
    table: dict[str, tuple[str, ...]] = {}
    for stratum_name in strata:
        table[stratum_name] = tuple(
            f"ENSG9SYNTH{stratum_name.upper()}{i:03d}" for i in range(n_controls)
        )
    return table
