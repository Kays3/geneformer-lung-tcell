#!/usr/bin/env bash
# E2 analysis (registration s.6), CPU only, after run_e2_isp.sh reports completed.
# 1 build_design_e2.py        design.json from committed pre-GPU inputs + the committed classifier gate
# 2 build_ovx_index.py x2     overexpress token-positive positions and count checks (panel+controls; null genes)
# 3 analyse.py (unchanged)    H2c per-gene statuses, lung rules file, host drift false (one host)
# 4 null_analysis.py (unch.)  H2b primary, LOO, bootstrap, no-op gate, validity check on the E2 controls
# 5 compare_panel_b.py        H2c primary reading against the LUAD rows
# Usage: run_e2_analysis.sh <repo_worktree> <data_root> <geneformer_root> <venv_python>
set -euo pipefail
REPO="$1"; DATA="$2"; GF="$3"; PY="$4"
E2="$REPO/pelka_crc_e2"; BD="$REPO/balanced_donor_luad"; OUT="$E2/results"; ISP="$DATA/isp/out"
TOK="$GF/geneformer/token_dictionary_gc104M.pkl"; MAN="$DATA/goals/donor_manifest.json"
mkdir -p "$OUT"
[ -f "$E2/phase4_results/classifier_gate.json" ] || { echo "classifier gate result not committed"; exit 3; }
[ "$(cat "$ISP/stop_reason.txt")" = "completed" ] || { echo "ISP did not complete; a partial arm is never analysed"; exit 3; }
sha256sum -c "$E2/registration/pins_isp.sha256"
"$PY" "$E2/scripts/build_design_e2.py" --base "$E2" --out "$OUT/design.json"
"$PY" "$BD/scripts/build_ovx_index.py" --phase6-root "$ISP" --design "$OUT/design.json" --manifest "$MAN" \
  --token-dict "$TOK" --out "$DATA/isp/ovx_index.json" | tee "$OUT/ovx_index_checks.json"
"$PY" "$BD/scripts/build_ovx_index.py" --phase6-root "$ISP" --design "$OUT/design.json" --manifest "$MAN" \
  --token-dict "$TOK" --genes-file "$E2/null/null_genes_100_estimable.json" --out "$DATA/isp/null_ovx_index.json" \
  | tee "$OUT/null_ovx_index_checks.json"
sha256sum "$DATA/isp/ovx_index.json" "$DATA/isp/null_ovx_index.json" > "$OUT/ovx_index_sha256.txt"
"$PY" "$BD/scripts/analyse.py" --phase6-root "$ISP" --design "$OUT/design.json" \
  --rules "$BD/registration/phase7_rules.json" --out "$OUT/panel_b" --host-drift false \
  --ovx-index "$DATA/isp/ovx_index.json"
set +e
"$PY" "$BD/scripts/null_analysis.py" --phase8-root "$ISP" --null-genes "$E2/null/null_genes_100_estimable.json" \
  --frozen "$E2/null/frozen_100_estimable.json" --noop-results "$DATA/isp/noop_spotchecks/results.jsonl" \
  --design "$OUT/design.json" --ovx-index "$DATA/isp/null_ovx_index.json" \
  --phase6-controls-root "$ISP" --controls-ovx-index "$DATA/isp/ovx_index.json" --out "$OUT/h2b_null_result.json"
echo "null_analysis exit $?" | tee "$OUT/h2b_exit.txt"
set -e
"$PY" "$E2/scripts/compare_panel_b.py" --luad-rows "$BD/phase7_results/outcome_rows.json" \
  --e2-rows "$OUT/panel_b/outcome_rows.json" --e2-per-donor "$OUT/panel_b/per_donor_adjusted.json" \
  --cohort-definition "$E2/cohort/COHORT_DEFINITION.json" --out "$OUT/h2c_result.json"
cp "$DATA/isp/noop_spotchecks/results.jsonl" "$OUT/noop_spotchecks.jsonl"
echo E2_ANALYSIS_DONE
