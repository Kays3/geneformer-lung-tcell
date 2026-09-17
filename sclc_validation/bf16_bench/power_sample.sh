#!/usr/bin/env bash
# GPU power/energy/temperature sampling wrapper for a bf16-bench arm.
#
# Samples `nvidia-smi --query-gpu=power.draw,temperature.gpu --format=csv -l 1`
# in the background for the duration of the wrapped command, then integrates
# mean-power x wall-time to Wh (same simple trapezoidal-by-mean approach the
# colleague's analysis/13_profile.py uses) and reports peak temperature.
#
# Usage:
#   power_sample.sh <label> <out_dir> -- <command> [args...]
#
# Writes <out_dir>/<label>.power.csv (raw samples) and
# <out_dir>/<label>.power.json (wall_s, mean/max/min power_w, energy_wh,
# peak temp_c, exit_code).
#
# Note (this host, GB10/thinkstation1): nvidia-smi reports
# "Memory-Usage: Not Supported" on this unified-memory GPU, so this script
# does not attempt a GPU memory column -- peak GPU memory instead comes from
# torch.cuda.max_memory_allocated(), recorded by the runner scripts
# themselves into each unit's completion marker (see dtype_cast.py /
# run_targeted_panel.py's peak_gpu_mem_gib field).
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "usage: $0 <label> <out_dir> -- <command> [args...]" >&2
    exit 2
fi
LABEL="$1"; OUT_DIR="$2"; shift 2
[[ "${1:-}" == "--" ]] && shift
if [[ $# -eq 0 ]]; then
    echo "error: no command given after --" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
CSV="$OUT_DIR/${LABEL}.power.csv"
SUMMARY="$OUT_DIR/${LABEL}.power.json"

echo "timestamp,power.draw [W],temperature.gpu" > "$CSV"
nvidia-smi --query-gpu=timestamp,power.draw,temperature.gpu --format=csv,noheader -l 1 >> "$CSV" &
SAMPLER_PID=$!
trap 'kill "$SAMPLER_PID" 2>/dev/null || true' EXIT

START=$(date +%s.%N)
set +e
"$@"
RC=$?
set -e
END=$(date +%s.%N)

kill "$SAMPLER_PID" 2>/dev/null || true
wait "$SAMPLER_PID" 2>/dev/null || true
trap - EXIT

WALL=$(awk -v a="$START" -v b="$END" 'BEGIN{printf "%.3f", b-a}')

python3 - "$CSV" "$SUMMARY" "$LABEL" "$WALL" "$RC" <<'PYEOF'
import csv, json, sys

csv_path, summary_path, label, wall_s, rc = sys.argv[1:6]
wall_s, rc = float(wall_s), int(rc)

powers, temps = [], []
with open(csv_path) as f:
    for row in csv.reader(f):
        if len(row) < 3 or row[0].strip() == "timestamp":
            continue
        try:
            powers.append(float(row[1].strip().split()[0]))
            temps.append(float(row[2].strip()))
        except (ValueError, IndexError):
            continue

mean_power = sum(powers) / len(powers) if powers else None
energy_wh = mean_power * (wall_s / 3600.0) if mean_power is not None else None

summary = {
    "label": label,
    "wall_s": round(wall_s, 2),
    "exit_code": rc,
    "n_samples": len(powers),
    "gpu_power_w": {
        "mean": round(mean_power, 3) if mean_power is not None else None,
        "max": max(powers) if powers else None,
        "min": min(powers) if powers else None,
    },
    "gpu_energy_wh": round(energy_wh, 4) if energy_wh is not None else None,
    "gpu_temp_c_max": max(temps) if temps else None,
    "csv": csv_path,
}
with open(summary_path, "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
PYEOF

exit "$RC"
