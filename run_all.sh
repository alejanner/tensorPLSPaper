#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export TENSORPLS_SOURCE_DIR="${TENSORPLS_SOURCE_DIR:-$here/../TensorPLS}"
[[ -f "$TENSORPLS_SOURCE_DIR/DESCRIPTION" ]] || { echo 'TensorPLS source not found' >&2; exit 1; }
run="${REPRO_RUN_DIR:-${REPRO_BASE:-$here/outputs}/run_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$(dirname "$run")"
mkdir "$run"
export GSE68310_DATA="$run/GSE68310_data"
export GSE182616_DATA="$run/GSE182616_data"
export GSE68310_RESULTS="$run/GSE68310_analysis"
export GSE182616_RESULTS="$run/GSE182616_analysis"
mkdir -p "$GSE68310_DATA/downloads" "$GSE182616_DATA" "$run/logs"

for script in 01_prepare_GSE68310.R 02_analyze_GSE68310.R 03_prepare_GSE182616.R 04_analyze_GSE182616.R; do
  echo " $script"
  Rscript --vanilla "$here/$script" 2>&1 | tee "$run/logs/${script%.R}.log"
done
for dataset in GSE68310 GSE182616; do
  cat "$run/${dataset}_analysis/COMPLETED.txt" "$run/${dataset}_analysis/reference_comparison.csv"
done