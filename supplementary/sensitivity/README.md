# GSE68310 sensitivity analysis

`run_gse68310_sensitivity.R` evaluates GSE68310 by subgroup and across enrollment seasons, using the data and main nested-CV results produced by `run_all.sh`.

## Analysis

**Subgroup performance.** The script reuses the main analysis's outer held-out predictions and calculates AUROC, sensitivity, specificity and balanced accuracy separately for the 2009 and 2010 enrollment groups, monoinfections and coinfections. No models are refitted. Sensitivity and specificity use a fixed threshold of 0.5; AUROC intervals use 2,000 outcome-stratified bootstrap resamples of the saved predictions.

**Cross-season validation.** Enrollment year defines two outer test sets: train on 2009 and test on 2010, then reverse the years. For each training year, `make_tensorpls_folds()` generates five stratified inner folds using training-year subjects only. These are not the folds from the main 5 × 5 nested CV. Within the training year, TensorPLS tunes H = 1–10 using inner Q²Y, with Tucker3 imputation (ranks 4, 4, 3) and fixed `VIP2D_intersect_M2` selection at the 99th percentile. The selected pipeline is refitted on the full training year and evaluated on the held-out year. Test-year Q²Y uses the corresponding training-year mean outcome as its baseline.

## Run

After `run_all.sh` has finished, set the paths to the paper repository, the TensorPLS source, and the main run output. For example, if the main run was created with `REPRO_RUN_DIR="$HOME/outputs/main_analysis"`:

```bash
P="$HOME/TensorPLS_reproduction/tensorPLSPaper"
SOURCE="$HOME/TensorPLS_reproduction/TensorPLS"
RUN="$HOME/outputs/main_analysis"

TENSORPLS_SOURCE_DIR="$SOURCE" \
GSE68310_DATA="$RUN/GSE68310_data" \
GSE68310_RESULTS="$RUN/GSE68310_analysis" \
GSE68310_SENSITIVITY_OUTPUT="$RUN/GSE68310_sensitivity" \
  Rscript --vanilla "$P/supplementary/sensitivity/run_gse68310_sensitivity.R"
```

These paths do not depend on the current working directory. Set `RUN` to the actual output directory if it differs from this example.

The sensitivity output directory must not exist before the run. The optional `GSE68310_SEASON_FOLDS` environment variable checks regenerated fold assignments against a saved RDS; it does not replace the generated folds.

## Results

The output includes `subgroup_metrics.csv` (subgroup performance), `cross_season_metrics.csv` (both transfer directions), `cross_season_predictions.csv`, `cross_season_folds.rds`, `cross_season_result.rds` and `integrity_checks.csv`. The script also saves the run settings and R session information.

The expression data were already processed before deposition in GEO. The train-only TensorPLS workflow does not establish that upstream preprocessing was independent across subjects.
