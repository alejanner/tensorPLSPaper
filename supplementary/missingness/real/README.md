# GSE68310 missing-visit reconstruction

This benchmark compares training-set mean imputation and Tucker3 reconstruction when entire visits are missing from held-out subjects. It uses the GSE68310 longitudinal expression data and the five folds saved by the main TensorPLS analysis.

## Inputs

Run `01_prepare_GSE68310.R` and `02_analyze_GSE68310.R` from the main paper repository first. The saved analysis directory must contain:

```text
GSE68310_data/X_incomplete.rds
GSE68310_analysis/nested_cv_folds.rds
```

The tensor contains 109 subjects, 47,162 probes and four visits (Day0, Day2, Day4 and Day6). The original 21 missing visits are kept missing and are never used as evaluation targets.

Install [TensorPLS](https://github.com/alejanner/TensorPLS) and its dependencies before running the scripts. By default, the scripts load the installed package. To use a local source checkout instead, set `TENSORPLS_SOURCE_DIR` and install `pkgload`.

## Reconstruction benchmark

Both imputers are fitted on training subjects only. For each held-out fold, the script masks 5%, 10%, 20%, 30% and 40% of the observed visits, using seeds 123, 124 and 125. Each test subject retains at least one observed visit.

Tucker3 uses ranks (4, 4, 3), ridge 0.001, HOSVD initialisation, tolerance 1e-6 and a maximum of 2,000 iterations. Reconstruction errors are calculated only for the artificially masked values. No outcome labels or PLS-DA models are used.

## Reproducing the benchmark

From the root of the paper repository, set the directory of the main analysis and a new output directory:

```bash
PRIMARY_RUN="/path/to/analysis_run"
BENCHMARK="/path/to/new_reconstruction_results"
FIGURES="/path/to/new_figure_S1"

Rscript --vanilla supplementary/missingness/real/01_gse68310_reconstruction.R "$PRIMARY_RUN" "$BENCHMARK"
Rscript --vanilla supplementary/missingness/real/02_gse68310_subject_analysis.R "$BENCHMARK"
Rscript --vanilla supplementary/missingness/real/03_plot_gse68310_reconstruction.R "$BENCHMARK" "$FIGURES"
```

`PRIMARY_RUN` is read-only. `BENCHMARK` and `FIGURES` must not exist before their respective commands. The plotting script saves `Figure_S1_pooled_RMSE_bootstrap.pdf` in `FIGURES`. Add `--preflight` to either command to check inputs before running. The reconstruction script also accepts an optional third argument pointing to saved historical benchmark results for a direct numerical comparison.

## Output files

The first script creates:

| File | Contents |
| --- | --- |
| `reconstruction_metrics.csv` | RMSE, Pearson correlation and squared errors for each fold, seed, masking level and method |
| `pooled_RMSE.csv` | RMSE pooled across folds and seeds for each masking level and method |
| `fit_summary.csv` | Tucker convergence information for each training fold |
| `fit_fold_1.rds` to `fit_fold_5.rds` | Saved training-mean and Tucker imputers |
| `masks/` | The 75 saved visit masks |
| `outer_folds.rds`, `input_provenance.rds` | Fold assignments and input checksums |
| `sessionInfo.txt`, `PASS.txt` | R environment and execution checks |

The second script reuses the saved imputers and adds:

| File | Contents |
| --- | --- |
| `subject_level_errors.csv` | Mean and Tucker reconstruction errors for each subject–mask evaluation |
| `subject_summary_overall.csv` | Overall subject-level summary |
| `subject_summary_by_level.csv` | Results by masking percentage |
| `subject_summary_by_remaining_visits.csv` | Results by number of visits still observed |

## Reference results

The verified execution produced five converged Tucker fits, 75 visit masks and 150 method-level evaluations. Pooled RMSE was:

| Masked visits | Training mean | Tucker3 |
| --- | ---: | ---: |
| 5% | 0.1162 | 0.1543 |
| 10% | 0.1157 | 0.1779 |
| 20% | 0.1145 | 0.1956 |
| 30% | 0.1126 | 0.1838 |
| 40% | 0.1133 | 0.1845 |

Tucker3 had lower individual RMSE in 497 of 896 subject–mask evaluations (55.47%), involving 108 distinct subjects. The 5% of evaluations with the largest Tucker squared errors contributed 53.23% of its total squared error.

Mean imputation therefore had lower pooled RMSE at all masking levels, although Tucker3 performed better for more than half of the individual evaluations. These results concern reconstruction of missing visits, not predictive classification.
