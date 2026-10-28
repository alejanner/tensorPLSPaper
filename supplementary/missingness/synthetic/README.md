# Synthetic missing-visit reconstruction

This benchmark compares training-set mean imputation and Tucker3 reconstruction on simulated longitudinal data. Three scenarios are considered, with three independently generated datasets per scenario. Each dataset contains 109 subjects, 1,000 features and four time points. No external expression data are required.

## Simulation and evaluation

Scenario A has strong subject-level variation shared across features, B has weaker individual variation and more noise, and C has more heterogeneous trajectories across 20 feature modules. The simulation parameters are defined in `01_synthetic_reconstruction.R`.

Each dataset is split into five folds. Both imputers are fitted on training subjects only, and 10%, 20%, 30% or 40% of visits are masked in held-out subjects using three masking seeds. Each subject keeps at least one observed visit.

Tucker3 uses ranks (4, 4, 3), ridge 0.001, HOSVD initialisation, tolerance 1e-6 and at most 2,000 iterations. Reconstruction errors are calculated against both the simulated noisy measurements and the underlying noise-free signal. Comparisons between the two methods use only masks where Tucker reconstruction succeeded.

## Reproducing the benchmark

Install [TensorPLS](https://github.com/alejanner/TensorPLS) and its dependencies, including `ggplot2` for the figure. The benchmark uses two internal TensorPLS imputation helpers, so use the package version associated with this paper.

From the root of the paper repository, choose new output directories:

```bash
SYNTHETIC="/path/to/new_synthetic_results"
FIGURES="/path/to/new_synthetic_figures"

Rscript --vanilla supplementary/missingness/synthetic/01_synthetic_reconstruction.R full "$SYNTHETIC"
Rscript --vanilla supplementary/missingness/synthetic/02_plot_synthetic_reconstruction.R "$SYNTHETIC" "$FIGURES"
```

Neither output directory may exist before its command runs. The benchmark generates all simulated datasets; no prepared input tensor is needed. An optional `pilot` mode tests one fold of the first replicate in each scenario before the full run. To load a local TensorPLS source checkout instead of the installed package, set `TENSORPLS_SOURCE_DIR`.

## Output files

The benchmark creates:

| File | Contents |
| --- | --- |
| `datasets/` | Simulated datasets, latent signal and fold assignments |
| `fits/`, `masks/` | Saved imputers and visit masks |
| `fit_summary.csv` | Convergence and fitting information |
| `reconstruction_metrics.csv` | Errors for each scenario, dataset, fold, seed and masking level |
| `paired_RMSE_by_level.csv` | RMSE pooled over valid paired cases for each scenario and masking level |
| `paired_RMSE_by_replicate.csv` | RMSE for each of the three datasets within each scenario |
| `subject_level_errors.csv`, `subject_summary_*.csv` | Subject-level reconstruction errors and summaries |
| `sessionInfo.txt` | R environment |

The figure script reads the two `paired_RMSE` tables and creates `Figure_S2_pooled_and_three_datasets.pdf` and `.png`, alongside the plotting data and input checksums. It does not repeat the simulations or imputer fits.

Figure S2 shows the three independent dataset curves for each scenario and method, with thicker lines for the pooled RMSE. Pooled RMSE is calculated from the matched valid cases, not as the arithmetic mean of the three dataset RMSEs. The y-axis scale differs between scenarios.

In the reference full run, 43 of the 45 Tucker fits converged. Failed cases in scenario B were excluded from the paired comparison and are not treated as successful reconstructions. These results assess missing-visit reconstruction. Not classification performance.
