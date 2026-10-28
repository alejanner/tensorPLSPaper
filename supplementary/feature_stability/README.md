# Feature-selection stability

This workflow compares four VIP-based selection strategies (V, M1, M2 and V ∩ M2) in GSE68310 and GSE182616. `generate_nested_cv.R` runs the six supplementary nested CVs; the two V ∩ M2 results come from the main analyses. The stability scripts summarise and verify these saved results without fitting models.

## 1. Generate the nested CV results

Run `run_all.sh` first. Set `MAIN` to its completed output directory, and choose new paths for `VIP` and `TABLES`:

```bash
P="/absolute/path/to/tensorPLSPaper"
SOURCE="/absolute/path/to/TensorPLS"
MAIN="/absolute/path/to/completed_main_run"
VIP="/absolute/path/to/new_vip_benchmark"
TABLES="/absolute/path/to/new_stability_tables"

Rscript --vanilla "$P/supplementary/feature_stability/generate_nested_cv.R" \
  "$SOURCE" "$MAIN" "$VIP"
```

The generator reuses the saved tensors, outcomes and subject-level folds for the 5 × 5 nested CV. It evaluates H = 1–10 with a fixed 99th-percentile VIP threshold and Q²Y tuning, using Tucker3 for GSE68310 and training-mean imputation for GSE182616. The six RDS files are saved under `VIP/runs/<dataset>/<strategy>/nested_cv_result.rds`. Both `VIP` and `TABLES` must be new directories.

## 2. Run the analysis

Generate the stability tables from the eight saved nested-CV results:

```bash
Rscript --vanilla "$P/supplementary/feature_stability/01_feature_stability.R" \
  --input-root "$MAIN" \
  --benchmark-root "$VIP" \
  --output-dir "$TABLES"
```

The script checks that all five outer folds completed successfully, that probe identifiers are valid, and that the saved strategies use matching cross-validation partitions and settings. It stops with an error if an input is missing or inconsistent.

## 3. Verify the results

Using the same directories and output path, run:

```bash
Rscript --vanilla "$P/supplementary/feature_stability/02_verify_feature_stability.R" \
  --input-root "$MAIN" \
  --benchmark-root "$VIP" \
  --output-dir "$TABLES"
```

This second script independently recalculates probe-selection counts, frequencies, and pairwise Jaccard indices from the saved fold selections. It compares them with the generated CSV tables, checks available stability values stored in the RDS objects, and checks the numerical reference results below. It does not modify the input files or write additional output files.

A successful verification ends with:

```text
PASS: all eight dataset-strategy combinations. No models fitted and no files written.
```

If verification fails, inspect the error message and check that the supplied RDS files come from the matching analyses.

## 4. Output tables

The analysis creates five CSV files under the chosen output directory:

| File | Contents |
| --- | --- |
| `strategy_stability_summary.csv` | Median number of selected probes and mean pairwise Jaccard for each dataset and strategy; union and five-fold intersection sizes for V ∩ M2 |
| `strategy_fold_counts.csv` | Number of selected probes in each of the five outer folds |
| `strategy_jaccard_pairs.csv` | Jaccard indices for all ten distinct pairs of outer folds, for each dataset and strategy |
| `primary_selection_frequency.csv` | Selection count (0–5) and relative frequency for every probe under V ∩ M2 |
| `primary_frequency_distribution.csv` | Number of probes selected in exactly 0, 1, 2, 3, 4, or 5 folds under V ∩ M2 |

All calculations use **probe identifiers**, not gene-level aggregations. Probes not selected in any fold are retained in the primary frequency table with a count of zero.

## 5. Expected results

The table below shows the values reported for the saved benchmark runs. Jaccard values are rounded to six decimal places.

| Dataset | Strategy | Median selected probes | Mean pairwise Jaccard |
| --- | --- | ---: | ---: |
| GSE68310 | V | 1,485 | 0.564644 |
| GSE68310 | M1 | 500 | 0.664511 |
| GSE68310 | M2 | 1,304 | 0.566827 |
| GSE68310 | V ∩ M2 | 1,122 | 0.583822 |
| GSE182616 | V | 2,684 | 0.247875 |
| GSE182616 | M1 | 542 | 0.322826 |
| GSE182616 | M2 | 2,538 | 0.258883 |
| GSE182616 | V ∩ M2 | 2,194 | 0.243913 |

For the primary **V ∩ M2** strategy, the additional overlap counts are:

| Dataset | Probes selected in at least one fold | Probes selected in all five folds |
| --- | ---: | ---: |
| GSE68310 | 1,847 | 629 |
| GSE182616 | 6,084 | 368 |

The reference values are used **only by the verification script**, not to populate the output tables.

## 6. What the stability measures mean

For selected probe sets $S_a$ and $S_b$, the Jaccard index is:

$$
J(S_a,S_b) = \frac{|S_a \cap S_b|}{|S_a \cup S_b|}
$$

A value of 1 indicates identical selected sets; 0 indicates no overlap. The reported mean is calculated across the **ten distinct pairs** formed by five outer folds. For each probe, selection frequency is the number of folds in which it was selected divided by five.

These measures describe **agreement between feature sets selected from overlapping training samples**.
