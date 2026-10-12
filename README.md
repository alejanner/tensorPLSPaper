# TensorPLS — Paper Reproducibility

This repository contains four R scripts and a Bash runner to reproduce the main analyses of the TensorPLS manuscript on two public longitudinal transcriptomic datasets. Both analyses use 5 × 5 nested cross-validation.

## Installation and execution

Requirements: R, Git, Bash, `gzip`, internet access for GEO downloads, and the R packages used by TensorPLS.

Clone **both repositories into the same parent directory**:

```bash
mkdir -p ~/TensorPLS_reproduction
cd ~/TensorPLS_reproduction
git clone https://github.com/alejanner/TensorPLS.git
git clone https://github.com/alejanner/tensorPLSPaper.git
```

The directory structure should be:

```text
TensorPLS_reproduction/
├── TensorPLS/
└── tensorPLSPaper/
    ├── run_all.sh
    ├── 01_prepare_GSE68310.R
    ├── 02_analyze_GSE68310.R
    ├── 03_prepare_GSE182616.R
    └── 04_analyze_GSE182616.R
```

Install the R dependencies:

```bash
Rscript -e 'install.packages(c("remotes", "pkgload", "data.table", "digest"))'
Rscript -e 'remotes::install_github("alejanner/TensorPLS", upgrade = "never")'
```

Run the analyses:

```bash
cd ~/TensorPLS_reproduction/TensorPLS_paper_reproducibility
chmod +x run_all.sh
./run_all.sh
```

The runner loads TensorPLS from the adjacent source checkout, prepares the two GEO datasets and executes their analyses in order. Each execution writes its results and logs to a new directory under `outputs/`.

## Manuscript results

| Dataset | Q²Y | AUROC | Balanced accuracy | Full-cohort components (H) | Selected probes |
|---|---:|---:|---:|---:|---:|
| [GSE68310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68310) | 0.6697 | 0.9715 | 0.9309 | 5 | 1,132 |
| [GSE182616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182616) | 0.2939 | 0.8272 | 0.7564 | 2 | 2,175 |

Predictive metrics are based on outer held-out predictions. Component counts and selected probe counts refer to full-cohort models. Each analysis saves `reference_comparison.csv` and `COMPLETED.txt` to indicate whether its results match the reported values.

