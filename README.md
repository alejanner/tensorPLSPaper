# TensorPLS — Paper reproducibility

This repository contains the scripts used for the main and supplementary analyses in the TensorPLS manuscript. The R package is available separately at [TensorPLS](https://github.com/alejanner/TensorPLS).

The analyses use two public longitudinal transcriptomic datasets: [GSE68310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68310) and [GSE182616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182616).

## Reproducibility map

This repository contains:

```text
tensorPLSPaper/
├── run_all.sh                                → run the four main scripts
├── 01_prepare_GSE68310.R                   → prepare the GSE68310 tensor
├── 02_analyze_GSE68310.R                   → nested CV and descriptive models
├── 03_prepare_GSE182616.R                  → prepare GSE182616 expression and metadata
├── 04_analyze_GSE182616.R                  → tensor, nested CV and descriptive models
│
├── images/                                 → figures in the main paper
│   ├── 07_gse68310_panels_AB_paper.R
│   │   ├── Figure 2A → GSE68310_Scores_Paper.pdf
│   │   └── Figure 2B → GSE68310_Temporal_VIP_Paper.pdf
│   ├── 08_gse182616_panels_AB_paper.R
│   │   ├── Figure 3A → GSE182616_Scores_Paper.pdf
│   │   └── Figure 3B → GSE182616_Temporal_VIP_Paper.pdf
│   └── README.md                           → figure inputs and commands
│
└── supplementary/
    ├── biological_interpretation/          → mapping, enrichment and Figure 2C
    ├── feature_stability/                  → VIP strategy comparison (Table S6)
    ├── missingness/
    │   ├── real/                           → Figure S1
    │   └── synthetic/                      → Figure S2
    └── sensitivity/
        ├── run_gse68310_sensitivity.R
        └── README.md                       → subgroup and cross-season analyses
```

## Run the main analyses

You need R, Git, Bash, `gzip`, internet access for the GEO downloads, and the R packages required by TensorPLS.

Clone both repositories into the same parent directory:

```bash
mkdir -p ~/TensorPLS_reproduction
cd ~/TensorPLS_reproduction
git clone --branch v0.1.0 https://github.com/alejanner/TensorPLS.git
git clone https://github.com/alejanner/tensorPLSPaper.git
```

Install TensorPLS and the dependencies used by the main reproduction scripts:

```bash
Rscript -e 'install.packages(c("remotes", "pkgload", "data.table", "digest"))'
Rscript -e 'remotes::install_local("TensorPLS", upgrade = "never")'
```

Run the four preparation and analysis scripts in order:

```bash
cd ~/TensorPLS_reproduction/tensorPLSPaper
bash run_all.sh
```

The runner downloads and prepares the GEO data, performs subject-level nested cross-validation, and fits the full-cohort descriptive models. Each run creates a new directory under `outputs/` with the data, results, and logs.

The scripts can also be read separately: [GSE68310 preparation](01_prepare_GSE68310.R), [GSE68310 analysis](02_analyze_GSE68310.R), [GSE182616 preparation](03_prepare_GSE182616.R), and [GSE182616 analysis](04_analyze_GSE182616.R).

## Main figures and tables

| Figure or table | What it shows | Script |
| --- | --- | --- |
| Figure 1 | TensorPLS workflow | Schematic; no R plotting script |
| Figure 2A–B | GSE68310 scores and temporal M2 profile | [GSE68310 panels](images/07_gse68310_panels_AB_paper.R) |
| Figure 2C | Interferon and cell-cycle temporal M2 profiles | [Hallmark M2 comparison](supplementary/biological_interpretation/04_hallmark_IFN_cellcycle_M2.R) |
| Figure 3A–B | GSE182616 scores and temporal M2 profile | [GSE182616 panels](images/08_gse182616_panels_AB_paper.R) |
| Predictive performance table | Q²Y, AUROC and balanced accuracy | [GSE68310 analysis](02_analyze_GSE68310.R) and [GSE182616 analysis](04_analyze_GSE182616.R) |

See [images/README.md](images/README.md) for the figure inputs, commands, and output filenames. These scripts use saved fitted models rather than rerunning the analyses. Figures 2 and 3 describe full-cohort models; they are not held-out predictive results.

## Supplementary analyses

| Analysis | Code and instructions |
| --- | --- |
| Figure S1 — Real-data missing-visit reconstruction and bootstrap intervals | [Real-data missingness README](supplementary/missingness/real/README.md), including [Figure S1 script](supplementary/missingness/real/03_plot_gse68310_reconstruction.R) |
| Figure S2 — Synthetic missing-visit reconstruction | [Synthetic missingness README](supplementary/missingness/synthetic/README.md), including [Figure S2 script](supplementary/missingness/synthetic/02_plot_synthetic_reconstruction.R) |
| VIP strategy comparison and feature-selection stability | [Feature stability README](supplementary/feature_stability/README.md) |
| Probe-to-gene mapping, enrichment, and temporal interpretation | [Biological interpretation README](supplementary/biological_interpretation/README.md) |
| GSE68310 sensitivity analyses | [Sensitivity README](supplementary/sensitivity/README.md) |

To run the GSE68310 sensitivity analysis on a completed main run:

```bash
RUN="/path/to/outputs/run_YYYYMMDD_HHMMSS"
TENSORPLS_SOURCE_DIR="../TensorPLS" \
GSE68310_DATA="$RUN/GSE68310_data" \
GSE68310_RESULTS="$RUN/GSE68310_analysis" \
GSE68310_SENSITIVITY_OUTPUT="/path/to/new_sensitivity_results" \
  Rscript --vanilla supplementary/sensitivity/run_gse68310_sensitivity.R
```

## Reference results

| Dataset | Q²Y | AUROC | Balanced accuracy | Full-cohort components (H) | Selected probes |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSE68310 | 0.6697 | 0.9715 | 0.9309 | 5 | 1,132 |
| GSE182616 | 0.2939 | 0.8272 | 0.7564 | 2 | 2,175 |

Predictive metrics come from outer held-out predictions in 5 × 5 subject-level nested cross-validation. Component counts and selected-probe counts come from separate full-cohort descriptive fits. Each main analysis saves `reference_comparison.csv` and `COMPLETED.txt` to compare its results with the manuscript values.
