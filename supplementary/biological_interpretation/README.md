# Biological interpretation

Run the scripts in order, using a new output directory for each step. The saved TensorPLS models and previous analysis results are read-only inputs.

## 1. Probe-to-gene mapping

`01_probe_to_gene.R` maps the selected probes for GSE68310 and GSE182616 to Entrez IDs, builds the measured-gene background and exports temporal M2 values from the all-feature model.

**Inputs:** a local copy of [TensorPLS](https://github.com/alejanner/TensorPLS) and a completed analysis run from [tensorPLSPaper](https://github.com/alejanner/tensorPLSPaper), containing:

```text
GSE68310_data/X_incomplete.rds
GSE68310_data/downloads/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz
GSE68310_analysis/descriptive_models.rds
GSE68310_analysis/selected_probes.csv
GSE182616_data/feature_annotation.csv
GSE182616_analysis/X_incomplete.rds
GSE182616_analysis/descriptive_models.rds
GSE182616_analysis/selected_probes.csv
```

The main preparation scripts obtain or create the annotation files. No annotations are downloaded during mapping. R requires `pkgload`, `digest` and the TensorPLS dependencies.

```bash
Rscript --vanilla supplementary/biological_interpretation/01_probe_to_gene.R \
  /path/to/TensorPLS \
  /path/to/saved_run \
  /path/to/new_mapping_output
```

**Output:** `new_mapping_output/results/{GSE68310,GSE182616}/` contains `probe_gene_mapping.csv`, `selected_probe_gene_mapping.csv`, `selected_genes.csv`, `background_genes.csv`, `unmapped_probes.csv`, `ambiguous_probes.csv`, `mapping_audit_examples.csv`, `mapping_summary.csv`, and `selected_probe_M2_temporal_profiles.csv`. The root contains `annotation_provenance.csv` and `mapping_sessionInfo.txt`.

## 2. Entrez ID validation

`01b_validate_entrez.R` checks the IDs from step 1 against `org.Hs.eg.db`. It keeps recognized Entrez IDs, saves excluded IDs and checks the expected gene counts.

**Inputs:** `selected_genes.csv` and `background_genes.csv` for both datasets, in the mapping output above. Requires Bioconductor packages `AnnotationDbi` and `org.Hs.eg.db` **3.23.1** for exact reproduction.

```bash
Rscript --vanilla supplementary/biological_interpretation/01b_validate_entrez.R \
  /path/to/new_mapping_output \
  /path/to/new_validation_output
```

Optionally, add `/path/to/previous_biology_results` as a third argument to compare the new selected/background ID sets against existing `results/<dataset>/validated_*.csv` files.

**Output:** `new_validation_output/{GSE68310,GSE182616}/` contains `validated_selected_genes.csv`, `validated_background_genes.csv`, `unrecognised_selected_Entrez.csv`, and `unrecognised_background_Entrez.csv`. The output root contains `validation_summary.csv` and `sessionInfo.txt`.

The expected validated counts are **921 selected / 23,306 background** for GSE68310 and **1,651 selected / 20,859 background** for GSE182616. The scripts do not overwrite existing output directories.

## 3. Functional enrichment

`02_enrichment.R` calls `TensorPLS::tensorpls_enrichment()` for both datasets. The enrichment methods are implemented in the package, not repeated here. This script only reads the gene lists, runs the function and saves its results.

**Study background and Hallmark annotation.** `01_probe_to_gene.R` builds the study background from all measured probes with an unambiguous Entrez mapping, not only probes selected by TensorPLS. `01b_validate_entrez.R` retains recognized IDs (23,306 background genes for GSE68310; 20,859 for GSE182616). Independently, `00_generate_hallmark.R` retrieves the 50 human Hallmark gene sets from MSigDB 2026.1.Hs using `msigdbr` 26.1.1 and creates the TERM2GENE table (7,322 term–Entrez gene pairs). This table defines gene-set membership; it does **not** generate the study background. The Hallmark collection is from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) (CC BY 4.0).

`02_enrichment.R` supplies the validated selected genes as `genes`, the validated study background as `universe`, and TERM2GENE as the Hallmark annotation. The **effective Hallmark universe** is the intersection of the study background with genes represented in TERM2GENE:

$$
B_{\mathrm{eff}} = B_{\mathrm{study}} \cap G_{\mathrm{Hallmark}}.
$$

Here, $G_{\mathrm{Hallmark}}$ is the union of genes annotated in the Hallmark gene sets. The effective universe contains 4,378 genes for GSE68310 and 4,373 for GSE182616; GO BP and Reactome have their own annotation-specific effective universes.

For each Hallmark gene set $T$, the one-sided hypergeometric test evaluates whether selected genes are over-represented:

$$
P(X \ge k) = \sum_{i=k}^{\min(M,n)}
\frac{\binom{M}{i}\binom{N-M}{n-i}}{\binom{N}{n}}.
$$

Here, $N=|B_{\mathrm{eff}}|$ is the effective background size, $M=|T\cap B_{\mathrm{eff}}|$ is the number of background genes in the gene set, $n=|S\cap B_{\mathrm{eff}}|$ is the number of selected genes in the effective background, and $k=|S\cap T\cap B_{\mathrm{eff}}|$ is their overlap, where $S$ denotes the validated selected genes. TERM2GENE therefore determines gene-set membership and contributes to the effective background used by the test. P-values are adjusted across tested gene sets using the Benjamini–Hochberg procedure.

To independently regenerate and verify the frozen Hallmark table, run from the paper repository root (requires `msigdbr` 26.1.1):

```bash
Rscript --vanilla supplementary/biological_interpretation/00_generate_hallmark.R \
  /tmp/Hallmark_regenerated.csv \
  supplementary/biological_interpretation/data/Hallmark_TERM2GENE.csv
```

**Inputs:** the TensorPLS package source, the `validation_results/` directory from step 2, and the frozen `data/Hallmark_TERM2GENE.csv` (MSigDB 2026.1.Hs). R needs `pkgload` and the optional TensorPLS enrichment dependencies (`AnnotationDbi`, `org.Hs.eg.db`, `clusterProfiler`, `ReactomePA`, and their dependencies).

```bash
Rscript --vanilla supplementary/biological_interpretation/02_enrichment.R \
  /path/to/TensorPLS \
  /path/to/validation_results \
  supplementary/biological_interpretation/data/Hallmark_TERM2GENE.csv \
  'MSigDB 2026.1.Hs' \
  /path/to/new_enrichment_output
```

**Output:** `new_enrichment_output/results/{GSE68310,GSE182616}/` contains the complete and significant GO BP and Reactome tables (`*_results.csv`, `*_significant_FDR005.csv`), database-specific effective universes (`*_actual_universe.csv`), fitted RDS objects (`*_enrichment.rds`), and the full `tensorpls_enrichment.rds` result. Hallmark files are under each dataset's `Hallmark/` directory, including `Hallmark_significant_gene_overlaps.csv` for temporal interpretation. At the output root, `enrichment_summary.csv`, `input_md5.csv`, `sessionInfo.txt` and `STATUS.txt` record the results and input provenance.

For these reference datasets, the expected numbers of significant terms (BH-adjusted P < 0.05) are 306/55/7 for GSE68310 and 0/2/1 for GSE182616 (GO BP/Reactome/Hallmark). The script stops if these counts or the number of tested terms differ.

## 4. Interferon and cell-cycle temporal profiles

`04_hallmark_IFN_cellcycle_M2.R` reproduces Figure 2C for GSE68310. It compares temporal M2 profiles of selected probes associated with Hallmark interferon alpha/gamma responses and cell-cycle gene sets (E2F targets and G2M checkpoint). Genes shared between the two groups are excluded, and the remaining selected probes form the reference group.

The script uses M2 values from the full-cohort all-feature model. It does not refit models or rerun enrichment.

**Inputs:** the mapping output from step 1 and the enrichment output from step 3. Requires `ggplot2`.

```bash
Rscript --vanilla supplementary/biological_interpretation/04_hallmark_IFN_cellcycle_M2.R \
  /path/to/new_mapping_output \
  /path/to/new_enrichment_output \
  /path/to/new_figure_output
```

The output directory must not already exist.

**Output:** `M2_IFN_vs_cellcycle.pdf` and `M2_IFN_vs_cellcycle.png` reproduce Figure 2C. The script also saves `Hallmark_terms.csv`, `gene_groups.csv`, `probe_groups.csv`, `temporal_summary.csv` and `sessionInfo.txt`.

The temporal profiles describe model-assigned variable importance, not gene expression or pathway activation.
