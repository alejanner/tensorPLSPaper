args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("How to use: Rscript 01_probe_to_gene.R PACKAGE_SOURCE_DIR SAVED_RUN_DIR NEW_OUTPUT_DIR",
    call. = FALSE)
}

source_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
run_dir <- normalizePath(path.expand(args[[2L]]), mustWork = TRUE)
requested_output <- path.expand(args[[3L]])
output_parent <- normalizePath(dirname(requested_output), mustWork = TRUE)
root <- file.path(output_parent, basename(requested_output))

inside_dir <- function(path, dir) identical(path, dir) || startsWith(path, paste0(dir, "/"))
if (inside_dir(output_parent, source_dir) || inside_dir(output_parent, run_dir)) {
  stop("Output must be outside the package source and saved analysis run.", call. = FALSE)
}
if (file.exists(root) || dir.exists(root)) {
  stop("Output directory already exists: ", root, call. = FALSE)
}

inputs <- c(
  file.path(run_dir, "GSE68310_analysis", c("descriptive_models.rds", "selected_probes.csv")),
  file.path(run_dir, "GSE68310_data", "X_incomplete.rds"),
  file.path(run_dir, "GSE68310_data/downloads/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"),
  file.path(run_dir, "GSE182616_analysis", c("descriptive_models.rds", "selected_probes.csv", "X_incomplete.rds")),
  file.path(run_dir, "GSE182616_data", "feature_annotation.csv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing analysis input(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required.", call. = FALSE)
}
pkgload::load_all(source_dir, export_all = FALSE, quiet = TRUE)
stopifnot(normalizePath(getNamespaceInfo(asNamespace("TensorPLS"), "path")) == normalizePath(source_dir))
if (!dir.create(root, showWarnings = FALSE)) {
  stop("Cannot create output directory: ", root, call. = FALSE)
}
write_table <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
parse_ids <- function(value) {
  if (is.na(value) || trimws(value) %in% c("", "NA", "---", "-", "0", "NULL")) return(character())
  unique(trimws(strsplit(value, "[[:space:];,|/]+")[[1L]]))
}
provenance <- list()
for (dataset in c("GSE68310", "GSE182616")) {
  out <- file.path(root, "results", dataset)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out, "mapping_summary.csv"))) stop("Mapping outputs already exist: ", out)
  analysis_dir <- file.path(run_dir, paste0(dataset, "_analysis"))
  model_path <- file.path(analysis_dir, "descriptive_models.rds")
  selected_path <- file.path(analysis_dir, "selected_probes.csv")
  tensor_path <- file.path(run_dir, if (dataset == "GSE68310") "GSE68310_data" else "GSE182616_analysis", "X_incomplete.rds")
  saved <- readRDS(model_path)
  selected <- read.csv(selected_path, colClasses = "character")$feature
  X <- readRDS(tensor_path)
  features <- dimnames(X)[[2L]]
  times <- dimnames(X)[[3L]]
  expected <- if (dataset == "GSE68310") c(109L, 47162L, 4L) else c(66L, 50599L, 6L)
  stopifnot(identical(dim(X), expected), !anyDuplicated(features),
    identical(features, saved$model_all$feature_names),
    identical(times, saved$model_all$time_names),
    identical(selected, saved$selection$selected_features),
    identical(selected, saved$model_selected$feature_names), all(selected %in% features),
    length(selected) == if (dataset == "GSE68310") 1132L else 2175L)
  rm(X)
  gc()
  if (dataset == "GSE68310") {
    platform <- "GPL10558"
    annotation_path <- file.path(run_dir, "GSE68310_data/downloads/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz")
    con <- gzfile(annotation_path, "rt")
    lines <- readLines(con, warn = FALSE)
    close(con)
    begin <- match("[Probes]", lines)
    end <- match("[Controls]", lines)
    annotation <- read.delim(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"),
      quote = "", comment.char = "", colClasses = "character", check.names = FALSE)
    stopifnot(digest::digest(file = annotation_path, algo = "sha256") ==
      "932e2f110bd30ac70ef4d7dbac2f87f8088f982bf9aa596732c7477ab50d43c2")
    controls_end <- match("[Columns]", lines)
    controls <- read.delim(text = paste(lines[(end + 1L):(controls_end - 1L)], collapse = "\n"),
      quote = "", comment.char = "", colClasses = "character", check.names = FALSE)
    ercc <- unique(controls$Probe_Id[grepl("^ERCC-", controls$Reporter_Group_Name)])
    stopifnot(length(ercc) == 92L, !any(features %in% ercc))
    key <- "Probe_Id"
    gene_column <- "Entrez_Gene_ID"
    symbol_column <- "Symbol"
    annotation_version <- "Illumina ContentVersion 4.0; 15/4/2010; R2_15002873_B"
  } else {
    platform <- "GPL17077"
    annotation_path <- file.path(run_dir, "GSE182616_data/feature_annotation.csv")
    annotation <- read.csv(annotation_path, colClasses = "character", check.names = FALSE)
    key <- "ID"
    gene_column <- "LOCUSLINK_ID"
    symbol_column <- "GENE_SYMBOL"
    stopifnot(all(annotation$retained[match(features, annotation$ID)] == "TRUE"))
    annotation_version <- "Local GPL17077-derived table; SOFT last update 2018-01-09, hg19:GRCh37:Feb2009"
  }
  stopifnot(!anyDuplicated(annotation[[key]]))
  index <- match(features, annotation[[key]])
  raw <- annotation[index, , drop = FALSE]
  rownames(raw) <- NULL
  ids <- lapply(raw[[gene_column]], parse_ids)
  valid <- vapply(ids, function(x) length(x) == 1L && grepl("^[1-9][0-9]*$", x), logical(1))
  ambiguous <- lengths(ids) > 1L
  status <- ifelse(is.na(index), "probe_not_in_annotation",
    ifelse(ambiguous, "multiple_gene_identifiers",
      ifelse(lengths(ids) == 0L, "no_gene_identifier",
        ifelse(valid, "unambiguous_annotation_ID", "invalid_gene_identifier"))))
  mapping <- data.frame(probe_id = features, selected = features %in% selected,
    gene_id = vapply(seq_along(ids), function(i) if (valid[i]) ids[[i]] else NA_character_, character(1)),
    gene_symbol = raw[[symbol_column]], mapping_status = status,
    current_Entrez_validation = "PENDING_org.Hs.eg.db", stringsAsFactors = FALSE)
  names(raw) <- paste0("annotation_", names(raw))
  mapping <- cbind(mapping, raw)
  stopifnot(nrow(mapping) == length(features), !anyDuplicated(mapping$probe_id))
  selected_mapping <- mapping[match(selected, mapping$probe_id), , drop = FALSE]
  background <- unique(mapping$gene_id[valid])
  genes <- unique(selected_mapping$gene_id[!is.na(selected_mapping$gene_id)])
  stopifnot(all(genes %in% background), identical(selected_mapping$probe_id, selected))
  write_table(mapping, file.path(out, "probe_gene_mapping.csv"))
  write_table(selected_mapping, file.path(out, "selected_probe_gene_mapping.csv"))
  write_table(mapping[!valid & !ambiguous, ], file.path(out, "unmapped_probes.csv"))
  write_table(mapping[ambiguous, ], file.path(out, "ambiguous_probes.csv"))
  write_table(data.frame(gene_id = genes), file.path(out, "selected_genes.csv"))
  write_table(data.frame(gene_id = background), file.path(out, "background_genes.csv"))
  categories <- interaction(mapping$selected, mapping$mapping_status, drop = TRUE)
  examples <- unique(unlist(lapply(split(seq_len(nrow(mapping)), categories), head, 5L)))
  examples <- unique(c(examples, seq_len(min(20L, nrow(mapping)))))
  write_table(mapping[examples, ], file.path(out, "mapping_audit_examples.csv"))
  counts <- table(mapping$gene_id[valid])
  summary <- data.frame(dataset, platform, tensor_probes = length(features),
    selected_probes = length(selected), annotation_matches = sum(!is.na(index)),
    unambiguous_annotation_probes = sum(valid), unmapped_probes = sum(!valid & !ambiguous),
    ambiguous_probes = sum(ambiguous), background_genes = length(background),
    genes_with_multiple_probes = sum(counts > 1L),
    selected_mapped_probes = sum(!is.na(selected_mapping$gene_id)),
    selected_mapping_percent = 100 * mean(!is.na(selected_mapping$gene_id)),
    selected_unique_genes = length(genes), selected_genes_outside_background = length(setdiff(genes, background)))
  write_table(summary, file.path(out, "mapping_summary.csv"))
  vips <- TensorPLS::tensorpls_vips(saved$model_all)
  temporal <- selected_mapping[, c("probe_id", "gene_id", "gene_symbol", "mapping_status")]
  temporal <- cbind(temporal, as.data.frame(vips$M2_mean[selected, , drop = FALSE], check.names = FALSE))
  write_table(temporal, file.path(out, "selected_probe_M2_temporal_profiles.csv"))
  provenance[[dataset]] <- data.frame(dataset, platform, annotation_path, annotation_version,
    source_url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", platform),
    annotation_sha256 = digest::digest(file = annotation_path, algo = "sha256"))
  print(summary, row.names = FALSE)
  rm(saved, vips, mapping, raw, annotation)
  gc()
}
write_table(do.call(rbind, provenance), file.path(root, "annotation_provenance.csv"))
writeLines(capture.output(sessionInfo()), file.path(root, "mapping_sessionInfo.txt"))
