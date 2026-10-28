#!/usr/bin/env Rscript

# Reproduce GO BP, Reactome and Hallmark enrichment using TensorPLS.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript 02_enrichment.R PACKAGE_SOURCE VALIDATION_DIR",
    "HALLMARK_TERM2GENE_CSV HALLMARK_VERSION NEW_OUTPUT_DIR"
  ), call. = FALSE)
}

source_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
validation_dir <- normalizePath(path.expand(args[[2L]]), mustWork = TRUE)
hallmark_file <- normalizePath(path.expand(args[[3L]]), mustWork = TRUE)
hallmark_version <- args[[4L]]
requested_out <- path.expand(args[[5L]])
out_parent <- normalizePath(dirname(requested_out), mustWork = TRUE)
out_dir <- file.path(out_parent, basename(requested_out))

inside <- function(path, directory) {
  identical(path, directory) || startsWith(path, paste0(directory, "/"))
}
if (!nzchar(hallmark_version) || file.exists(out_dir)) {
  stop("Supply a Hallmark release label and a new output directory.", call. = FALSE)
}
if (any(vapply(c(source_dir, validation_dir, dirname(hallmark_file)),
               function(x) inside(out_dir, x), logical(1)))) {
  stop("Output must be outside all input directories.", call. = FALSE)
}

read_ids <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  tab <- read.csv(path, colClasses = "character", check.names = FALSE)
  if (ncol(tab) != 1L || !nrow(tab)) stop("Expected one nonempty gene column: ", path)
  ids <- tab[[1L]]
  if (anyNA(ids) || any(!grepl("^[1-9][0-9]*$", ids)) || anyDuplicated(ids)) {
    stop("Invalid or duplicated Entrez IDs in ", path, call. = FALSE)
  }
  ids
}
write_ids <- function(ids, path) {
  write.csv(data.frame(gene_id = as.character(ids)), path, row.names = FALSE)
}

datasets <- c("GSE68310", "GSE182616")
input_files <- hallmark_file
for (dataset in datasets) {
  input_files <- c(input_files, file.path(validation_dir, dataset,
    c("validated_selected_genes.csv", "validated_background_genes.csv")))
}
if (any(!file.exists(input_files))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "),
       call. = FALSE)
}

raw <- read.csv(hallmark_file, colClasses = "character", check.names = FALSE)
term_name <- intersect(c("term", "gs_name", "pathway", "ID"), names(raw))
gene_name <- intersect(c("gene", "ncbi_gene", "gene_id", "ENTREZID", "entrez_gene"), names(raw))
if (length(term_name) && length(gene_name)) {
  hallmark <- raw[, c(term_name[[1L]], gene_name[[1L]]), drop = FALSE]
} else if (ncol(raw) == 2L) {
  hallmark <- raw
} else {
  stop("Cannot identify Hallmark term and Entrez columns.", call. = FALSE)
}
names(hallmark) <- c("term", "gene")
hallmark <- unique(hallmark)
if (!nrow(hallmark) || anyNA(hallmark) ||
    any(!grepl("^HALLMARK_", hallmark$term)) ||
    any(!grepl("^[1-9][0-9]*$", hallmark$gene)) ||
    length(unique(hallmark$term)) != 50L) {
  stop("Invalid Hallmark TERM2GENE table or unexpected number of gene sets.", call. = FALSE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.", call. = FALSE)
pkgload::load_all(source_dir, export_all = FALSE, quiet = TRUE)
stopifnot(identical(normalizePath(getNamespaceInfo("TensorPLS", "path")), source_dir))

expected <- list(
  GSE68310 = list(GO_BP = c(4629L, 306L, 17807L),
                  Reactome = c(1390L, 55L, 11168L),
                  Hallmark = c(50L, 7L, 4378L)),
  GSE182616 = list(GO_BP = c(5178L, 0L, 16877L),
                   Reactome = c(1530L, 2L, 11089L),
                   Hallmark = c(50L, 1L, 4373L))
)

if (!dir.create(out_dir)) stop("Cannot create: ", out_dir, call. = FALSE)
writeLines("RUNNING", file.path(out_dir, "STATUS.txt"))
summary_rows <- list()

for (dataset in datasets) {
  input_dir <- file.path(validation_dir, dataset)
  genes <- read_ids(file.path(input_dir, "validated_selected_genes.csv"))
  universe <- read_ids(file.path(input_dir, "validated_background_genes.csv"))
  if (!all(genes %in% universe)) stop("Selected genes outside universe: ", dataset)

  enrichment <- TensorPLS::tensorpls_enrichment(
    genes = genes, universe = universe, id_type = "ENTREZID",
    databases = c("GO_BP", "Reactome", "Hallmark"),
    hallmark_term2gene = hallmark, hallmark_version = hallmark_version
  )
  if (!setequal(enrichment$converted$genes, genes) ||
      !setequal(enrichment$converted$universe, universe)) {
    stop("Entrez recognition differs from the validation step: ", dataset)
  }

  base <- file.path(out_dir, "results", dataset)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  saveRDS(enrichment, file.path(base, "tensorpls_enrichment.rds"))

  for (database in names(enrichment$results)) {
    result <- enrichment$results[[database]]
    counts <- c(nrow(result$full), nrow(result$significant),
                length(result$effective_universe))
    if (result$status != "tested" ||
        !identical(as.integer(counts), expected[[dataset]][[database]])) {
      stop(dataset, " / ", database, ": unexpected tested/significant/universe counts: ",
           paste(counts, collapse = "/"), call. = FALSE)
    }

    target <- if (database == "Hallmark") file.path(base, "Hallmark") else base
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    write.csv(result$full, file.path(target, paste0(database, "_results.csv")), row.names = FALSE)
    write.csv(result$significant,
              file.path(target, paste0(database, "_significant_FDR005.csv")), row.names = FALSE)
    write_ids(result$effective_universe,
              file.path(target, paste0(database, "_actual_universe.csv")))
    saveRDS(result$native, file.path(target, paste0(database, "_enrichment.rds")))

    if (database == "Hallmark") {
      overlaps <- data.frame(ID = character(), gene_id = character())
      if (nrow(result$significant)) {
        overlaps <- do.call(rbind, lapply(seq_len(nrow(result$significant)), function(i) {
          data.frame(ID = result$significant$ID[[i]],
                     gene_id = strsplit(result$significant$geneID[[i]], "/", fixed = TRUE)[[1L]])
        }))
      }
      write.csv(overlaps, file.path(target, "Hallmark_significant_gene_overlaps.csv"),
                row.names = FALSE)
    }

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = dataset, database = database, selected = length(genes),
      background = length(universe), tested = counts[[1L]], significant = counts[[2L]],
      effective_universe = length(result$effective_universe)
    )
    message("PASS ", dataset, " / ", database, ": ", counts[[1L]], " tested, ",
            counts[[2L]], " significant")
  }
}

write.csv(do.call(rbind, summary_rows), file.path(out_dir, "enrichment_summary.csv"),
          row.names = FALSE)
write.csv(data.frame(path = input_files, md5 = unname(tools::md5sum(input_files))),
          file.path(out_dir, "input_md5.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines("PASS", file.path(out_dir, "STATUS.txt"))
message("Enrichment results saved in: ", normalizePath(out_dir))
