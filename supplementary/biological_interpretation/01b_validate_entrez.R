#!/usr/bin/env Rscript

# Validate mapped Entrez IDs against the recorded human annotation database.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    "Usage: Rscript 01b_validate_entrez.R MAPPING_OUTPUT_DIR NEW_OUTPUT_DIR [REFERENCE_BIOLOGY_DIR]",
    call. = FALSE
  )
}

mapping_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
requested_output <- path.expand(args[[2L]])
output_parent <- normalizePath(dirname(requested_output), mustWork = TRUE)
output_dir <- file.path(output_parent, basename(requested_output))
reference_dir <- if (length(args) == 3L) {
  normalizePath(path.expand(args[[3L]]), mustWork = TRUE)
} else {
  mapping_dir
}

inside_dir <- function(path, dir) {
  identical(path, dir) || startsWith(path, paste0(dir, "/"))
}
if (inside_dir(output_parent, mapping_dir) || inside_dir(output_parent, reference_dir)) {
  stop("Output must be outside the mapping and reference directories.", call. = FALSE)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output directory already exists: ", output_dir, call. = FALSE)
}

for (dataset in c("GSE68310", "GSE182616")) {
  input_dir <- file.path(mapping_dir, "results", dataset)
  inputs <- file.path(input_dir, c("selected_genes.csv", "background_genes.csv"))
  if (!all(file.exists(inputs))) {
    stop("Missing mapping input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "), call. = FALSE)
  }
  if (length(args) == 3L) {
    refs <- file.path(reference_dir, "results", dataset,
                      c("validated_selected_genes.csv", "validated_background_genes.csv"))
    if (!all(file.exists(refs))) {
      stop("Missing reference file(s): ", paste(refs[!file.exists(refs)], collapse = ", "), call. = FALSE)
    }
  }
}

for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing Bioconductor package: ", pkg, call. = FALSE)
  }
}

expected_orgdb <- "3.23.1"
actual_orgdb <- as.character(utils::packageVersion("org.Hs.eg.db"))
if (actual_orgdb != expected_orgdb) {
  stop(
    "Exact replay requires org.Hs.eg.db ", expected_orgdb,
    "; installed: ", actual_orgdb,
    ". Install the recorded version before comparing results.",
    call. = FALSE
  )
}

known_ids <- unique(AnnotationDbi::keys(
  org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID"
))

read_gene_set <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  dat <- utils::read.csv(
    path, colClasses = "character", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!nrow(dat) || !ncol(dat)) stop("Empty input: ", path, call. = FALSE)

  normalized_names <- gsub("[^a-z0-9]", "", tolower(names(dat)))
  candidates <- which(normalized_names %in% c(
    "entrezid", "entrezgeneid", "geneid", "geneentrezid", "ncbigeneid"
  ))
  if (ncol(dat) == 1L) {
    column <- 1L
  } else if (length(candidates) == 1L) {
    column <- candidates
  } else {
    stop(
      "Cannot identify a unique Entrez ID column in ", path,
      ". Columns: ", paste(names(dat), collapse = ", "),
      call. = FALSE
    )
  }

  ids <- trimws(dat[[column]])
  if (anyNA(ids) || any(!grepl("^[1-9][0-9]*$", ids))) {
    stop("Invalid or missing Entrez IDs in: ", path, call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("Input should already contain unique IDs: ", path, call. = FALSE)
  }
  list(ids = ids, column = names(dat)[[column]])
}

write_gene_set <- function(ids, path, column) {
  output <- data.frame(ids, stringsAsFactors = FALSE)
  names(output) <- column
  utils::write.csv(output, file = path, row.names = FALSE)
}

expectations <- list(
  GSE68310 = c(selected_before = 990L, background_before = 30466L,
               selected_after = 921L, background_after = 23306L),
  GSE182616 = c(selected_before = 1672L, background_before = 21755L,
                selected_after = 1651L, background_after = 20859L)
)

if (!dir.create(output_dir, showWarnings = FALSE)) {
  stop("Cannot create output directory: ", output_dir, call. = FALSE)
}

checks <- list()
for (dataset in names(expectations)) {
  input_dir <- file.path(mapping_dir, "results", dataset)
  historical_selected <- read_gene_set(file.path(input_dir, "selected_genes.csv"))
  historical_background <- read_gene_set(file.path(input_dir, "background_genes.csv"))

  if (!all(historical_selected$ids %in% historical_background$ids)) {
    stop("Selected genes are not a subset of background: ", dataset,
         call. = FALSE)
  }

  selected <- historical_selected$ids[historical_selected$ids %in% known_ids]
  background <- historical_background$ids[historical_background$ids %in% known_ids]
  excluded_selected <- historical_selected$ids[!(historical_selected$ids %in% known_ids)]
  excluded_background <- historical_background$ids[!(historical_background$ids %in% known_ids)]

  actual <- c(
    selected_before = length(historical_selected$ids),
    background_before = length(historical_background$ids),
    selected_after = length(selected),
    background_after = length(background)
  )
  if (!identical(unname(as.integer(actual)), unname(expectations[[dataset]]))) {
    stop(
      dataset, " counts differ from the audited reference. Got ",
      paste(names(actual), actual, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  stopifnot(!anyDuplicated(selected), !anyDuplicated(background),
            all(selected %in% background))

  ref_sel_path <- file.path(reference_dir, "results", dataset, "validated_selected_genes.csv")
  ref_bg_path <- file.path(reference_dir, "results", dataset, "validated_background_genes.csv")
  sel_column <- historical_selected$column
  bg_column <- historical_background$column

  if (file.exists(ref_sel_path)) {
    reference <- read_gene_set(ref_sel_path)
    if (!setequal(selected, reference$ids)) {
      stop("Selected IDs disagree with saved reference: ", dataset,
           call. = FALSE)
    }
    sel_column <- reference$column
  }
  if (file.exists(ref_bg_path)) {
    reference <- read_gene_set(ref_bg_path)
    if (!setequal(background, reference$ids)) {
      stop("Background IDs disagree with saved reference: ", dataset,
           call. = FALSE)
    }
    bg_column <- reference$column
  }

  out <- file.path(output_dir, dataset)
  targets <- file.path(out, c(
    "validated_selected_genes.csv", "validated_background_genes.csv",
    "unrecognised_selected_Entrez.csv", "unrecognised_background_Entrez.csv"
  ))
  if (any(file.exists(targets))) {
    stop("Output files already exist in ", out,
         ". Choose a fresh output directory.", call. = FALSE)
  }
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  write_gene_set(selected, targets[[1L]], sel_column)
  write_gene_set(background, targets[[2L]], bg_column)
  write_gene_set(excluded_selected, targets[[3L]], sel_column)
  write_gene_set(excluded_background, targets[[4L]], bg_column)

  checks[[dataset]] <- data.frame(
    dataset = dataset,
    historical_selected = length(historical_selected$ids),
    historical_background = length(historical_background$ids),
    validated_selected = length(selected),
    validated_background = length(background),
    excluded_selected = length(excluded_selected),
    excluded_background = length(excluded_background),
    orgdb_version = actual_orgdb
  )
  message("PASS ", dataset, ": ", length(selected), " selected / ",
          length(background), " background Entrez IDs")
}

utils::write.csv(
  do.call(rbind, checks), file.path(output_dir, "validation_summary.csv"),
  row.names = FALSE
)
utils::capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Reproduced files saved in: ", normalizePath(output_dir))
