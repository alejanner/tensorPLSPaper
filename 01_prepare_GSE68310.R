# GSE68310: prepare deposited expression data and the tensor.

# Paths and dependencies
data_dir <- Sys.getenv("GSE68310_DATA", "")
source_dir <- Sys.getenv("TENSORPLS_SOURCE_DIR", "")
if (!nzchar(data_dir) || !nzchar(source_dir)) stop("Set GSE68310_DATA and TENSORPLS_SOURCE_DIR.")
source_dir <- normalizePath(source_dir, mustWork = TRUE)
project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
if ("TensorPLS" %in% loadedNamespaces()) stop("Start a fresh R session.")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The pkgload package is required.")
pkgload::load_all(source_dir, export_all = FALSE)
stopifnot(identical(normalizePath(getNamespaceInfo("TensorPLS", "path")), source_dir))

download_dir <- path.expand(Sys.getenv(
  "GSE68310_DOWNLOADS", file.path(data_dir, "downloads")
))

required <- c("data.table", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the required R packages: ", paste(missing, collapse = ", "))
}

prepared <- file.path(data_dir, c("X_incomplete.rds", "Y.rds", "subjects.csv", "provenance.rds"))
if (any(file.exists(prepared))) {
  stop("Prepared outputs already exist in ", data_dir, ". Choose a new GSE68310_DATA directory.")
}

dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(1800, getOption("timeout")))

# Download and verify GEO files 
base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68310/"
sources <- data.frame(
  file = c(
    "GSE68310_series_matrix.txt.gz",
    "GSE68310_SubjectPhenotypes1.txt.gz",
    "GSE68310_RAW.tar"
  ),
  url = paste0(base_url, c("matrix/", "suppl/", "suppl/"), c(
    "GSE68310_series_matrix.txt.gz",
    "GSE68310_SubjectPhenotypes1.txt.gz",
    "GSE68310_RAW.tar"
  )),
  SHA256 = c(
    "3ba7898be1c6c279fb2779b1e7bf99b887374fccf4221a0c7bcb7b2cfae312de",
    "9285b14adef0338f41334bb1a75e43d7b4ffa58130eb63583cd4361ee64c6cbe",
    "19dbc61b4da9ac807799a96081f4eae405079dad108b82a9d2bdb74678cdb948"
  )
)

verified_download <- function(url, destination, expected_sha256) {
  if (file.exists(destination)) {
    if (digest::digest(file = destination, algo = "sha256") != expected_sha256) {
      stop("SHA256 mismatch for existing file: ", destination)
    }
    return(invisible(destination))
  }
  
  partial <- paste0(destination, ".partial")
  if (file.exists(partial)) unlink(partial)
  on.exit(unlink(partial), add = TRUE)
  
  status <- utils::download.file(url, partial, method = "libcurl", mode = "wb")
  if (status != 0L || digest::digest(file = partial, algo = "sha256") != expected_sha256) {
    stop("Download failed or SHA256 mismatch: ", url)
  }
  if (!file.rename(partial, destination)) stop("Cannot save: ", destination)
  invisible(destination)
}

source_paths <- file.path(download_dir, sources$file)
for (i in seq_len(2L)) {
  verified_download(sources$url[i], source_paths[i], sources$SHA256[i])
}

series_file <- source_paths[1L]
phenotype_file <- source_paths[2L]
raw_archive <- source_paths[3L]

# Read and verify GEO sample metadata 
connection <- gzfile(series_file, "rt")
header <- character()
repeat {
  line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(line)) {
    close(connection)
    stop("The GEO Series Matrix header is incomplete.")
  }
  header <- c(header, line)
  if (line == "!series_matrix_table_begin") break
}
close(connection)

fields <- list()
for (line in header[startsWith(header, "!Sample_")]) {
  entry <- read.delim(
    text = line, header = FALSE, quote = '"',
    colClasses = "character", check.names = FALSE
  )
  key <- entry[[1L]]
  values <- unlist(entry[-1L], use.names = FALSE)
  if (key == "!Sample_characteristics_ch1") {
    characteristic <- sub(":.*$", "", values)
    stopifnot(length(unique(characteristic)) == 1L)
    key <- characteristic[1L]
    values <- sub("^[^:]+: *", "", values)
  }
  fields[[key]] <- values
}

samples <- data.frame(
  sample = fields[["!Sample_geo_accession"]],
  title = fields[["!Sample_title"]],
  platform = fields[["!Sample_platform_id"]],
  subject_id = fields[["subject id"]],
  time_point = fields[["time point"]],
  gender = fields[["gender"]],
  infection = fields[["infection"]]
)

subjects <- data.table::fread(
  cmd = paste("gzip -dc", shQuote(phenotype_file)), data.table = FALSE
)
stopifnot(
  nrow(samples) == 880L,
  nrow(subjects) == 133L,
  !anyDuplicated(samples$sample),
  !anyDuplicated(subjects$SUBJECT_ID),
  setequal(samples$subject_id, subjects$SUBJECT_ID),
  all(samples$platform == "GPL10558"),
  all(samples$title == paste("WholeBloodRNA", samples$subject_id, samples$time_point, sep = "_")),
  !anyDuplicated(samples[c("subject_id", "time_point")])
)

subject_index <- match(samples$subject_id, subjects$SUBJECT_ID)
expected_infection <- tolower(gsub("_", " ", subjects$Virus_by_PCR[subject_index]))
expected_infection[expected_infection == "none detected"] <-
  "our tests did not detect one of the viruses sought"
stopifnot(
  all(expected_infection == tolower(samples$infection)),
  all(c(F = "female", M = "male")[subjects$Gender[subject_index]] == samples$gender)
)

# Select PCR-positive subjects and four acute visits 
subjects$influenza <- grepl("influenza_", subjects$Virus_by_PCR)
subjects$coinfection <- grepl("_and_", subjects$Virus_by_PCR)
subjects$primary_class <- ifelse(
  subjects$influenza, "Influenza",
  ifelse(subjects$Virus_by_PCR == "none_detected", "Excluded_negative", "Other_virus")
)
write.csv(subjects, file.path(data_dir, "all_subjects.csv"), row.names = FALSE)

subjects <- subjects[subjects$primary_class != "Excluded_negative", ]
subjects <- subjects[order(subjects$SUBJECT_ID), ]
ids <- subjects$SUBJECT_ID
times <- c("Day0", "Day2", "Day4", "Day6")
samples <- samples[
  samples$subject_id %in% ids & samples$time_point %in% times, , drop = FALSE
]
stopifnot(
  length(ids) == 109L,
  sum(subjects$primary_class == "Influenza") == 73L,
  sum(subjects$primary_class == "Other_virus") == 36L,
  nrow(samples) == 415L
)

# Identify ERCC controls from the Illumina manifest 
manifest_name <- "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"
manifest_file <- file.path(download_dir, manifest_name)
if (!file.exists(manifest_file)) {
  verified_download(sources$url[3L], raw_archive, sources$SHA256[3L])
  utils::untar(raw_archive, files = manifest_name, exdir = download_dir)
}
if (file.exists(raw_archive)) {
  stopifnot(digest::digest(file = raw_archive, algo = "sha256") == sources$SHA256[3L])
}
stopifnot(
  file.exists(manifest_file),
  digest::digest(file = manifest_file, algo = "sha256") ==
    "932e2f110bd30ac70ef4d7dbac2f87f8088f982bf9aa596732c7477ab50d43c2"
)
sources$verified_locally <- file.exists(source_paths)
write.csv(sources, file.path(data_dir, "source_manifest.csv"), row.names = FALSE)

connection <- gzfile(manifest_file, "rt")
manifest_lines <- readLines(connection, warn = FALSE)
close(connection)
first <- match("[Controls]", manifest_lines)
last <- match("[Columns]", manifest_lines)
stopifnot(!anyNA(c(first, last)), first < last)

controls <- data.table::fread(
  text = paste(manifest_lines[(first + 1L):(last - 1L)], collapse = "\n"),
  sep = "\t", quote = "", colClasses = "character", data.table = FALSE
)
ercc <- unique(controls$Probe_Id[grepl("^ERCC-", controls$Reporter_Group_Name)])

#Read deposited expression values
expression <- data.table::fread(
  cmd = paste("gzip -dc", shQuote(series_file)),
  skip = length(header), nrows = 47254L,
  data.table = FALSE, check.names = FALSE
)
stopifnot(
  nrow(expression) == 47254L,
  ncol(expression) == 881L,
  !anyDuplicated(expression[[1L]]),
  identical(names(expression)[-1L], fields[["!Sample_geo_accession"]])
)

ercc <- ercc[ercc %in% expression[[1L]]]
stopifnot(length(ercc) == 92L)
expression <- expression[!expression[[1L]] %in% ercc, , drop = FALSE]
probes <- expression[[1L]]
values <- as.matrix(expression[samples$sample])
stopifnot(length(probes) == 47162L, is.numeric(values), all(is.finite(values)))

input <- as.data.frame(t(values), check.names = FALSE)
names(input) <- probes
input$subject_id <- samples$subject_id
input$time_point <- samples$time_point

#Construct the incomplete subjects x probes x time tensor
X <- TensorPLS::prepare_omics(
  input, id_col = "subject_id", time_col = "time_point",
  transpose = "never", feature_cols = probes,
  subjects = ids, time_points = times,
  coercion_mode = "custom", id_as = "character", time_as = "character"
)
stopifnot(
  identical(dim(X), c(109L, 47162L, 4L)),
  identical(dimnames(X)[[1L]], ids),
  identical(dimnames(X)[[2L]], probes),
  identical(dimnames(X)[[3L]], times),
  sum(is.na(X)) == 990402L
)

for (i in seq_len(nrow(samples))) {
  stopifnot(identical(
    as.numeric(X[samples$subject_id[i], , samples$time_point[i]]),
    as.numeric(values[, i])
  ))
}

visit_counts <- apply(!is.na(X), c(1L, 3L), sum)
stopifnot(
  all(visit_counts %in% c(0L, 47162L)),
  sum(visit_counts == 0L) == 21L
)

Y <- matrix(
  as.integer(subjects$primary_class == "Influenza"),
  ncol = 1L, dimnames = list(ids, "Influenza")
)
stopifnot(identical(rownames(Y), dimnames(X)[[1L]]))

#Save prepared data and provenance
saveRDS(X, file.path(data_dir, "X_incomplete.rds"))
saveRDS(Y, file.path(data_dir, "Y.rds"))
write.csv(subjects, file.path(data_dir, "subjects.csv"), row.names = FALSE)
write.csv(samples, file.path(data_dir, "sample_mapping.csv"), row.names = FALSE)
write.csv(data.frame(probe = ercc), file.path(data_dir, "excluded_ERCC.csv"), row.names = FALSE)

visit_availability <- data.frame(subject_id = ids, visit_counts > 0L, check.names = FALSE)
names(visit_availability) <- c("subject_id", times)
write.csv(visit_availability, file.path(data_dir, "visit_availability.csv"), row.names = FALSE)

saveRDS(list(
  accession = "GSE68310",
  expression = "Deposited normalized series matrix, no additional normalization",
  outcome = "Influenza-positive versus other PCR-confirmed viruses, coinfections retained",
  exclusion = "none_detected subjects and 92 ERCC probes only",
  times = times,
  dimensions = dim(X),
  missing_cells = sum(is.na(X)),
  input_provenance = "unknown",
  sources = sources,
  manifest = list(file = manifest_name, SHA256 = digest::digest(file = manifest_file, algo = "sha256"))
), file.path(data_dir, "provenance.rds"))

capture.output(sessionInfo(), file = file.path(data_dir, "sessionInfo.txt"))
cat("Prepared:", paste(dim(X), collapse = " x "), "\n")
cat("Missing cells:", sum(is.na(X)), "\n")
cat("Saved to:", data_dir, "\n")

# Check that the newly prepared tensor matches the original input.
reference_dir <- Sys.getenv("GSE68310_REFERENCE_DATA", "")
if (nzchar(reference_dir) && file.exists(file.path(reference_dir, "X_incomplete.rds"))) {
  reference_X <- readRDS(file.path(reference_dir, "X_incomplete.rds"))
  stopifnot(identical(X, reference_X))
  rm(reference_X)
  stopifnot(identical(Y, readRDS(file.path(reference_dir, "Y.rds"))))
  cat("GSE68310 input matches the previously saved X and Y exactly.\n")
}
