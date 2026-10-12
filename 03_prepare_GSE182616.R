# GSE182616: prepare GEO expression, sample metadata and probe annotation.
# Requires the data.table package and the gzip command 

project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table with install.packages('data.table').")
}

out <- Sys.getenv("GSE182616_DATA", "")
if (!nzchar(out)) stop("Set GSE182616_DATA.")
if (file.exists(file.path(out, "expression_as_deposited.rds"))) {
  stop("Prepared outputs already exist: ", out)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)

series <- file.path(out, "GSE182616_series_matrix.txt.gz")
platform <- file.path(out, "GPL17077.soft")

series_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/",
  "GSE182616/matrix/GSE182616_series_matrix.txt.gz"
)
platform_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
  "acc=GPL17077&targ=self&form=text&view=full"
)

# Download the complete GEO matrix, never reusing a truncated gzip file.
gzip_valid <- function(path) {
  file.exists(path) &&
    system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE) == 0L
}

if (!gzip_valid(series)) {
  partial <- paste0(series, ".partial")
  if (file.exists(partial)) file.remove(partial)
  download.file(series_url, partial, method = "libcurl", mode = "wb")
  if (!gzip_valid(partial)) stop("Series Matrix download is incomplete.")
  if (file.exists(series)) file.remove(series)
  stopifnot(file.rename(partial, series))
}

if (!file.exists(platform)) {
  download.file(platform_url, platform, method = "libcurl", mode = "wb")
}

# Read the sample descriptions in the GEO Series Matrix header.
con <- gzfile(series, open = "rt")
header <- character()
repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) {
    close(con)
    stop("The GEO Series Matrix header is incomplete.")
  }
  header <- c(header, line)
  if (line == "!series_matrix_table_begin") break
}
close(con)

fields <- list()
for (line in header[startsWith(header, "!Sample_")]) {
  entry <- read.delim(
    text = line, header = FALSE, quote = '"',
    colClasses = "character", check.names = FALSE
  )
  key <- entry[[1]]
  values <- unlist(entry[-1], use.names = FALSE)
  if (key == "!Sample_characteristics_ch1") {
    names_found <- sub(":.*$", "", values)
    stopifnot(length(unique(names_found)) == 1L)
    key <- names_found[1]
    values <- trimws(sub("^[^:]+:", "", values))
  }
  fields[[key]] <- values
}

to_numeric <- function(values) {
  values[values %in% c("NA", "")] <- NA_character_
  as.numeric(values)
}

metadata <- data.frame(
  gsm = fields[["!Sample_geo_accession"]],
  title = fields[["!Sample_title"]],
  time_label = fields[["time point"]],
  age = to_numeric(fields[["age"]]),
  gender = fields[["gender"]],
  tbsa = to_numeric(fields[["tbsa"]]),
  mortality = fields[["mortality"]],
  stringsAsFactors = FALSE
)

# The title prefix identifies the patient; hr and d denote hours and days.
normalize_hours <- function(labels) {
  labels <- tolower(gsub("[[:space:]_]", "", labels))
  if (!all(grepl("^(hr|d)[0-9]+$", labels))) stop("Unknown time label.")
  hours <- as.numeric(sub("^(hr|d)", "", labels))
  hours[startsWith(labels, "d")] <- hours[startsWith(labels, "d")] * 24
  hours
}

stopifnot(!anyDuplicated(metadata$gsm))
stopifnot(all(grepl("^[0-9]+\\.[0-9]{3}\\.[0-9]{3}_", metadata$title)))
metadata$subject_token <- sub("\\..*$", "", metadata$title)
metadata$subject_id <- sprintf("%04d", as.integer(metadata$subject_token))
metadata$time_hours <- normalize_hours(metadata$time_label)
stopifnot(identical(
  normalize_hours(sub("^[^_]+_", "", metadata$title)),
  metadata$time_hours
))

# Check that repeated visits refer to one clinically consistent patient.
stopifnot(!anyNA(metadata$tbsa), all(metadata$tbsa >= 0 & metadata$tbsa <= 100))
stopifnot(!anyDuplicated(metadata[c("subject_id", "time_hours")]))
for (patient in split(metadata, metadata$subject_id)) {
  for (field in c("tbsa", "age", "gender", "mortality")) {
    values <- unique(patient[[field]])
    values <- values[!is.na(values) & nzchar(as.character(values))]
    if (length(values) > 1L) stop("Conflicting patient metadata: ", patient$subject_id[1])
  }
}
write.csv(metadata, file.path(out, "sample_metadata_audited.csv"),
          row.names = FALSE, na = "NA")

# Read GPL17077 and mark technical control probes.
platform_lines <- readLines(platform, warn = FALSE)
start <- match("!platform_table_begin", platform_lines)
end <- match("!platform_table_end", platform_lines)
if (anyNA(c(start, end))) stop("The GPL17077 file is incomplete.")

annotation <- data.table::fread(
  text = paste(platform_lines[(start + 1L):(end - 1L)], collapse = "\n"),
  sep = "\t", quote = "", colClasses = "character", data.table = FALSE
)
stopifnot(!anyDuplicated(annotation$ID))
controls <- annotation$SPOT_ID == "CONTROL"
stopifnot(identical(controls, annotation$CONTROL_TYPE %in% c("pos", "neg")))
stopifnot(all(annotation$CONTROL_TYPE[!controls] == "FALSE"))
annotation$retained <- !controls
annotation$filter_reason <- ifelse(
  controls, "Documented GPL control probe", "Biological probe; no outcome filtering"
)
write.csv(annotation, file.path(out, "feature_annotation.csv"),
          row.names = FALSE, na = "NA")

# Convert deposited GEO values into a probe-by-sample numeric matrix.
tab <- data.table::fread(
  cmd = paste("gzip -dc", shQuote(series)),
  skip = '"ID_REF"', nrows = nrow(annotation),
  data.table = FALSE, check.names = FALSE
)
stopifnot(identical(names(tab)[-1], metadata$gsm))
stopifnot(!anyDuplicated(tab[[1]]), setequal(tab[[1]], annotation$ID))

expression <- as.matrix(tab[-1])
rownames(expression) <- tab[[1]]
stopifnot(is.numeric(expression), !any(is.infinite(expression)))
saveRDS(expression, file.path(out, "expression_as_deposited.rds"))

# Summarize the observations available for the later TensorPLS analysis.
times <- c(0, 2, 4, 8, 12, 24)
visits <- metadata[metadata$time_hours %in% times, , drop = FALSE]

cat("\nGEO expression dimensions (probes x GSM samples):\n")
print(dim(expression))
cat("\nProbe types (FALSE = control, TRUE = biological):\n")
print(table(annotation$retained))
cat("\nSelected early-window visits:", nrow(visits), "\n")
cat("Patients with early visits:", length(unique(visits$subject_id)), "\n")
cat("\nSaved input files to:", out, "\n")

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
