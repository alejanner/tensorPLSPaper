# Rebuild the human Hallmark TERM2GENE table from MSigDB 2026.1.Hs.
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% 1:2)) {
  stop("Usage: Rscript 00_generate_hallmark.R OUTPUT_CSV [FROZEN_CSV]")
}
output <- path.expand(args[1])
if (file.exists(output)) stop("Output already exists: ", output)
if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr.")
stopifnot(as.character(packageVersion("msigdbr")) == "26.1.1")

source <- msigdbr::msigdbr(
  db_species = "HS", species = "Homo sapiens", collection = "H"
)
stopifnot(
  identical(unique(source$db_version), "2026.1.Hs"),
  all(source$gs_collection == "H"),
  !anyNA(source$ncbi_gene),
  identical(as.character(source$ncbi_gene), as.character(source$db_ncbi_gene))
)

result <- unique(data.frame(
  term = as.character(source$gs_name),
  gene = as.character(source$ncbi_gene)
))
stopifnot(
  length(unique(result$term)) == 50L,
  all(grepl("^HALLMARK_", result$term)),
  all(grepl("^[1-9][0-9]*$", result$gene))
)

if (length(args) == 2L) {
  reference <- path.expand(args[2])
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest.")
  if (!file.exists(reference)) stop("Frozen CSV not found: ", reference)
  expected_sha <- "e83e3a01936a4dd9bde6f870dffe8f1f92f679565ec1e7fae4427efeb3f71628"
  stopifnot(identical(digest::digest(file = reference, algo = "sha256"), expected_sha))
  frozen <- read.csv(reference, colClasses = "character", check.names = FALSE)
  stopifnot(identical(names(frozen), c("term", "gene")))
  ids <- function(x) sort(paste(x$term, x$gene, sep = "\t"))
  stopifnot(identical(ids(result), ids(frozen)))
  cat("PASS: all term-gene pairs match the validated CSV.\n")
}

write.csv(result, output, row.names = FALSE)
cat("Gene sets:", length(unique(result$term)), " | Pairs:", nrow(result), "\n")
cat("Saved:", output, "\n")
