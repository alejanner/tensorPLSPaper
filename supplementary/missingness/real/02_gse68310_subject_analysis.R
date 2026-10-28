# Subject-level diagnostics from saved reconstruction imputers and masks.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 02_gse68310_subject_analysis.R RECONSTRUCTION_OUTPUT [PRIMARY_RUN] [--preflight]"
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0) }
preflight <- "--preflight" %in% args
args <- args[args != "--preflight"]
if (length(args) < 1L || length(args) > 2L) stop(usage)
run_dir <- normalizePath(path.expand(args[1]), mustWork = TRUE)
required <- file.path(run_dir, c("PASS.txt", "input_provenance.rds", "outer_folds.rds",
  "reconstruction_metrics.csv", "pooled_RMSE.csv", sprintf("fit_fold_%d.rds", 1:5)))
if (!all(file.exists(required))) stop("Missing benchmark input: ", paste(required[!file.exists(required)], collapse = ", "))
if (file.exists(file.path(run_dir, "subject_level_errors.csv"))) stop("Subject-level results already exist.")
provenance <- readRDS(file.path(run_dir, "input_provenance.rds"))
primary_dir <- if (length(args) == 2L) normalizePath(path.expand(args[2]), mustWork = TRUE) else provenance$primary_dir
inputs <- file.path(primary_dir, c("GSE68310_data/X_incomplete.rds", "GSE68310_analysis/nested_cv_folds.rds"))
if (!all(file.exists(inputs))) stop("Missing primary input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
stopifnot(identical(unname(tools::md5sum(inputs)), unname(provenance$input_md5)))
source_dir <- Sys.getenv("TENSORPLS_SOURCE_DIR", "")
project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (nzchar(project_lib)) .libPaths(c(path.expand(project_lib), .libPaths()))
if (nzchar(source_dir)) {
  source_dir <- normalizePath(path.expand(source_dir), mustWork = TRUE)
  pkgload::load_all(source_dir, export_all = FALSE, quiet = TRUE)
  stopifnot(identical(normalizePath(getNamespaceInfo("TensorPLS", "path")), source_dir))
} else {
  library(TensorPLS)
}
for (name in c("fit_tucker_imputer", "predict_tucker_imputer", ".nested_imputer", ".nested_impute"))
  if (!exists(name, envir = asNamespace("TensorPLS"), inherits = FALSE))
    stop("The loaded TensorPLS version does not provide ", name)
cat("TensorPLS:", getNamespaceInfo("TensorPLS", "path"),
    "version", as.character(packageVersion("TensorPLS")), "\n")
X <- readRDS(file.path(primary_dir, "GSE68310_data/X_incomplete.rds"))
folds <- readRDS(file.path(primary_dir, "GSE68310_analysis/nested_cv_folds.rds"))$outer
cases <- read.csv(file.path(run_dir, "reconstruction_metrics.csv"))
stopifnot(identical(names(folds), dimnames(X)[[1]]), nrow(cases) == 150L)
stopifnot(identical(folds, readRDS(file.path(run_dir, "outer_folds.rds"))))
expected_masks <- expand.grid(fold = 1:5, seed = c(123L, 124L, 125L), percent = c(5, 10, 20, 30, 40))
mask_paths <- file.path(run_dir, "masks", with(expected_masks,
  sprintf("fold_%d_seed_%d_mask_%02d_mask.rds", fold, seed, percent)))
if (!all(file.exists(mask_paths))) stop("Missing saved masks: ", paste(mask_paths[!file.exists(mask_paths)], collapse = ", "))
if (preflight) { cat("Preflight passed; no imputation or output creation.\n"); quit(status = 0) }
mean_apply <- getFromNamespace(".nested_impute", "TensorPLS")

rows <- list()
for (fold in 1:5) {
  test <- X[folds == fold, , , drop = FALSE]
  fit <- readRDS(file.path(run_dir, sprintf("fit_fold_%d.rds", fold)))
  stopifnot(isTRUE(fit$tucker$convergence$converged),
            identical(fit$tucker$training_ids, dimnames(X)[[1]][folds != fold]))
  observed_visits <- apply(!is.na(test), c(1, 3), any)

  for (seed in c(123L, 124L, 125L)) for (level in c(.05, .10, .20, .30, .40)) {
    tag <- sprintf("fold_%d_seed_%d_mask_%02d", fold, seed, round(level * 100))
    saved <- readRDS(file.path(run_dir, "masks", paste0(tag, "_mask.rds")))
    stopifnot(identical(saved$observed_visits, observed_visits),
              identical(saved$subject_ids, dimnames(test)[[1]]))
    hidden <- saved$mask$hidden
    evaluation <- array(FALSE, dim(test), dimnames(test))
    for (ti in seq_len(dim(test)[3])) evaluation[hidden[, ti], , ti] <- TRUE
    stopifnot(!any(evaluation & is.na(test)))

    masked <- test
    masked[evaluation] <- NA_real_
    mean_completed <- array(
      mean_apply(fit$mean, matrix(masked, dim(masked)[1])),
      dim(masked), dimnames(masked)
    )
    tucker_completed <- TensorPLS::predict_tucker_imputer(fit$tucker, masked)
    stopifnot(identical(mean_completed[!is.na(masked)], masked[!is.na(masked)]),
              identical(tucker_completed[!is.na(masked)], masked[!is.na(masked)]),
              all(is.finite(mean_completed[evaluation])),
              all(is.finite(tucker_completed[evaluation])))

    # Reproduce both case-level errors before calculating per-subject errors.
    for (method in c("training_mean", "tucker")) {
      completed <- if (method == "tucker") tucker_completed else mean_completed
      error <- completed[evaluation] - test[evaluation]
      saved_case <- cases[cases$fold == fold & cases$seed == seed &
                            abs(cases$level - level) < 1e-12 & cases$method == method, ]
      stopifnot(nrow(saved_case) == 1L,
                length(error) == saved_case$n,
                abs(sum(error^2) - saved_case$squared_error) < 1e-8,
                abs(sqrt(mean(error^2)) - saved_case$RMSE) < 1e-10)
    }

    remaining <- rowSums(observed_visits & !hidden)
    for (i in which(rowSums(hidden) > 0)) {
      target <- evaluation[i, , ]
      truth <- test[i, , ][target]
      e_mean <- mean_completed[i, , ][target] - truth
      e_tucker <- tucker_completed[i, , ][target] - truth
      rows[[length(rows) + 1L]] <- data.frame(
        subject_id = dimnames(test)[[1]][i], fold = fold, seed = seed,
        level = level, hidden_visits = sum(hidden[i, ]),
        remaining_visits = remaining[i], n = length(truth),
        SSE_mean = sum(e_mean^2), SSE_tucker = sum(e_tucker^2),
        RMSE_mean = sqrt(mean(e_mean^2)),
        RMSE_tucker = sqrt(mean(e_tucker^2))
      )
    }
  }
  cat("Fold", fold, "verified\n")
}

subjects <- do.call(rbind, rows)
subjects$delta_RMSE <- subjects$RMSE_tucker - subjects$RMSE_mean
subjects$tucker_better <- subjects$delta_RMSE < 0
stopifnot(nrow(subjects) == 896L, length(unique(subjects$subject_id)) == 108L)

summarise <- function(z) data.frame(
  evaluations = nrow(z), distinct_subjects = length(unique(z$subject_id)),
  tucker_wins = sum(z$tucker_better),
  tucker_win_percent = 100 * mean(z$tucker_better),
  median_RMSE_mean = median(z$RMSE_mean),
  median_RMSE_tucker = median(z$RMSE_tucker),
  median_delta_RMSE = median(z$delta_RMSE),
  p90_RMSE_tucker = unname(quantile(z$RMSE_tucker, .90)),
  p95_RMSE_tucker = unname(quantile(z$RMSE_tucker, .95)),
  pooled_RMSE_mean = sqrt(sum(z$SSE_mean) / sum(z$n)),
  pooled_RMSE_tucker = sqrt(sum(z$SSE_tucker) / sum(z$n))
)
overall <- data.frame(summarise(subjects))
top_n <- ceiling(.05 * nrow(subjects))
overall$top_5_percent_tucker_SSE_share <- 100 * sum(
  head(sort(subjects$SSE_tucker, decreasing = TRUE), top_n)
) / sum(subjects$SSE_tucker)
by_level <- do.call(rbind, lapply(split(subjects, subjects$level), function(z)
  data.frame(level = z$level[1], summarise(z))))
by_remaining <- do.call(rbind, lapply(split(subjects, subjects$remaining_visits), function(z)
  data.frame(remaining_visits = z$remaining_visits[1], summarise(z))))

# The cell-pooled errors must equal the already reproduced benchmark results.
pooled <- read.csv(file.path(run_dir, "pooled_RMSE.csv"))
for (i in seq_len(nrow(by_level))) {
  z <- by_level[i, ]
  for (method in c("training_mean", "tucker")) {
    expected <- pooled$pooled_RMSE[pooled$level == z$level & pooled$method == method]
    actual <- if (method == "tucker") z$pooled_RMSE_tucker else z$pooled_RMSE_mean
    stopifnot(length(expected) == 1L, abs(expected - actual) < 1e-10)
  }
}

write.csv(subjects, file.path(run_dir, "subject_level_errors.csv"), row.names = FALSE)
write.csv(overall, file.path(run_dir, "subject_summary_overall.csv"), row.names = FALSE)
write.csv(by_level, file.path(run_dir, "subject_summary_by_level.csv"), row.names = FALSE)
write.csv(by_remaining, file.path(run_dir, "subject_summary_by_remaining_visits.csv"), row.names = FALSE)

cat("\nVerified:", nrow(subjects), "subject-mask evaluations across",
    length(unique(subjects$subject_id)), "distinct subjects\n")
cat("Tucker better:", sum(subjects$tucker_better), "/", nrow(subjects), "\n")
cat("Top 5% of subject-mask evaluations: ",
    round(overall$top_5_percent_tucker_SSE_share, 2),
    "% of Tucker SSE\n\n", sep = "")
print(by_level, row.names = FALSE, digits = 4)
cat("\nBy remaining visits:\n")
print(by_remaining, row.names = FALSE, digits = 4)
cat("\nResults saved to:", run_dir, "\n")
