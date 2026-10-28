# Held-out visit reconstruction using the paper reproduction tensor and folds.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 01_gse68310_reconstruction.R PRIMARY_RUN OUTPUT [HISTORICAL_FULL] [--preflight]"
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0) }
preflight <- "--preflight" %in% args
args <- args[args != "--preflight"]
if (length(args) < 2L || length(args) > 3L) stop(usage)
primary_dir <- normalizePath(path.expand(args[1]), mustWork = TRUE)
out_dir <- path.expand(args[2])
reference_dir <- if (length(args) == 3L) normalizePath(path.expand(args[3]), mustWork = TRUE) else ""
inputs <- file.path(primary_dir, c("GSE68310_data/X_incomplete.rds", "GSE68310_analysis/nested_cv_folds.rds"))
if (!all(file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
if (file.exists(out_dir)) stop("Output already exists: ", out_dir)
if (nzchar(reference_dir) && !file.exists(file.path(reference_dir, "reconstruction_metrics.csv")))
  stop("Missing historical reconstruction_metrics.csv")
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
old <- if (nzchar(reference_dir)) read.csv(file.path(reference_dir, "reconstruction_metrics.csv")) else NULL

stopifnot(
  identical(dim(X), c(109L, 47162L, 4L)),
  identical(names(folds), dimnames(X)[[1]]),
  setequal(folds, 1:5),
  sum(is.na(X)) == 990402L
)
if (!is.null(old)) stopifnot(nrow(old) == 150L)
if (preflight) { cat("Preflight passed; no fitting or output creation.\n"); quit(status = 0) }
dir.create(out_dir, recursive = TRUE)
cat("New output:", out_dir, "\n")
saveRDS(list(primary_dir = primary_dir, input_md5 = tools::md5sum(inputs)),
        file.path(out_dir, "input_provenance.rds"))
saveRDS(folds, file.path(out_dir, "outer_folds.rds"))
dir.create(file.path(out_dir, "masks"))

levels <- c(.05, .10, .20, .30, .40)
seeds <- c(123L, 124L, 125L)
rng <- c("Mersenne-Twister", "Inversion", "Rejection")
tucker_control <- list(
  ranks = c(4L, 4L, 3L), ridge = .001, tolerance = 1e-6,
  max_iter = 2000L, initialization = "hosvd",
  projection_tolerance = 1e-10, seed = 123L
)

# Generate visit masks while retaining at least one observed visit.
make_mask <- function(observed_visits, level, seed) {
  do.call(RNGkind, as.list(rng))
  set.seed(seed)
  candidates <- which(observed_visits, arr.ind = TRUE)
  permutation <- sample.int(nrow(candidates))
  remaining <- rowSums(observed_visits)
  target <- as.integer(floor(level * nrow(candidates) + .5))
  hidden <- matrix(FALSE, nrow(observed_visits), ncol(observed_visits),
                   dimnames = dimnames(observed_visits))
  for (index in permutation) {
    if (sum(hidden) == target) break
    subject <- candidates[index, 1]
    time <- candidates[index, 2]
    if (remaining[subject] <= 1L) next
    hidden[subject, time] <- TRUE
    remaining[subject] <- remaining[subject] - 1L
  }
  stopifnot(sum(hidden) == target, all(remaining >= 1L))
  list(hidden = hidden, target = target, achieved = sum(hidden),
       seed = seed, rng = rng, permutation = permutation)
}

results <- list()
fit_summary <- list()
mean_fit_fn <- getFromNamespace(".nested_imputer", "TensorPLS")
mean_apply <- getFromNamespace(".nested_impute", "TensorPLS")

for (fold in 1:5) {
  train <- X[folds != fold, , , drop = FALSE]
  test <- X[folds == fold, , , drop = FALSE]
  cat("\nFold", fold, ": fitting Tucker and training mean from scratch\n")

  mean_fit <- mean_fit_fn(matrix(train, dim(train)[1]), "training_mean")
  tucker_fit <- do.call(TensorPLS::fit_tucker_imputer,
                        c(list(X_train = train), tucker_control))
  stopifnot(isTRUE(tucker_fit$convergence$converged),
            identical(tucker_fit$training_ids, dimnames(train)[[1]]),
            identical(tucker_fit$X_train_imputed[!is.na(train)], train[!is.na(train)]))
  fit_summary[[fold]] <- data.frame(
    fold = fold, iterations = tucker_fit$convergence$iterations,
    converged = tucker_fit$convergence$converged,
    training_RMSE = tucker_fit$convergence$observed_RMSE
  )

  # The completed training tensor is not needed for applying this frozen fit.
  tucker_fit$X_train_imputed <- NULL
  saveRDS(list(mean = mean_fit, tucker = tucker_fit),
          file.path(out_dir, sprintf("fit_fold_%d.rds", fold)))

  observed_visits <- apply(!is.na(test), c(1, 3), any)
  for (seed in seeds) for (level in levels) {
    mask <- make_mask(observed_visits, level, seed)
    tag <- sprintf("fold_%d_seed_%d_mask_%02d", fold, seed, round(100 * level))
    if (nzchar(reference_dir)) {
      saved_mask <- readRDS(file.path(reference_dir, paste0(tag, "_mask.rds")))
      stopifnot(identical(mask, saved_mask$mask),
                identical(observed_visits, saved_mask$observed_visits))
    }
    saveRDS(list(mask = mask, observed_visits = observed_visits,
                 subject_ids = dimnames(test)[[1]], time_names = dimnames(test)[[3]]),
            file.path(out_dir, "masks", paste0(tag, "_mask.rds")))

    evaluation <- array(FALSE, dim(test), dimnames(test))
    for (ti in seq_len(dim(test)[3])) {
      evaluation[mask$hidden[, ti], , ti] <- TRUE
    }
    stopifnot(!any(evaluation & is.na(test)),
              sum(evaluation) == mask$achieved * dim(test)[2])
    masked <- test
    masked[evaluation] <- NA_real_
    truth <- test[evaluation]

    mean_completed <- array(
      mean_apply(mean_fit, matrix(masked, dim(masked)[1])),
      dim(masked), dimnames(masked)
    )
    projection <- TensorPLS::predict_tucker_imputer(
      tucker_fit, masked, diagnostics = TRUE
    )
    stopifnot(!any(projection$diagnostics$rank_deficient))

    predictions <- list(training_mean = mean_completed, tucker = projection$X)
    for (method in names(predictions)) {
      completed <- predictions[[method]]
      stopifnot(identical(completed[!is.na(masked)], masked[!is.na(masked)]),
                all(is.finite(completed[evaluation])))
      estimate <- completed[evaluation]
      error <- estimate - truth
      results[[length(results) + 1L]] <- data.frame(
        fold = fold, seed = seed, level = level, method = method,
        n = length(error), squared_error = sum(error^2),
        RMSE = sqrt(mean(error^2)), Pearson = cor(truth, estimate)
      )
    }
    cat(sprintf("  seed %d, masked %d%% of observed visits\n",
                seed, round(level * 100)))
  }
  rm(train, test, mean_fit, tucker_fit)
  gc(verbose = FALSE)
}

results <- do.call(rbind, results)
fit_summary <- do.call(rbind, fit_summary)
write.csv(results, file.path(out_dir, "reconstruction_metrics.csv"), row.names = FALSE)
write.csv(fit_summary, file.path(out_dir, "fit_summary.csv"), row.names = FALSE)

pooled <- do.call(rbind, lapply(
  split(results, interaction(results$method, results$level, drop = TRUE)),
  function(z) data.frame(method = z$method[1], level = z$level[1],
                         cases = nrow(z),
                         pooled_RMSE = sqrt(sum(z$squared_error) / sum(z$n)))
))
write.csv(pooled, file.path(out_dir, "pooled_RMSE.csv"), row.names = FALSE)

# Optional verification against historical case-level results.
case_id <- function(z) paste(z$fold, z$seed, z$level, z$method)
stopifnot(nrow(results) == 150L, !anyDuplicated(case_id(results)))
if (!is.null(old)) {
  stopifnot(!anyDuplicated(case_id(old)),
            setequal(case_id(results), case_id(old)), all(old$status == "ok"))
  old <- old[match(case_id(results), case_id(old)), ]
  comparison <- data.frame(
    fold = results$fold, seed = results$seed, level = results$level,
    method = results$method, RMSE_old = old$RMSE,
    RMSE_new = results$RMSE,
    absolute_RMSE_difference = abs(results$RMSE - old$RMSE)
  )
  write.csv(comparison, file.path(out_dir, "RMSE_comparison.csv"), row.names = FALSE)

  stopifnot(identical(as.integer(results$n), as.integer(old$n)),
            all(abs(results$squared_error - old$squared_error) < 1e-8),
            all(abs(results$Pearson - old$Pearson) < 1e-10),
            all(comparison$absolute_RMSE_difference < 1e-10))
}
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("\nPASS: 5 Tucker fits, 75 masks, 150 evaluations.\n")
print(pooled, row.names = FALSE)
writeLines(if (is.null(old)) "PASS: 5 Tucker fits; 75 masks; 150 evaluations; reference comparison not requested." else
  "PASS: 5 Tucker fits; 75 masks; 150 evaluations; reference comparison passed.", file.path(out_dir, "PASS.txt"))
