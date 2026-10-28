args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)

source_dir <- normalizePath(args[1], mustWork = TRUE)
main_dir <- normalizePath(args[2], mustWork = TRUE)
out_dir <- path.expand(args[3])

if (file.exists(out_dir)) stop("Output directory already exists.")

lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))

pkgload::load_all(source_dir, export_all = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

strategies <- c(V = "VIP2D", M1 = "M1", M2 = "M2")

tucker_control <- list(
  ranks = c(4L, 4L, 3L),
  tolerance = 1e-6,
  max_iter = 2000L,
  ridge = 0.001,
  seed = 123L,
  initialization = "hosvd",
  projection_tolerance = 1e-10
)

for (dataset in c("GSE68310", "GSE182616")) {

  data_dir <- file.path(main_dir,
    if (dataset == "GSE68310") "GSE68310_data" else "GSE182616_analysis"
  )

  analysis_dir <- file.path(main_dir, paste0(dataset, "_analysis"))

  X <- readRDS(file.path(data_dir, "X_incomplete.rds"))
  Y <- readRDS(file.path(data_dir, "Y.rds"))
  folds <- readRDS(file.path(analysis_dir, "nested_cv_folds.rds"))

  for (name in names(strategies)) {

    message("Running ", dataset, " / ", name)

    options <- list(
      X = X,
      Y = Y,
      folds = folds,
      outer_folds = 5L,
      inner_folds = 5L,
      seed = 123L,
      ncomp = 1:10,
      feature_strategy = strategies[[name]],
      vip_thr = 99,
      strategy_mode = "fixed",
      imputation = if (dataset == "GSE68310") "tucker" else "training_mean",
      input_provenance = "unknown",
      selection_metric = "Q2Y",
      keep_models = FALSE,
      verbose = TRUE
    )

    if (dataset == "GSE68310")
      options$tucker_control <- tucker_control

    cv <- do.call(TensorPLS::nested_cv_tensorpls, options)

    stopifnot(
      !length(cv$errors),
      all(cv$outer_fold_metrics$status == "ok"),
      identical(cv$folds, folds)
    )

    destination <- file.path(out_dir, "runs", dataset, name)
    dir.create(destination, recursive = TRUE)

    saveRDS(cv, file.path(destination, "nested_cv_result.rds"))

    message("Saved: ", destination)

    rm(cv)
    gc(verbose = FALSE)
  }

  rm(X, Y, folds)
  gc(verbose = FALSE)
}