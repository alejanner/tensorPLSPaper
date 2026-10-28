# GSE68310: primary nested CV and full-cohort descriptive fit.
source_dir <- normalizePath(Sys.getenv("TENSORPLS_SOURCE_DIR"), mustWork = TRUE)
data_dir <- normalizePath(Sys.getenv("GSE68310_DATA"), mustWork = TRUE)
out_dir <- Sys.getenv("GSE68310_RESULTS")
if (!nzchar(out_dir) || dir.exists(out_dir)) stop("Choose a new GSE68310_RESULTS directory.")

project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The pkgload package is required.")
pkgload::load_all(source_dir, export_all = FALSE)
stopifnot(identical(normalizePath(getNamespaceInfo("TensorPLS", "path")), source_dir))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

X <- readRDS(file.path(data_dir, "X_incomplete.rds"))
Y <- readRDS(file.path(data_dir, "Y.rds"))
ids <- dimnames(X)[[1]]
stopifnot(identical(dim(X), c(109L, 47162L, 4L)),
          identical(dimnames(X)[[3]], c("Day0", "Day2", "Day4", "Day6")),
          identical(rownames(Y), ids), sum(is.na(X)) == 990402L,
          sum(Y == 1L) == 73L, sum(Y == 0L) == 36L)

seed <- 123L
components <- 1:10
strategy <- "VIP2D_intersect_M2"
threshold <- 99
tucker_control <- list(
  ranks = c(4L, 4L, 3L), tolerance = 1e-6, max_iter = 2000L,
  ridge = 0.001, seed = seed, initialization = "hosvd",
  projection_tolerance = 1e-10
)

# Use the original five-by-five folds, checking that they can be regenerated.
generated <- TensorPLS::make_tensorpls_folds(
  Y, outer_folds = 5L, inner_folds = 5L, seed = seed, subject_ids = ids
)
reference_folds <- Sys.getenv("GSE68310_REFERENCE_FOLDS", "")
folds <- generated
if (nzchar(reference_folds) && file.exists(reference_folds)) {
  folds <- readRDS(reference_folds)
  stopifnot(identical(folds, generated))
  cat("GSE68310: regenerated folds match the original folds.\n")
} else {
  warning("Original GSE68310 folds unavailable; using generated folds.")
}

dir.create(out_dir, recursive = TRUE)
saveRDS(folds, file.path(out_dir, "nested_cv_folds.rds"))

# Predictive evaluation: preprocessing and selection remain inside the training folds.
cv <- TensorPLS::nested_cv_tensorpls(
  X = X, Y = Y, folds = folds, outer_folds = 5L, inner_folds = 5L,
  seed = seed, ncomp = components, feature_strategy = strategy,
  vip_thr = threshold, strategy_mode = "fixed", imputation = "tucker",
  tucker_control = tucker_control, input_provenance = "unknown",
  selection_metric = "Q2Y", keep_models = FALSE, verbose = TRUE
)
saveRDS(cv, file.path(out_dir, "nested_cv_result.rds"))
stopifnot(!length(cv$errors), all(cv$outer_fold_metrics$status == "ok"),
          nrow(cv$outer_predictions) == length(ids),
          !anyDuplicated(cv$outer_predictions$subject_id),
          all(is.finite(cv$outer_predictions$Y_pred)))

metrics <- c(
  Q2Y = cv$Q2Y,
  AUC = cv$classification_metrics[["AUC"]],
  balanced_accuracy = cv$classification_metrics[["balanced_accuracy"]],
  sensitivity = cv$classification_metrics[["sensitivity"]],
  specificity = cv$classification_metrics[["specificity"]]
)
write.csv(data.frame(metric = names(metrics), value = unname(metrics)),
          file.path(out_dir, "nested_cv_metrics.csv"), row.names = FALSE)
write.csv(cv$outer_fold_metrics, file.path(out_dir, "outer_fold_metrics.csv"), row.names = FALSE)
pred <- cv$outer_predictions
write.csv(pred[, !vapply(pred, is.list, logical(1)), drop = FALSE],
          file.path(out_dir, "outer_predictions.csv"), row.names = FALSE)

# Compare against the original subject-level predictions when available.
reference_cv <- Sys.getenv("GSE68310_REFERENCE_CV", "")
if (nzchar(reference_cv) && file.exists(reference_cv)) {
  old <- readRDS(reference_cv)
  stopifnot(identical(old$folds, folds), setequal(old$outer_predictions$subject_id, pred$subject_id))
  old_pred <- old$outer_predictions$Y_pred[match(pred$subject_id, old$outer_predictions$subject_id)]
  max_difference <- max(abs(pred$Y_pred - old_pred))
  write.csv(data.frame(max_absolute_OOF_score_difference = max_difference),
            file.path(out_dir, "reference_prediction_comparison.csv"), row.names = FALSE)
  cat("GSE68310: max absolute difference from saved OOF scores:", max_difference, "\n")
  rm(old)
}

# Full-cohort model: descriptive interpretation, not predictive validation.
development <- TensorPLS::tune_tensorpls_components(
  X = X, Y = Y, ncomp = components, feature_strategy = strategy,
  vip_thr = threshold,
  cv_folds = 5L, folds = folds$outer, seed = seed,
  imputation = "tucker", tucker_control = tucker_control,
  input_provenance = "unknown", selection_metric = "Q2Y", verbose = TRUE
)
H <- development$selected_ncomp
saveRDS(development, file.path(out_dir, "component_tuning.rds"))

fit <- do.call(TensorPLS::fit_tucker_imputer,
               c(list(X_train = X), tucker_control))
X_complete <- fit$X_train_imputed
stopifnot(all(is.finite(X_complete)),
          identical(X_complete[!is.na(X)], X[!is.na(X)]))

model_all <- TensorPLS::fit_tensorpls_pls(X_complete, Y, ncomp = H)
vips_all <- TensorPLS::tensorpls_vips(model_all)
selection <- TensorPLS::vip_feature_sets(vips_all, strategy = strategy, thr = threshold)
X_selected <- X_complete[, selection$selected_indices, , drop = FALSE]
model_selected <- TensorPLS::fit_tensorpls_pls(X_selected, Y, ncomp = H)
vips_selected <- TensorPLS::tensorpls_vips(model_selected)

saveRDS(fit$convergence, file.path(out_dir, "tucker_convergence.rds"))
saveRDS(list(H = H, model_all = model_all, model_selected = model_selected,
             selection = selection, vips_selected = vips_selected),
        file.path(out_dir, "descriptive_models.rds"))
write.csv(data.frame(feature = selection$selected_features),
          file.path(out_dir, "selected_probes.csv"), row.names = FALSE)

# Compare the numerical results reported in the manuscript.
expected <- c(Q2Y = 0.6697, AUC = 0.9715, balanced_accuracy = 0.9309,
              sensitivity = 0.9452, specificity = 0.9167)
checks <- data.frame(
  quantity = c(names(expected), "full_cohort_H", "selected_probes"),
  expected = c(unname(expected), 5, 1132),
  observed = c(unname(metrics[names(expected)]), H, length(selection$selected_features))
)
checks$matches <- round(checks$observed, 4) == round(checks$expected, 4)
write.csv(checks, file.path(out_dir, "reference_comparison.csv"), row.names = FALSE)

inputs <- c(file.path(data_dir, c("X_incomplete.rds", "Y.rds")),
            list.files(file.path(source_dir, "R"), pattern = "\\.[Rr]$",
                       full.names = TRUE, recursive = TRUE),
            file.path(source_dir, "DESCRIPTION"))
write.csv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
          file.path(out_dir, "source_and_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines(c(paste("TensorPLS source:", source_dir), "Seed: 123; five outer and five inner folds",
             "Tucker ranks: (4,4,3); H candidates: 1:10; fixed VIP 99th percentile",
             "Nested CV = held-out performance; full-cohort fit = descriptive"),
           file.path(out_dir, "run_info.txt"))
writeLines(if (all(checks$matches)) "PASS: manuscript rounded metrics reproduced"
           else "DIFFERENCES: inspect reference_comparison.csv",
           file.path(out_dir, "COMPLETED.txt"))
print(checks, row.names = FALSE)
cat("GSE68310 results:", out_dir, "\n")