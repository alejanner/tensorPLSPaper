# GSE182616:  the primary nested CV and full-cohort descriptive fit.
source_dir <- normalizePath(Sys.getenv("TENSORPLS_SOURCE_DIR"), mustWork = TRUE)
data_dir <- normalizePath(Sys.getenv("GSE182616_DATA"), mustWork = TRUE)
out_dir <- Sys.getenv("GSE182616_RESULTS")
if (!nzchar(out_dir) || dir.exists(out_dir)) stop("Choose a new GSE182616_RESULTS directory.")

project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
if (!requireNamespace("pkgload", quietly = TRUE)) stop("The pkgload package is required.")
pkgload::load_all(source_dir, export_all = FALSE)
stopifnot(identical(normalizePath(getNamespaceInfo("TensorPLS", "path")), source_dir))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

expression <- readRDS(file.path(data_dir, "expression_as_deposited.rds"))
metadata <- read.csv(file.path(data_dir, "sample_metadata_audited.csv"),
                     colClasses = c(subject_id = "character", subject_token = "character"))
annotation <- read.csv(file.path(data_dir, "feature_annotation.csv"), check.names = FALSE)

times <- c(0, 2, 4, 8, 12, 24)
visits <- metadata[metadata$time_hours %in% times, , drop = FALSE]
ids <- sort(unique(visits$subject_id))
features <- rownames(expression)[rownames(expression) %in% annotation$ID[annotation$retained]]
stopifnot(!anyDuplicated(visits[c("subject_id", "time_hours")]),
          length(features) == 50599L)
input <- data.frame(visits[c("subject_id", "time_hours")],
                    t(expression[features, visits$gsm, drop = FALSE]), check.names = FALSE)
X <- TensorPLS::prepare_omics(
  input, id_col = "subject_id", time_col = "time_hours", transpose = "never",
  feature_cols = features, subjects = ids, time_points = times,
  coercion_mode = "custom", id_as = "character", time_as = "numeric",
  features_numeric = FALSE
)
names(dimnames(X)) <- NULL
subject_metadata <- metadata[match(ids, metadata$subject_id), , drop = FALSE]
Y <- factor(ifelse(subject_metadata$tbsa >= 20, "Severe", "LessSevere"),
            levels = c("LessSevere", "Severe"))
names(Y) <- ids
rm(expression, input)
stopifnot(identical(dim(X), c(66L, 50599L, 6L)),
          identical(dimnames(X)[[3]], as.character(times)),
          sum(is.na(X)) == 91L * 50599L,
          identical(as.integer(table(Y)), c(32L, 34L)),
          identical(dimnames(X)[[1]], names(Y)))

reference_tensor <- Sys.getenv("GSE182616_REFERENCE_TENSOR", "")
if (nzchar(reference_tensor) && file.exists(reference_tensor)) {
  original_X <- readRDS(reference_tensor)
  stopifnot(identical(dim(X), dim(original_X)),
            identical(dimnames(X), dimnames(original_X)),
            identical(as.vector(X), as.vector(original_X)))
  rm(original_X)
  cat("GSE182616 tensor matches the saved reference exactly.\n")
}

seed <- 123L
components <- 1:10
strategy <- "VIP2D_intersect_M2"
threshold <- 99

# Recover the original folds, checking the regenerated assignments.
generated <- TensorPLS::make_tensorpls_folds(
  Y, outer_folds = 5L, inner_folds = 5L, seed = seed, subject_ids = ids
)
reference_folds <- Sys.getenv("GSE182616_REFERENCE_FOLDS", "")
folds <- generated
if (nzchar(reference_folds) && file.exists(reference_folds)) {
  folds <- readRDS(reference_folds)
  stopifnot(identical(folds, generated))
  cat("GSE182616: regenerated folds match the original folds.\n")
} else {
  warning("Original GSE182616 folds unavailable; using generated folds.")
}

dir.create(out_dir, recursive = TRUE)
saveRDS(X, file.path(out_dir, "X_incomplete.rds"), compress = FALSE)
saveRDS(Y, file.path(out_dir, "Y.rds"))
saveRDS(folds, file.path(out_dir, "nested_cv_folds.rds"))

# Predictive evaluation with training-only mean imputation.
cv <- TensorPLS::nested_cv_tensorpls(
  X = X, Y = Y, folds = folds, outer_folds = 5L, inner_folds = 5L,
  seed = seed, ncomp = components, feature_strategy = strategy,
  vip_thr = threshold, strategy_mode = "fixed", imputation = "training_mean",
  input_provenance = "unknown", selection_metric = "Q2Y",
  keep_models = FALSE, verbose = TRUE
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

reference_cv <- Sys.getenv("GSE182616_REFERENCE_CV", "")
if (nzchar(reference_cv) && file.exists(reference_cv)) {
  old <- readRDS(reference_cv)
  stopifnot(identical(old$folds, folds), setequal(old$outer_predictions$subject_id, pred$subject_id))
  old_pred <- old$outer_predictions$Y_pred[match(pred$subject_id, old$outer_predictions$subject_id)]
  max_difference <- max(abs(pred$Y_pred - old_pred))
  write.csv(data.frame(max_absolute_OOF_score_difference = max_difference),
            file.path(out_dir, "reference_prediction_comparison.csv"), row.names = FALSE)
  cat("GSE182616: max absolute difference from saved OOF scores:", max_difference, "\n")
  rm(old)
}

# Full-cohort interpretation is separate from held-out predictive assessment.
development <- TensorPLS::tune_tensorpls_components(
  X = X, Y = Y, ncomp = components, feature_strategy = strategy,
  vip_thr = threshold,
  cv_folds = 5L, folds = folds$outer, seed = seed,
  imputation = "training_mean", input_provenance = "unknown",
  selection_metric = "Q2Y", verbose = TRUE
)
H <- development$selected_ncomp
saveRDS(development, file.path(out_dir, "component_tuning.rds"))

X_complete <- TensorPLS::ImputeTensor(X, method = "mean")
stopifnot(all(is.finite(X_complete)),
          identical(X_complete[!is.na(X)], X[!is.na(X)]))
model_all <- TensorPLS::fit_tensorpls_pls(X_complete, Y, ncomp = H)
vips_all <- TensorPLS::tensorpls_vips(model_all)
selection <- TensorPLS::vip_feature_sets(vips_all, strategy = strategy, thr = threshold)
X_selected <- X_complete[, selection$selected_indices, , drop = FALSE]
model_selected <- TensorPLS::fit_tensorpls_pls(X_selected, Y, ncomp = H)
vips_selected <- TensorPLS::tensorpls_vips(model_selected)

saveRDS(list(H = H, model_all = model_all, model_selected = model_selected,
             selection = selection, vips_selected = vips_selected),
        file.path(out_dir, "descriptive_models.rds"))
write.csv(data.frame(feature = selection$selected_features),
          file.path(out_dir, "selected_probes.csv"), row.names = FALSE)

expected <- c(Q2Y = 0.2939, AUC = 0.8272, balanced_accuracy = 0.7564,
              sensitivity = 0.7941, specificity = 0.7188)
checks <- data.frame(
  quantity = c(names(expected), "full_cohort_H", "selected_probes"),
  expected = c(unname(expected), 2, 2175),
  observed = c(unname(metrics[names(expected)]), H, length(selection$selected_features))
)
checks$matches <- round(checks$observed, 4) == round(checks$expected, 4)
write.csv(checks, file.path(out_dir, "reference_comparison.csv"), row.names = FALSE)

inputs <- c(file.path(data_dir, c("expression_as_deposited.rds",
                                    "sample_metadata_audited.csv", "feature_annotation.csv")),
            list.files(file.path(source_dir, "R"), pattern = "\\.[Rr]$",
                       full.names = TRUE, recursive = TRUE),
            file.path(source_dir, "DESCRIPTION"))
write.csv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
          file.path(out_dir, "source_and_input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines(c(paste("TensorPLS source:", source_dir), "Seed: 123; five outer and five inner folds",
             "Mean imputation: training-only in nested CV; H candidates: 1:10; VIP 99th percentile",
             "Nested CV = held-out performance; full-cohort fit = descriptive"),
           file.path(out_dir, "run_info.txt"))
writeLines(if (all(checks$matches)) "PASS: manuscript rounded metrics reproduced"
           else "DIFFERENCES: inspect reference_comparison.csv",
           file.path(out_dir, "COMPLETED.txt"))
print(checks, row.names = FALSE)
cat("GSE182616 results:", out_dir, "\n")