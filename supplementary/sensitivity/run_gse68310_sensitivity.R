# Fixed analysis settings.
config <- list(
  seed = 123L,
  ncomp = 1:10,
  inner_folds = 5L,
  feature_strategy = "VIP2D_intersect_M2",
  vip_thr = 99,
  threshold = 0.5,
  bootstrap_replicates = 2000L,
  tucker = list(
    ranks = c(4L, 4L, 3L),
    tolerance = 1e-6,
    max_iter = 2000L,
    ridge = 0.001,
    initialization = "hosvd",
    projection_tolerance = 1e-10,
    seed = 123L
  )
)

sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

write_table <- function(x, filename, output_dir) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}

load_analysis_inputs <- function(source_dir, data_dir, primary_dir, output_dir,
                                 project_lib, config) {
  if ("TensorPLS" %in% loadedNamespaces()) {
    stop("Start a fresh R session.")
  }
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  primary_dir <- normalizePath(primary_dir, mustWork = TRUE)
  if (nzchar(project_lib)) {
    .libPaths(c(path.expand(project_lib), .libPaths()))
  }
  pkgload::load_all(source_dir, export_all = FALSE)
  namespace <- asNamespace("TensorPLS")
  loaded_source <- normalizePath(getNamespaceInfo(namespace, "path"))
  stopifnot(identical(loaded_source, source_dir))
  for (name in c("nested_cv_tensorpls", "fit_tucker_imputer",
                 "predict_tucker_imputer")) {
    function_namespace <- environment(getExportedValue("TensorPLS", name))
    stopifnot(identical(function_namespace, namespace))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required.")
  }
  if (file.exists(output_dir)) {
    stop("Output already exists; choose a fresh directory.")
  }
  inputs <- c(
    X = file.path(data_dir, "X_incomplete.rds"),
    Y = file.path(data_dir, "Y.rds"),
    clinical = file.path(
      data_dir, "downloads/GSE68310_SubjectPhenotypes1.txt.gz"
    ),
    primary = file.path(primary_dir, "nested_cv_result.rds"),
    primary_folds = file.path(primary_dir, "nested_cv_folds.rds")
  )
  if (!all(file.exists(inputs))) {
    stop(
      "Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", ")
    )
  }
  protected <- c(
    inputs,
    list.files(file.path(source_dir, "R"), full.names = TRUE, recursive = TRUE),
    file.path(source_dir, c("DESCRIPTION", "NAMESPACE"))
  )
  protected_hashes <- vapply(protected, sha256, character(1))
  X <- readRDS(inputs["X"])
  Y <- readRDS(inputs["Y"])
  ids <- dimnames(X)[[1]]
  y <- as.numeric(Y[, 1])
  stopifnot(
    identical(dim(X), c(109L, 47162L, 4L)),
    sum(is.na(X)) == 990402L,
    identical(dimnames(X)[[3]], c("Day0", "Day2", "Day4", "Day6")),
    identical(rownames(Y), ids),
    !anyDuplicated(ids),
    all(y %in% 0:1),
    sum(y == 1) == 73L
  )
  clinical <- read.delim(gzfile(inputs["clinical"]), check.names = FALSE)
  stopifnot(!anyDuplicated(clinical$SUBJECT_ID))
  clinical <- clinical[match(ids, clinical$SUBJECT_ID), ]
  stopifnot(
    identical(clinical$SUBJECT_ID, ids),
    !anyNA(clinical$Virus_by_PCR),
    !any(clinical$Virus_by_PCR == "none_detected")
  )
  diagnoses <- strsplit(clinical$Virus_by_PCR, "_and_", fixed = TRUE)
  influenza <- vapply(diagnoses, function(z) {
    any(z %in% c("influenza_A_virus", "influenza_B_virus"))
  }, logical(1))
  stopifnot(all(y == as.integer(influenza)))
  coinfection <- lengths(diagnoses) > 1L
  year <- clinical$Enrollment_year
  stopifnot(
    all(year %in% c(2009, 2010)),
    all(
      table(factor(year, levels = c(2009, 2010)), factor(y, levels = 0:1)) ==
        matrix(c(25L, 11L, 31L, 42L), 2)
    )
  )
  primary <- readRDS(inputs["primary"])
  primary_folds <- readRDS(inputs["primary_folds"])
  stopifnot(
    identical(primary$folds, primary_folds),
    primary$settings$imputation == "tucker",
    all(primary$outer_fold_metrics$status == "ok"),
    !length(primary$errors)
  )
  controls <- config$tucker
  for (name in names(controls)) {
    stopifnot(isTRUE(all.equal(
      primary$settings$tucker_control[[name]], controls[[name]]
    )))
  }
  predictions <- primary$outer_predictions
  stopifnot(
    nrow(predictions) == length(ids),
    !anyDuplicated(predictions$subject_id),
    setequal(predictions$subject_id, ids)
  )
  predictions <- predictions[match(ids, predictions$subject_id), ]
  stopifnot(
    all(predictions$Y_true == y),
    all(is.finite(predictions$Y_pred)),
    all(predictions$outer_fold == primary_folds$outer[ids])
  )
  list(
    X = X, Y = Y, ids = ids, y = y, year = year, coinfection = coinfection,
    primary = primary, predictions = predictions, source_dir = source_dir,
    data_dir = data_dir, primary_dir = primary_dir,
    protected = protected, protected_hashes = protected_hashes
  )
}

make_seasonal_folds <- function(Y, ids, year, y, season_folds_file, config) {
  season <- match(year, c(2009, 2010))
  inner_folds <- lapply(1:2, function(f) {
    train <- which(season != f)
    TensorPLS::make_tensorpls_folds(
      Y[train, , drop = FALSE],
      outer_folds = config$inner_folds,
      inner_folds = config$inner_folds,
      seed = config$seed,
      subject_ids = ids[train]
    )$outer
  })
  folds <- list(
    outer = stats::setNames(season, ids),
    inner = stats::setNames(inner_folds, c("1", "2")),
    seed = config$seed
  )
  if (nzchar(season_folds_file)) {
    stopifnot(identical(folds, readRDS(season_folds_file)))
  }
  stopifnot(identical(names(folds$outer), ids), all(folds$outer == season))
  for (fold in 1:2) {
    training_ids <- ids[folds$outer != fold]
    test_ids <- ids[folds$outer == fold]
    inner <- folds$inner[[as.character(fold)]]
    stopifnot(
      !length(intersect(training_ids, test_ids)),
      identical(names(inner), training_ids),
      setequal(inner, seq_len(config$inner_folds))
    )
    for (j in seq_len(config$inner_folds)) {
      inner_training_classes <- y[match(names(inner)[inner != j], ids)]
      inner_validation_classes <- y[match(names(inner)[inner == j], ids)]
      stopifnot(
        length(unique(inner_training_classes)) == 2L,
        length(unique(inner_validation_classes)) == 2L
      )
    }
  }
  folds
}

calculate_metrics <- function(truth, score, config) {
  positive <- truth == 1
  stopifnot(all(is.finite(score)), all(c(0, 1) %in% truth))
  sensitivity <- mean(score[positive] >= config$threshold)
  specificity <- mean(score[!positive] < config$threshold)
  auc_numerator <- sum(rank(score)[positive]) -
    sum(positive) * (sum(positive) + 1) / 2
  auc <- auc_numerator / (sum(positive) * sum(!positive))
  c(
    AUROC = auc,
    balanced_accuracy = (sensitivity + specificity) / 2,
    sensitivity = sensitivity,
    specificity = specificity
  )
}

summarize_subgroup <- function(truth, score, config) {
  estimate <- calculate_metrics(truth, score, config)
  set.seed(config$seed)
  strata <- split(seq_along(truth), truth)
  bootstrap <- replicate(config$bootstrap_replicates, {
    index <- unlist(lapply(strata, function(z) {
      z[sample.int(length(z), length(z), replace = TRUE)]
    }), use.names = FALSE)
    calculate_metrics(truth[index], score[index], config)["AUROC"]
  })
  interval <- quantile(bootstrap, c(0.025, 0.975), names = FALSE)
  data.frame(
    N = length(truth), Influenza = sum(truth == 1), Other = sum(truth == 0),
    as.list(estimate), AUROC_low = interval[1], AUROC_high = interval[2]
  )
}

initialize_outputs <- function(data, folds, output_dir, config) {
  source_dir <- data$source_dir
  data_dir <- data$data_dir
  primary_dir <- data$primary_dir
  protected <- data$protected
  protected_hashes <- data$protected_hashes
  controls <- config$tucker
  dir.create(output_dir, recursive = TRUE)
  write_table(
    data.frame(path = protected, SHA256 = unname(protected_hashes)),
    "input_source_sha256.csv", output_dir
  )
  git <- function(...) {
    system2(
      "git", c("-C", shQuote(source_dir), ...), stdout = TRUE, stderr = TRUE
    )
  }
  git_root <- git("rev-parse", "--show-toplevel")
  if (is.null(attr(git_root, "status")) &&
      identical(normalizePath(git_root), source_dir)) {
    git_info <- c(
      paste("Git commit:", git("rev-parse", "HEAD")),
      "Git status --porcelain:", git("status", "--porcelain")
    )
  } else {
    git_info <- c(
      "Git commit: unavailable (source is not a Git repository root).",
      "Dirty working-tree status: unavailable."
    )
  }
  writeLines(c(
    paste("Started:", Sys.time()),
    paste("TensorPLS source:", source_dir),
    paste("TensorPLS version:", packageVersion("TensorPLS")),
    git_info,
    paste("GSE68310 data:", data_dir),
    paste("Primary analysis:", primary_dir),
    "Seasonal folds: regenerated within training year (seed 123)",
    "H=1:10; Q2Y; VIP2D_intersect_M2; percentile=99",
    "Input provenance: unknown"
  ), file.path(output_dir, "run_info.txt"))
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  saveRDS(folds, file.path(output_dir, "cross_season_folds.rds"))
  saveRDS(controls, file.path(output_dir, "tucker_controls.rds"))
}

evaluate_subgroups <- function(data, config) {
  year <- data$year
  coinfection <- data$coinfection
  y <- data$y
  predictions <- data$predictions
  groups <- list(
    Year2009 = year == 2009, Year2010 = year == 2010,
    Monoinfections = !coinfection, Coinfections = coinfection
  )
  subgroups <- do.call(rbind, lapply(names(groups), function(name) {
    index <- groups[[name]]
    cbind(
      subgroup = name,
      summarize_subgroup(y[index], predictions$Y_pred[index], config)
    )
  }))
  subgroups
}

run_cross_season_validation <- function(data, folds, output_dir, config) {
  X <- data$X
  Y <- data$Y
  controls <- config$tucker
  started <- proc.time()[["elapsed"]]
  result <- withCallingHandlers(
    TensorPLS::nested_cv_tensorpls(
      X, Y, outer_folds = 2L, inner_folds = config$inner_folds,
      seed = config$seed, folds = folds, ncomp = config$ncomp,
      feature_strategy = config$feature_strategy,
      vip_thr = config$vip_thr, strategy_mode = "fixed",
      imputation = "tucker", tucker_control = controls,
      input_provenance = "unknown", selection_metric = "Q2Y",
      keep_models = TRUE, verbose = TRUE
    ),
    warning = function(w) {
      cat(
        conditionMessage(w), "\n",
        file = file.path(output_dir, "warnings.txt"), append = TRUE
      )
    }
  )
  saveRDS(result, file.path(output_dir, "cross_season_result.rds"))
  write_table(
    data.frame(elapsed_seconds = proc.time()[["elapsed"]] - started),
    "runtime.csv", output_dir
  )
  predictions <- result$outer_predictions
  write_table(
    predictions[, !vapply(predictions, is.list, logical(1)), drop = FALSE],
    "cross_season_predictions.csv", output_dir
  )
  saveRDS(result$convergence, file.path(output_dir, "convergence.rds"))
  result
}

validate_cross_season_result <- function(result, data, folds, config) {
  X <- data$X
  ids <- data$ids
  y <- data$y
  year <- data$year
  predictions <- result$outer_predictions
  stopifnot(
    identical(result$folds, folds),
    !length(result$errors),
    nrow(result$outer_fold_metrics) == 2L,
    all(result$outer_fold_metrics$status == "ok"),
    nrow(predictions) == 109L,
    !anyDuplicated(predictions$subject_id),
    setequal(predictions$subject_id, ids),
    all(is.finite(predictions$Y_pred))
  )
  season_metrics <- list()
  checks <- list()
  for (fold in 1:2) {
    train <- which(folds$outer != fold)
    test <- which(folds$outer == fold)
    model <- result$models[[fold]]
    imputer <- model$imputer
    stopifnot(
      identical(imputer$training_ids, ids[train]),
      isTRUE(imputer$convergence$converged)
    )
    inner <- folds$inner[[as.character(fold)]]
    expected_ids <- c(
      list(ids[train]),
      lapply(seq_len(config$inner_folds), function(j) {
        names(inner)[inner != j]
      })
    )
    diagnostics <- result$inner_tuning_results[[fold]]$imputation_diagnostics
    stopifnot(length(diagnostics) >= config$inner_folds)
    for (diagnostic in diagnostics) {
      stopifnot(
        isTRUE(diagnostic$convergence$converged),
        sum(vapply(expected_ids, function(z) {
          identical(z, diagnostic$training_ids)
        }, logical(1))) == 1L
      )
    }
    frozen <- digest::digest(imputer)
    held_out <- X[test, , , drop = FALSE]
    completed <- TensorPLS::predict_tucker_imputer(imputer, held_out)
    stopifnot(
      all(is.finite(completed)),
      identical(held_out[!is.na(held_out)], completed[!is.na(held_out)]),
      identical(frozen, digest::digest(imputer))
    )
    training_matrix <- matrix(X[train, , , drop = FALSE], length(train))
    expected_center <- matrix(
      colMeans(training_matrix, na.rm = TRUE), dim(X)[2], dim(X)[3]
    )
    expected_center[!is.finite(expected_center)] <- 0
    center_difference <- max(abs(
      expected_center - imputer$preprocessing$center
    ))
    stopifnot(center_difference < 1e-12)
    fold_predictions <- predictions[match(ids[test], predictions$subject_id), ]
    stopifnot(
      identical(fold_predictions$subject_id, ids[test]),
      all(fold_predictions$outer_fold == fold),
      all(fold_predictions$Y_true == y[test]),
      max(abs(fold_predictions$training_mean - mean(y[train]))) < 1e-12
    )
    prediction_sse <- sum((fold_predictions$Y_true - fold_predictions$Y_pred)^2)
    baseline_sse <- sum(
      (fold_predictions$Y_true - fold_predictions$training_mean)^2
    )
    q2 <- 1 - prediction_sse / baseline_sse
    tuning <- result$inner_tuning_results[[fold]]$tuning
    valid <- which(is.finite(tuning$selection_score))
    chosen <- valid[order(
      -tuning$selection_score[valid], tuning$ncomp[valid]
    )[1]]
    stopifnot(
      result$selected_ncomp_per_outer_fold[fold] == tuning$ncomp[chosen]
    )
    fold_metrics <- summarize_subgroup(
      fold_predictions$Y_true, fold_predictions$Y_pred, config
    )
    stopifnot(
      abs(q2 - result$outer_fold_metrics$Q2Y[fold]) < 1e-12,
      abs(fold_metrics$AUROC - result$outer_fold_metrics$AUC[fold]) < 1e-12,
      abs(fold_metrics$balanced_accuracy -
            result$outer_fold_metrics$balanced_accuracy[fold]) < 1e-12
    )
    season_metrics[[fold]] <- cbind(
      direction = paste(
        unique(year[train]), "to", unique(year[test]), sep = "_"
      ),
      fold_metrics, Q2Y = q2, H = tuning$ncomp[chosen],
      selected_probes = unique(fold_predictions$n_selected_features)
    )
    checks[[fold]] <- data.frame(
      outer_fold = fold,
      converged = imputer$convergence$converged,
      iterations = imputer$convergence$iterations,
      final_delta = imputer$convergence$final_delta,
      training_center_difference = center_difference,
      observed_preserved = TRUE,
      frozen_imputer = TRUE
    )
  }
  list(
    metrics = do.call(rbind, season_metrics), checks = do.call(rbind, checks)
  )
}

save_analysis_results <- function(validation, result, data, output_dir) {
  protected <- data$protected
  protected_hashes <- data$protected_hashes
  write_table(validation$metrics, "cross_season_metrics.csv", output_dir)
  write_table(validation$checks, "integrity_checks.csv", output_dir)
  stopifnot(
    identical(protected_hashes, vapply(protected, sha256, character(1))),
    identical(
      readRDS(file.path(output_dir, "cross_season_folds.rds")), result$folds
    )
  )
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  writeLines(
    paste("Completed:", Sys.time()), file.path(output_dir, "COMPLETED.txt")
  )
}

# Use the inputs and primary results produced by run_all.sh.
source_dir <- Sys.getenv("TENSORPLS_SOURCE_DIR", unset = "../TensorPLS")
data_dir <- Sys.getenv("GSE68310_DATA", unset = "GSE68310_data")
primary_dir <- Sys.getenv(
  "GSE68310_RESULTS", unset = file.path(dirname(data_dir), "GSE68310_analysis")
)
output_dir <- Sys.getenv(
  "GSE68310_SENSITIVITY_OUTPUT",
  unset = file.path(dirname(data_dir), "GSE68310_sensitivity")
)
season_folds_file <- Sys.getenv("GSE68310_SEASON_FOLDS", unset = "")
project_lib <- Sys.getenv("TENSORPLS_R_LIBRARY", "")

data <- load_analysis_inputs(
  source_dir, data_dir, primary_dir, output_dir, project_lib, config
)
folds <- make_seasonal_folds(
  data$Y, data$ids, data$year, data$y, season_folds_file, config
)
initialize_outputs(data, folds, output_dir, config)

subgroups <- evaluate_subgroups(data, config)
write_table(subgroups, "subgroup_metrics.csv", output_dir)

result <- run_cross_season_validation(data, folds, output_dir, config)
validation <- validate_cross_season_result(result, data, folds, config)
save_analysis_results(validation, result, data, output_dir)

print(subgroups)
print(validation$metrics)
