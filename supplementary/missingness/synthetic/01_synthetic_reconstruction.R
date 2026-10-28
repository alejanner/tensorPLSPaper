# Synthetic held-out visit reconstruction using complete training tensors.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 01_synthetic_reconstruction.R pilot|full NEW_OUTPUT_DIR [--preflight]"
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0) }
preflight <- "--preflight" %in% args
args <- args[args != "--preflight"]
if (length(args) != 2L) stop(usage, call. = FALSE)
mode <- args[1L]
if (!mode %in% c("pilot", "full")) stop(usage, call. = FALSE)
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

out <- path.expand(args[2L])
if (file.exists(out)) stop("Output already exists: ", out)
reference <- Sys.getenv("SYNTHETIC_ORIGINAL_ROOT", "")
if (nzchar(reference)) reference <- normalizePath(path.expand(reference), mustWork = TRUE)
if (preflight) { cat("Preflight passed; no simulation or fitting.\n"); quit(status = 0) }
dir.create(out, recursive = TRUE)
cat("Output:", out, "\n")

config <- list(dimensions = c(109L, 1000L, 4L), time = seq(0, 1, length.out = 4),
  levels = c(.1, .2, .3, .4), mask_seeds = c(5001L, 5002L, 5003L),
  fold_seed = 123L, rng = c("Mersenne-Twister", "Inversion", "Rejection"),
  baseline_sd = 1, mean_slope_sd = .2, mean_response_sd = .3,
  sensitivity_range = c(.5, 1.5),
  tucker = list(ranks = c(4L, 4L, 3L), ridge = .001, tolerance = 1e-6,
    max_iter = 2000L, initialization = "hosvd",
    projection_tolerance = 1e-10, seed = 123L))
scenarios <- data.frame(scenario = c("A", "B", "C"),
  description = c("Strong structured signal", "Weak individual signal",
    "Heterogeneous trajectories"),
  modules = c(1L, 1L, 20L), offset_sd = c(.4, .032, .4),
  amplitude_sd = c(1, .08, 1), slope_sd = c(.5, .04, .5),
  peak_sd = c(.04, .04, .25), width_min = c(.35, .35, .12),
  width_max = c(.35, .35, .35), noise_sd = c(.3, .6, .3))
seed_table <- expand.grid(replicate = 1:3, scenario = scenarios$scenario,
  stringsAsFactors = FALSE)
seed_table$simulation_seed <- match(seed_table$scenario, scenarios$scenario) *
  1000L + seed_table$replicate

set_rng <- function(seed) {
  do.call(RNGkind, as.list(config$rng))
  set.seed(seed)
}
# Generate longitudinal datasets and held-out visit masks.
generate_dataset <- function(parameters, seed) {
  set_rng(seed)
  n <- config$dimensions[1]
  p <- config$dimensions[2]
  k <- config$dimensions[3]
  modules <- parameters$modules
  baseline <- rnorm(p, sd = config$baseline_sd)
  mean_slope <- rnorm(p, sd = config$mean_slope_sd)
  mean_response <- rnorm(p, sd = config$mean_response_sd)
  feature_means <- outer(baseline, rep(1, k)) +
    outer(mean_slope, config$time) + outer(mean_response, sin(pi * config$time))
  offset_loading <- runif(p, config$sensitivity_range[1], config$sensitivity_range[2])
  response_loading <- runif(p, config$sensitivity_range[1], config$sensitivity_range[2])
  slope_loading <- runif(p, config$sensitivity_range[1], config$sensitivity_range[2])
  module <- sample(rep(seq_len(modules), length.out = p))
  offsets <- matrix(rnorm(n * modules, sd = parameters$offset_sd), n, modules)
  amplitudes <- matrix(rnorm(n * modules, sd = parameters$amplitude_sd), n, modules)
  slopes <- matrix(rnorm(n * modules, sd = parameters$slope_sd), n, modules)
  peaks <- matrix(.5 + rnorm(n * modules, sd = parameters$peak_sd), n, modules)
  widths <- matrix(runif(n * modules, parameters$width_min, parameters$width_max),
    n, modules)
  latent <- array(0, config$dimensions, dimnames = list(
    sprintf("Subject%03d", seq_len(n)), sprintf("Feature%04d", seq_len(p)),
    paste0("Time", seq_len(k))))
  for (t in seq_len(k)) {
    response <- amplitudes * exp(-.5 * ((config$time[t] - peaks) / widths)^2)
    deviations <- sweep(offsets[, module, drop = FALSE], 2, offset_loading, "*") +
      sweep(response[, module, drop = FALSE], 2, response_loading, "*") +
      sweep(slopes[, module, drop = FALSE], 2, slope_loading, "*") * config$time[t]
    latent[, , t] <- sweep(deviations, 2, feature_means[, t], "+")
  }
  noise <- array(rnorm(length(latent), sd = parameters$noise_sd), dim(latent))
  X <- latent + noise
  dimnames(X) <- dimnames(latent)
  list(X = X, latent = latent, feature_time_means = feature_means,
    module = module, offsets = offsets, amplitudes = amplitudes, slopes = slopes,
    peaks = peaks, widths = widths,
    loadings = list(offset = offset_loading, response = response_loading,
      slope = slope_loading), parameters = parameters, seed = seed)
}

make_folds <- function(ids) {
  set_rng(config$fold_seed)
  setNames(sample(rep(1:5, length.out = length(ids))), ids)
}

make_mask <- function(n, k, level, seed, ids, times) {
  set_rng(seed)
  candidates <- which(matrix(TRUE, n, k), arr.ind = TRUE)
  permutation <- sample.int(nrow(candidates))
  remaining <- rep(k, n)
  target <- as.integer(floor(level * n * k + .5))
  hidden <- matrix(FALSE, n, k, dimnames = list(ids, times))
  for (index in permutation) {
    if (sum(hidden) == target) break
    subject <- candidates[index, 1]
    time <- candidates[index, 2]
    if (remaining[subject] <= 1L) next
    hidden[subject, time] <- TRUE
    remaining[subject] <- remaining[subject] - 1L
  }
  stopifnot(all(remaining >= 1L))
  list(hidden = hidden, target = target, achieved = sum(hidden), seed = seed,
    rng = config$rng, permutation = permutation)
}


mean_fit_fn <- getFromNamespace(".nested_imputer", "TensorPLS")
mean_apply <- getFromNamespace(".nested_impute", "TensorPLS")

fit_rows <- list()
case_rows <- list()
subject_rows <- list()
selection <- seed_table[if (mode == "pilot") seed_table$replicate == 1L else TRUE, ]

for (r in seq_len(nrow(selection))) {
  meta <- selection[r, ]
  label <- paste0(meta$scenario, "_replicate_", meta$replicate)
  params <- scenarios[scenarios$scenario == meta$scenario, ]
  dataset <- generate_dataset(params, meta$simulation_seed)
  X <- dataset$X
  latent <- dataset$latent
  folds <- make_folds(dimnames(X)[[1]])
  stopifnot(identical(dim(X), config$dimensions), !anyNA(X),
            all(is.finite(X)), all(folds %in% 1:5))

  data_dir <- file.path(out, "datasets", label)
  fits_dir <- file.path(out, "fits", label)
  mask_dir <- file.path(out, "masks", label)
  for (d in c(data_dir, fits_dir, mask_dir))
    dir.create(d, recursive = TRUE)
  saveRDS(dataset, file.path(data_dir, "dataset.rds"))
  saveRDS(folds, file.path(data_dir, "folds.rds"))

  if (nzchar(reference)) {
    stopifnot(identical(dataset, readRDS(file.path(reference, "datasets", label, "dataset.rds"))),
              identical(folds, readRDS(file.path(reference, "datasets", label, "folds.rds"))))
  }

  for (fold in if (mode == "pilot") 1L else 1:5) {
    train <- X[folds != fold, , , drop = FALSE]
    test <- X[folds == fold, , , drop = FALSE]
    truth_latent <- latent[folds == fold, , , drop = FALSE]
    mean_fit <- mean_fit_fn(matrix(train, dim(train)[1]), "training_mean")

    error_message <- ""
    start <- proc.time()[["elapsed"]]
    tucker <- tryCatch(
      do.call(TensorPLS::fit_tucker_imputer,
              c(list(X_train = train), config$tucker)),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    seconds <- proc.time()[["elapsed"]] - start
    converged <- !is.null(tucker) && isTRUE(tucker$convergence$converged)
    if (!is.null(tucker))
      stopifnot(identical(tucker$training_ids, dimnames(train)[[1]]),
                identical(tucker$X_train_imputed, train))

    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      scenario = meta$scenario, replicate = meta$replicate, fold = fold,
      converged = converged,
      iterations = if (is.null(tucker)) NA_integer_ else tucker$convergence$iterations,
      final_delta = if (is.null(tucker)) NA_real_ else tucker$convergence$final_delta,
      fit_seconds = seconds, error = error_message
    )
    saveRDS(list(mean = mean_fit, tucker = tucker),
      file.path(fits_dir, paste0("fold_", fold, ".rds")))
    cat(label, "fold", fold, "Tucker converged:", converged, "\n")

    for (seed in if (mode == "pilot") config$mask_seeds[1] else config$mask_seeds) {
      for (level in if (mode == "pilot") c(.1, .3) else config$levels) {
        mask <- make_mask(dim(test)[1], dim(test)[3], level, seed,
                          dimnames(test)[[1]], dimnames(test)[[3]])
        stopifnot(mask$achieved == mask$target, all(rowSums(!mask$hidden) >= 1L))
        tag <- sprintf("fold_%d_seed_%d_mask_%02d", fold, seed, round(level * 100))
        saveRDS(mask, file.path(mask_dir, paste0(tag, "_mask.rds")))
        if (nzchar(reference))
          stopifnot(identical(mask, readRDS(file.path(reference, "full", label,
                                                       paste0(tag, "_mask.rds")))))

        evaluation <- array(FALSE, dim(test), dimnames(test))
        for (ti in seq_len(dim(test)[3]))
          evaluation[mask$hidden[, ti], , ti] <- TRUE
        masked <- test
        masked[evaluation] <- NA_real_
        stopifnot(sum(evaluation) == mask$target * dim(test)[2],
                  all(is.na(masked[evaluation])),
                  all(!is.na(masked[!evaluation])))

        mean_pred <- array(
          mean_apply(mean_fit, matrix(masked, dim(masked)[1])),
          dim(masked), dimnames(masked)
        )
        stopifnot(identical(mean_pred[!evaluation], test[!evaluation]),
                  all(is.finite(mean_pred[evaluation])))

        status <- if (converged) "ok" else "non_converged"
        tucker_pred <- NULL
        if (converged) {
          prediction <- tryCatch(
            TensorPLS::predict_tucker_imputer(tucker, masked, diagnostics = TRUE),
            error = function(e) {
              status <<- paste0("prediction_failed: ", conditionMessage(e))
              NULL
            }
          )
          if (!is.null(prediction)) {
            if (any(prediction$diagnostics$rank_deficient)) status <- "rank_deficient"
            if (!identical(prediction$X[!evaluation], test[!evaluation]) ||
                !all(is.finite(prediction$X[evaluation]))) status <- "invalid_prediction"
            if (status == "ok") tucker_pred <- prediction$X
          }
        }

        n <- sum(evaluation)
        cases <- data.frame(scenario = meta$scenario, replicate = meta$replicate,
          fold = fold, seed = seed, level = level, hidden_visits = mask$achieved,
          n = n, tucker_status = status)
        for (target in c("noisy", "latent")) {
          truth <- if (target == "noisy") test else truth_latent
          em <- mean_pred[evaluation] - truth[evaluation]
          cases[[paste0("mean_SSE_", target)]] <- sum(em^2)
          cases[[paste0("mean_RMSE_", target)]] <- sqrt(mean(em^2))
          if (status == "ok") {
            et <- tucker_pred[evaluation] - truth[evaluation]
            cases[[paste0("tucker_SSE_", target)]] <- sum(et^2)
            cases[[paste0("tucker_RMSE_", target)]] <- sqrt(mean(et^2))
          } else {
            cases[[paste0("tucker_SSE_", target)]] <- NA_real_
            cases[[paste0("tucker_RMSE_", target)]] <- NA_real_
          }
        }
        case_rows[[length(case_rows) + 1L]] <- cases

        for (i in which(rowSums(mask$hidden) > 0L)) {
          index <- evaluation[i, , ]
          for (target in c("noisy", "latent")) {
            truth <- if (target == "noisy") test else truth_latent
            em <- mean_pred[i, , ][index] - truth[i, , ][index]
            et <- if (status == "ok") tucker_pred[i, , ][index] - truth[i, , ][index]
            subject_rows[[length(subject_rows) + 1L]] <- data.frame(
              scenario = meta$scenario, replicate = meta$replicate,
              fold = fold, seed = seed, level = level,
              subject_id = dimnames(test)[[1]][i],
              hidden_visits = sum(mask$hidden[i, ]),
              remaining_visits = sum(!mask$hidden[i, ]),
              target = target, n = length(em),
              mean_SSE = sum(em^2), mean_RMSE = sqrt(mean(em^2)),
              tucker_SSE = if (status == "ok") sum(et^2) else NA_real_,
              tucker_RMSE = if (status == "ok") sqrt(mean(et^2)) else NA_real_,
              tucker_status = status
            )
          }
        }
      }
    }
    # Save completed cases after every fold.
    write.csv(do.call(rbind, fit_rows), file.path(out, "fit_summary.csv"), row.names = FALSE)
    write.csv(do.call(rbind, case_rows), file.path(out, "reconstruction_metrics.csv"), row.names = FALSE)
  }
}

fits <- do.call(rbind, fit_rows)
cases <- do.call(rbind, case_rows)
subjects <- do.call(rbind, subject_rows)
write.csv(subjects, file.path(out, "subject_level_errors.csv"), row.names = FALSE)

# Paired RMSE: Mean and Tucker are evaluated on exactly the same valid masks.
summarise_paired <- function(z) {
  good <- z[z$tucker_status == "ok", , drop = FALSE]
  if (!nrow(good)) stop("No valid paired cases")
  n <- sum(good$n)
  data.frame(cases = nrow(z), paired_cases = nrow(good),
    failed_tucker = nrow(z) - nrow(good),
    mean_RMSE_noisy = sqrt(sum(good$mean_SSE_noisy) / n),
    tucker_RMSE_noisy = sqrt(sum(good$tucker_SSE_noisy) / n),
    mean_RMSE_latent = sqrt(sum(good$mean_SSE_latent) / n),
    tucker_RMSE_latent = sqrt(sum(good$tucker_SSE_latent) / n))
}

by_level <- do.call(rbind, lapply(split(cases,
  interaction(cases$scenario, cases$level, drop = TRUE)), function(z) {
    cbind(z[1, c("scenario", "level")], summarise_paired(z))
  }))
by_replicate <- do.call(rbind, lapply(split(cases,
  interaction(cases$scenario, cases$replicate, cases$level, drop = TRUE)), function(z) {
    cbind(z[1, c("scenario", "replicate", "level")], summarise_paired(z))
  }))
write.csv(by_level, file.path(out, "paired_RMSE_by_level.csv"), row.names = FALSE)
write.csv(by_replicate, file.path(out, "paired_RMSE_by_replicate.csv"), row.names = FALSE)

valid_subjects <- subjects[subjects$tucker_status == "ok", , drop = FALSE]
subject_stats <- function(z) data.frame(
  evaluations = nrow(z),
  tucker_wins = sum(z$tucker_RMSE < z$mean_RMSE),
  tucker_win_percent = 100 * mean(z$tucker_RMSE < z$mean_RMSE),
  median_delta_RMSE = median(z$tucker_RMSE - z$mean_RMSE),
  p95_tucker_RMSE = unname(quantile(z$tucker_RMSE, 0.95))
)
by_subject_level <- do.call(rbind, lapply(split(valid_subjects,
  interaction(valid_subjects$scenario, valid_subjects$level, valid_subjects$target,
              drop = TRUE)), function(z)
    cbind(z[1, c("scenario", "level", "target")], subject_stats(z))))
by_remaining <- do.call(rbind, lapply(split(valid_subjects,
  interaction(valid_subjects$scenario, valid_subjects$remaining_visits,
              valid_subjects$target, drop = TRUE)), function(z)
    cbind(z[1, c("scenario", "remaining_visits", "target")], subject_stats(z))))
write.csv(by_subject_level, file.path(out, "subject_summary_by_level.csv"), row.names = FALSE)
write.csv(by_remaining, file.path(out, "subject_summary_by_remaining_visits.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("\nNew Tucker fits:", nrow(fits),
    "\nConverged:", sum(fits$converged),
    "\nMasks:", nrow(cases),
    "\nPaired valid:", sum(cases$tucker_status == "ok"), "\n")
print(by_level, row.names = FALSE)

# An optional reference check compares the reproduced overall case scores.
if (nzchar(reference) && mode == "full") {
  original <- read.csv(file.path(reference, "full", "reconstruction_metrics.csv"))
  original <- original[original$time == "all", , drop = FALSE]
  new <- do.call(rbind, lapply(c("training_mean", "tucker"), function(method) {
    do.call(rbind, lapply(c("noisy", "latent"), function(target) {
      z <- cases[, c("scenario", "replicate", "fold", "seed", "level", "n", "tucker_status")]
      z$method <- method
      z$target <- target
      z$SSE <- cases[[paste0(if (method == "tucker") "tucker_SSE_" else "mean_SSE_", target)]]
      z$RMSE <- cases[[paste0(if (method == "tucker") "tucker_RMSE_" else "mean_RMSE_", target)]]
      z
    }))
  }))
  keys <- c("scenario", "replicate", "fold", "seed", "level", "method", "target")
  check <- merge(new, original[, c(keys, "status", "n", "SSE", "RMSE")],
                 by = keys, suffixes = c("_new", "_old"))
  stopifnot(nrow(check) == nrow(new),
            all(check$n_new[check$method == "training_mean"] ==
                check$n_old[check$method == "training_mean"]))
  reproduced_ok <- check$method == "training_mean" | check$tucker_status == "ok"
  original_ok <- check$status == "ok"
  stopifnot(identical(reproduced_ok, original_ok))
  valid <- reproduced_ok & original_ok
  stopifnot(all(check$n_new[valid] == check$n_old[valid]))
  difference <- max(abs(check$RMSE_new[valid] - check$RMSE_old[valid]))
  write.csv(check, file.path(out, "reference_comparison.csv"), row.names = FALSE)
  cat("Reference max RMSE difference:", difference, "\n")
  if (difference <= 1e-8) writeLines("PASS: overall RMSE reproduced", file.path(out, "PASS.txt"))
  else warning("Reference case RMSEs differ: inspect reference_comparison.csv")
}
