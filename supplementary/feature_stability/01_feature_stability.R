# Summarize probe selections from saved, completed nested CV results.
strategies <- c(V = "VIP2D", M1 = "M1", M2 = "M2",
  V_intersect_M2 = "VIP2D_intersect_M2")

require_valid <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

parse_arguments <- function(args) {
  options <- list(`input-root` = NULL, `benchmark-root` = NULL, `output-dir` = NULL)
  if (identical(args, "--help")) {
    cat("Required: --input-root PATH --benchmark-root PATH --output-dir PATH\n")
    return(NULL)
  }
  require_valid(length(args) > 0L && length(args) %% 2L == 0L,
    "Supply all three named arguments, each followed by a path; use --help.")
  for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    require_valid(startsWith(args[i], "--") && key %in% names(options),
      paste("Unknown argument:", args[i]))
    require_valid(is.null(options[[key]]), paste("Repeated argument:", args[i]))
    require_valid(nzchar(args[i + 1L]) && !startsWith(args[i + 1L], "--"),
      paste("Missing path for", args[i]))
    options[[key]] <- path.expand(args[i + 1L])
  }
  require_valid(all(lengths(options) == 1L), "All three path arguments are required.")
  for (key in c("input-root", "benchmark-root")) {
    require_valid(dir.exists(options[[key]]), paste("Directory not found:", options[[key]]))
    options[[key]] <- normalizePath(options[[key]])
  }
  options
}

read_result <- function(options, dataset, strategy) {
  path <- if (strategy == "V_intersect_M2") {
    file.path(options[["input-root"]], paste0(dataset, "_analysis"), "nested_cv_result.rds")
  } else {
    file.path(options[["benchmark-root"]], "runs", dataset, strategy, "nested_cv_result.rds")
  }
  require_valid(file.exists(path), paste("Missing saved CV result:", path))
  readRDS(path)
}

validate_result <- function(result, reference, strategy, context) {
  check <- function(ok, message) require_valid(ok, paste(context, message, sep = ": "))
  folds <- as.character(1:5)
  metrics <- result$outer_fold_metrics
  sets <- result$selected_features_per_outer_fold
  features <- reference$stability$frequency$feature
  check(is.character(features) && length(features) > 0L && !anyNA(features) &&
    all(nzchar(trimws(features))) && !anyDuplicated(features), "Invalid probe universe.")
  check(is.data.frame(metrics) && nrow(metrics) == 5L &&
    setequal(metrics$fold, 1:5), "Five distinct outer folds are required.")
  check(all(metrics$status == "ok") && !length(result$errors), "Failed or incomplete outer fold.")
  check(identical(names(sets), folds), "Selected sets must be named 1 to 5 in fold order.")
  for (fold in folds) {
    ids <- sets[[fold]]
    check(is.character(ids) && length(ids) > 0L && !anyNA(ids) &&
      all(nzchar(trimws(ids))) && !anyDuplicated(ids) && all(ids %in% features),
      paste("Invalid, empty, duplicated or unknown probe IDs in fold", fold))
  }
  counts <- result$number_of_selected_features_per_fold
  check(identical(names(counts), folds) &&
    identical(as.numeric(counts), as.numeric(lengths(sets))), "Stored feature counts disagree.")
  check(identical(as.numeric(metrics$n_selected_features[match(1:5, metrics$fold)]),
    as.numeric(lengths(sets))), "Outer metric feature counts disagree.")
  check(identical(result$folds, reference$folds) &&
    identical(result$outer_fold_assignments, reference$folds$outer) &&
    identical(result$inner_fold_assignments, reference$folds$inner), "Outer/inner assignments differ.")
  outer <- result$folds$outer
  check(!anyNA(outer) && setequal(unname(outer), 1:5) &&
    length(names(outer)) == length(outer) && !anyDuplicated(names(outer)) &&
    all(nzchar(names(outer))), "Invalid subject-level outer assignments.")
  check(length(result$folds$inner) == 5L, "Five inner partitions are required.")
  for (i in 1:5) {
    inner <- result$folds$inner[[i]]
    check(identical(names(inner), names(outer)[outer != i]) &&
      !anyNA(inner) && setequal(unname(inner), 1:5), paste("Invalid inner partition", i))
  }
  config <- result$settings$configurations
  check(identical(result$settings$strategy_mode, "fixed") && nrow(config) == 10L &&
    all(config$strategy == strategies[[strategy]]) && all(config$thr == 99) &&
    identical(sort(config$ncomp), 1:10), "Unexpected strategy, percentile or component grid.")
  for (field in setdiff(names(reference$settings), "configurations")) {
    check(identical(result$settings[[field]], reference$settings[[field]]),
      paste("Settings differ from the primary analysis:", field))
  }
  predictions <- result$outer_predictions
  check(nrow(predictions) == length(outer) && !anyDuplicated(predictions$subject_id) &&
    setequal(predictions$subject_id, names(outer)) && all(is.finite(predictions$Y_pred)),
    "Missing or non-finite held-out predictions.")
  index <- match(names(outer), predictions$subject_id)
  reference_index <- match(names(outer), reference$outer_predictions$subject_id)
  check(identical(as.integer(predictions$outer_fold[index]), as.integer(outer)) &&
    identical(predictions$Y_true[index], reference$outer_predictions$Y_true[reference_index]),
    "Subject outcomes or outer assignments disagree.")
  if (!is.null(result$stability$frequency)) {
    saved <- result$stability$frequency
    check(!anyDuplicated(saved$feature) && setequal(saved$feature, features),
      "Saved stability probe universe differs.")
  }
  invisible(features)
}

extract_performance <- function(result, context) {
  scores <- result$classification_metrics
  require_valid(is.numeric(result$Q2Y) && length(result$Q2Y) == 1L &&
    is.finite(result$Q2Y), paste(context, "has no valid Q2Y."))
  require_valid(is.numeric(scores) &&
    all(c("AUC", "balanced_accuracy") %in% names(scores)) &&
    all(is.finite(scores[c("AUC", "balanced_accuracy")])),
    paste(context, "has incomplete classification metrics."))
  c(Q2Y = unname(result$Q2Y), AUROC = unname(scores["AUC"]),
    balanced_accuracy = unname(scores["balanced_accuracy"]))
}

compute_stability <- function(result, features) {
  sets <- result$selected_features_per_outer_fold
  membership <- vapply(sets, function(ids) features %in% ids, logical(length(features)))
  frequency <- data.frame(feature = features, count = rowSums(membership))
  frequency$frequency <- frequency$count / 5
  jaccard <- matrix(NA_real_, 5, 5, dimnames = list(names(sets), names(sets)))
  for (i in 1:5) for (j in 1:5) {
    jaccard[i, j] <- length(intersect(sets[[i]], sets[[j]])) /
      length(union(sets[[i]], sets[[j]]))
  }
  list(counts = lengths(sets), frequency = frequency, jaccard = jaccard)
}

build_tables <- function(options) {
  summaries <- counts <- pairs <- frequencies <- distributions <- list()
  for (dataset in c("GSE68310", "GSE182616")) {
    reference <- read_result(options, dataset, "V_intersect_M2")
    for (strategy in names(strategies)) {
      result <- if (strategy == "V_intersect_M2") reference else read_result(options, dataset, strategy)
      features <- validate_result(result, reference, strategy, paste(dataset, strategy))
      stability <- compute_stability(result, features)
      performance <- extract_performance(result, paste(dataset, strategy))
      primary <- strategy == "V_intersect_M2"
      key <- paste(dataset, strategy)
      summaries[[key]] <- data.frame(dataset, strategy,
        Q2Y = performance["Q2Y"], AUROC = performance["AUROC"],
        balanced_accuracy = performance["balanced_accuracy"],
        median_selected_features = median(stability$counts),
        mean_pairwise_jaccard = mean(stability$jaccard[upper.tri(stability$jaccard)]),
        union_features = if (primary) sum(stability$frequency$count > 0) else NA_integer_,
        intersection_features = if (primary) sum(stability$frequency$count == 5) else NA_integer_)
      counts[[key]] <- data.frame(dataset, strategy, fold = 1:5,
        selected_features = as.integer(stability$counts))
      indices <- which(upper.tri(stability$jaccard), arr.ind = TRUE)
      pairs[[key]] <- data.frame(dataset, strategy, fold1 = indices[, 1], fold2 = indices[, 2],
        jaccard = stability$jaccard[indices])
      if (primary) {
        frequencies[[dataset]] <- data.frame(dataset, stability$frequency)
        distributions[[dataset]] <- data.frame(dataset, folds_selected = 0:5,
          n_features = as.integer(table(factor(stability$frequency$count, levels = 0:5))))
      }
    }
  }
  tables <- list(strategy_stability_summary = summaries, strategy_fold_counts = counts,
    strategy_jaccard_pairs = pairs, primary_selection_frequency = frequencies,
    primary_frequency_distribution = distributions)
  lapply(tables, function(rows) { x <- do.call(rbind, rows); rownames(x) <- NULL; x })
}

main <- function() {
  options <- parse_arguments(commandArgs(trailingOnly = TRUE))
  if (is.null(options)) return(invisible(NULL))
  output <- options[["output-dir"]]
  require_valid(!file.exists(output), "Output directory already exists; choose a new directory.")
  tables <- build_tables(options)
  require_valid(dir.create(output, recursive = TRUE), "Cannot create output directory.")
  for (name in names(tables)) {
    write.csv(tables[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  }
  print(tables$strategy_stability_summary, row.names = FALSE)
  cat("Saved five CSV tables to", normalizePath(output), "\n")
}

if (sys.nframe() == 0L) main()
