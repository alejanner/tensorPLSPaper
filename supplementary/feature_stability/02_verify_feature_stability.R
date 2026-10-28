# Verify CSV values against saved probe sets and native stability objects.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(script)), "01_feature_stability.R"))

verify_tables <- function(options) {
  names <- c("strategy_stability_summary", "strategy_fold_counts", "strategy_jaccard_pairs",
    "primary_selection_frequency", "primary_frequency_distribution")
  tables <- setNames(lapply(names, function(name) {
    path <- file.path(options[["output-dir"]], paste0(name, ".csv"))
    require_valid(file.exists(path), paste("Missing CSV:", path))
    classes <- if (name == "primary_selection_frequency") c(feature = "character") else NA
    read.csv(path, colClasses = classes, check.names = FALSE)
  }), names)
  expected_columns <- list(
    c("dataset", "strategy", "Q2Y", "AUROC", "balanced_accuracy",
      "median_selected_features", "mean_pairwise_jaccard",
      "union_features", "intersection_features"),
    c("dataset", "strategy", "fold", "selected_features"),
    c("dataset", "strategy", "fold1", "fold2", "jaccard"),
    c("dataset", "feature", "count", "frequency"),
    c("dataset", "folds_selected", "n_features"))
  for (i in seq_along(tables)) {
    require_valid(identical(colnames(tables[[i]]), expected_columns[[i]]),
      paste("Unexpected CSV schema:", names[i]))
    require_valid(all(tables[[i]]$dataset %in% c("GSE68310", "GSE182616")), "Unknown dataset in CSV.")
  }
  require_valid(identical(vapply(tables[c(1, 2, 3, 5)], nrow, integer(1)),
    setNames(c(8L, 40L, 80L, 12L), names[c(1, 2, 3, 5)])), "Unexpected CSV row counts.")
  reference_medians <- c(1485, 500, 1304, 1122, 2684, 542, 2538, 2194)
  reference_jaccards <- c(.564644, .664511, .566827, .583822, .247875, .322826, .258883, .243913)
  reference_q2 <- c(.663965, .686649, .663517, .669742,
    .285090, .261601, .250226, .293936)
  reference_auc <- c(.969939, .974125, .968798, .971461,
    .823529, .826287, .803309, .827206)
  reference_ba <- c(.930936, .958524, .917047, .930936,
    .756434, .758272, .727941, .756434)
  reference_index <- 0L
  total_features <- 0L
  for (dataset in c("GSE68310", "GSE182616")) {
    reference <- read_result(options, dataset, "V_intersect_M2")
    total_features <- total_features + nrow(reference$stability$frequency)
    for (strategy in names(strategies)) {
      context <- paste(dataset, strategy)
      check <- function(ok, message) require_valid(ok, paste(context, message, sep = ": "))
      result <- if (strategy == "V_intersect_M2") reference else read_result(options, dataset, strategy)
      features <- validate_result(result, reference, strategy, context)
      sets <- result$selected_features_per_outer_fold
      # Reconstruct counts by tabulating the five saved selections together.
      frequency <- as.integer(table(factor(unlist(sets, use.names = FALSE), levels = features)))
      membership <- sapply(sets, function(ids) features %in% ids)
      shared <- crossprod(membership)
      union_counts <- outer(lengths(sets), lengths(sets), "+") - shared
      jaccard <- shared / union_counts
      selected_rows <- function(tab) tab[tab$dataset == dataset & tab$strategy == strategy, , drop = FALSE]
      counts <- selected_rows(tables$strategy_fold_counts)
      check(nrow(counts) == 5L && setequal(counts$fold, 1:5) &&
        all(counts$selected_features[match(1:5, counts$fold)] == lengths(sets)), "Fold-count CSV mismatch.")
      pairs <- selected_rows(tables$strategy_jaccard_pairs)
      pair_keys <- paste(pairs$fold1, pairs$fold2)
      check(nrow(pairs) == 10L && !anyDuplicated(pair_keys) &&
        all(pairs$fold1 %in% 1:5 & pairs$fold2 %in% 1:5 & pairs$fold1 < pairs$fold2),
        "Expected ten distinct fold pairs.")
      check(all(is.finite(pairs$jaccard)) &&
        max(abs(pairs$jaccard - jaccard[cbind(pairs$fold1, pairs$fold2)])) < 1e-12,
        "Pairwise Jaccard CSV mismatch.")
      summary <- selected_rows(tables$strategy_stability_summary)
      mean_jaccard <- mean(jaccard[upper.tri(jaccard)])
      check(nrow(summary) == 1L && summary$median_selected_features == median(lengths(sets)) &&
        abs(summary$mean_pairwise_jaccard - mean_jaccard) < 1e-12, "Summary CSV mismatch.")
      performance <- extract_performance(result, context)
      csv_performance <- unlist(summary[1, names(performance)], use.names = FALSE)
      check(all(is.finite(csv_performance)) &&
        max(abs(csv_performance - performance)) < 1e-12,
        "Predictive metrics CSV mismatch.")

      # Recalculate pooled predictive metrics from held-out predictions.
      pred <- result$outer_predictions
      check(all(c("Y_true", "Y_pred", "training_mean", "true_class", "predicted_class") %in%
        names(pred)), "Missing held-out prediction columns.")
      levels <- as.character(result$settings$class_levels)
      check(length(levels) == 2L && all(is.finite(pred$Y_true)) &&
        all(pred$Y_true %in% c(0, 1)) && all(is.finite(pred$Y_pred)) &&
        all(is.finite(pred$training_mean)) &&
        all(as.character(pred$true_class) %in% levels) &&
        all(pred$Y_true == as.integer(as.character(pred$true_class) == levels[2])),
        "Invalid binary outcomes in held-out predictions.")
      baseline_ss <- sum((pred$Y_true - pred$training_mean)^2)
      check(baseline_ss > 0, "Q2Y reference denominator is zero.")
      q2 <- 1 - sum((pred$Y_true - pred$Y_pred)^2) / baseline_ss
      positive <- pred$Y_true == 1
      np <- sum(positive)
      nn <- sum(!positive)
      check(np > 0 && nn > 0, "Both outcome classes are required.")
      auc <- (sum(rank(pred$Y_pred, ties.method = "average")[positive]) -
        np * (np + 1) / 2) / (np * nn)
      assigned <- ifelse(pred$Y_pred >= 0.5, levels[2], levels[1])
      check(all(as.character(pred$predicted_class) == assigned),
        "Saved predicted classes disagree with the 0.5 threshold.")
      ba <- (mean(assigned[positive] == levels[2]) +
        mean(assigned[!positive] == levels[1])) / 2
      check(max(abs(c(Q2Y = q2, AUROC = auc, balanced_accuracy = ba) - performance)) < 1e-10,
        "Predictive metrics disagree with held-out predictions.")
      if (!is.null(result$stability$jaccard)) {
        check(isTRUE(all.equal(jaccard, result$stability$jaccard, tolerance = 1e-12)),
          "Native Jaccard matrix mismatch.")
      } else cat(context, "native Jaccard unavailable; CSV verified from saved sets.\n")
      if (!is.null(result$stability$frequency)) {
        stored <- result$stability$frequency
        index <- match(features, stored$feature)
        check(all(frequency == stored$count[index]) &&
          max(abs(frequency / 5 - stored$frequency[index])) < 1e-12, "Native frequency mismatch.")
      }
      check(all(frequency %in% 0:5) && sum(frequency) == sum(lengths(sets)), "Frequency totals mismatch.")
      if (strategy == "V_intersect_M2") {
        tab <- tables$primary_selection_frequency
        tab <- tab[tab$dataset == dataset, ]
        check(nrow(tab) == length(features) && !anyDuplicated(tab$feature) &&
          setequal(tab$feature, features), "Frequency CSV probe identities mismatch.")
        index <- match(features, tab$feature)
        check(all(tab$count[index] == frequency) &&
          max(abs(tab$frequency[index] - frequency / 5)) < 1e-12, "Frequency CSV values mismatch.")
        distribution <- tables$primary_frequency_distribution
        distribution <- distribution[distribution$dataset == dataset, ]
        check(nrow(distribution) == 6L && setequal(distribution$folds_selected, 0:5) &&
          all(distribution$n_features[match(0:5, distribution$folds_selected)] == tabulate(frequency + 1L, 6L)),
          "Frequency-distribution CSV mismatch.")
        union_size <- length(Reduce(union, sets))
        intersection_size <- length(Reduce(intersect, sets))
        check(summary$union_features == union_size && summary$intersection_features == intersection_size,
          "Union/intersection summary mismatch.")
        expected <- if (dataset == "GSE68310") c(1847L, 629L) else c(6084L, 368L)
        check(identical(c(union_size, intersection_size), expected), "Paper union/intersection reference mismatch.")
      } else check(is.na(summary$union_features) && is.na(summary$intersection_features),
        "Union/intersection columns are reserved for the primary strategy.")
      reference_index <- reference_index + 1L
      check(median(lengths(sets)) == reference_medians[reference_index] &&
        abs(mean_jaccard - reference_jaccards[reference_index]) <= 5e-7 &&
        abs(performance["Q2Y"] - reference_q2[reference_index]) <= 5e-7 &&
        abs(performance["AUROC"] - reference_auc[reference_index]) <= 5e-7 &&
        abs(performance["balanced_accuracy"] - reference_ba[reference_index]) <= 5e-7,
        "Paper reference mismatch (metrics rounded to six decimals).")
      cat("PASS", context, "held-out metrics, stability, CSVs and paper references\n")
    }
  }
  require_valid(nrow(tables$primary_selection_frequency) == total_features, "Extra or missing frequency rows.")
  cat("PASS: all eight dataset-strategy combinations. No models fitted and no files written.\n")
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (!is.null(options)) verify_tables(options)
