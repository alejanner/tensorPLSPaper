# Plot the three independent datasets and pooled reconstruction RMSE.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 02_plot_synthetic_reconstruction.R FULL_BENCHMARK_DIR NEW_FIGURES_DIR"
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0) }
if (length(args) != 2L) stop(usage, call. = FALSE)

run_dir <- normalizePath(path.expand(args[1L]), mustWork = TRUE)
out_dir <- path.expand(args[2L])
if (file.exists(out_dir)) stop("Output directory already exists: ", out_dir, call. = FALSE)
output_parent <- normalizePath(dirname(out_dir), mustWork = FALSE)
if (identical(output_parent, run_dir) ||
    startsWith(output_parent, paste0(run_dir, .Platform$file.sep))) {
  stop("Choose a figures directory outside the benchmark results.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.", call. = FALSE)

input_files <- file.path(run_dir, c("paired_RMSE_by_replicate.csv", "paired_RMSE_by_level.csv"))
if (!all(file.exists(input_files))) {
  stop("Missing input: ", paste(input_files[!file.exists(input_files)], collapse = ", "), call. = FALSE)
}
by_dataset <- read.csv(input_files[1L], stringsAsFactors = FALSE)
pooled <- read.csv(input_files[2L], stringsAsFactors = FALSE)

required <- c("scenario", "level", "cases", "paired_cases", "failed_tucker",
              "mean_RMSE_noisy", "tucker_RMSE_noisy")
if (!all(c(required, "replicate") %in% names(by_dataset)) ||
    !all(required %in% names(pooled))) {
  stop("The benchmark tables do not have the expected columns.", call. = FALSE)
}

scenarios <- c("A", "B", "C")
levels <- c(.10, .20, .30, .40)
expected_datasets <- expand.grid(scenario = scenarios, replicate = 1:3,
                                 level = levels, stringsAsFactors = FALSE)
expected_pooled <- expand.grid(scenario = scenarios, level = levels,
                               stringsAsFactors = FALSE)
key_dataset <- function(x) paste(x$scenario, x$replicate, round(x$level * 100), sep = "|")
key_pooled <- function(x) paste(x$scenario, round(x$level * 100), sep = "|")

if (nrow(by_dataset) != nrow(expected_datasets) ||
    anyDuplicated(key_dataset(by_dataset)) ||
    !setequal(key_dataset(by_dataset), key_dataset(expected_datasets)) ||
    nrow(pooled) != nrow(expected_pooled) ||
    anyDuplicated(key_pooled(pooled)) ||
    !setequal(key_pooled(pooled), key_pooled(expected_pooled))) {
  stop("A complete full run with three datasets per scenario is required.", call. = FALSE)
}

for (tab in list(by_dataset, pooled)) {
  values <- as.matrix(tab[c("cases", "paired_cases", "failed_tucker",
                            "mean_RMSE_noisy", "tucker_RMSE_noisy")])
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(tab$paired_cases <= 0) || any(tab$cases != tab$paired_cases + tab$failed_tucker) ||
      any(tab$mean_RMSE_noisy < 0) || any(tab$tucker_RMSE_noisy < 0)) {
    stop("Non-finite RMSE values or invalid paired-case counts.", call. = FALSE)
  }
}

counts <- aggregate(cbind(cases, paired_cases, failed_tucker) ~ scenario + level,
                    data = by_dataset, FUN = sum)
check <- merge(counts, pooled, by = c("scenario", "level"),
               suffixes = c("_datasets", "_pooled"), sort = FALSE)
if (nrow(check) != nrow(pooled) ||
    !all(check$cases_datasets == check$cases_pooled) ||
    !all(check$paired_cases_datasets == check$paired_cases_pooled) ||
    !all(check$failed_tucker_datasets == check$failed_tucker_pooled)) {
  stop("Pooled and dataset-level paired-case counts disagree.", call. = FALSE)
}

make_long <- function(tab, dataset = FALSE) {
  keep <- data.frame(
    scenario = tab$scenario, level = as.numeric(tab$level) * 100,
    replicate = if (dataset) factor(tab$replicate, levels = 1:3) else NA,
    stringsAsFactors = FALSE
  )
  rbind(
    transform(keep, method = "Mean", RMSE = tab$mean_RMSE_noisy),
    transform(keep, method = "Tucker3", RMSE = tab$tucker_RMSE_noisy)
  )
}

individual <- make_long(by_dataset, dataset = TRUE)
summary <- make_long(pooled)
individual$scenario <- factor(individual$scenario, levels = scenarios)
summary$scenario <- factor(summary$scenario, levels = scenarios)
individual$method <- factor(individual$method, levels = c("Mean", "Tucker3"))
summary$method <- factor(summary$method, levels = c("Mean", "Tucker3"))

plot <- ggplot2::ggplot() +
  ggplot2::geom_line(data = individual,
    ggplot2::aes(x = level, y = RMSE, colour = method,
                 linetype = replicate, group = interaction(method, replicate)),
    linewidth = 0.8, alpha = 0.67) +
  ggplot2::geom_point(data = individual,
    ggplot2::aes(x = level, y = RMSE, colour = method),
    size = 1.25, alpha = 0.67) +
  ggplot2::geom_line(data = summary,
    ggplot2::aes(x = level, y = RMSE, colour = method, group = method),
    linewidth = 1.5) +
  ggplot2::geom_point(data = summary,
    ggplot2::aes(x = level, y = RMSE, colour = method),
    size = 2.3) +
  ggplot2::facet_wrap(~scenario, nrow = 1L, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = c(10, 20, 30, 40)) +
  ggplot2::scale_colour_manual(values = c(Mean = "#E56F68", Tucker3 = "#00AEBB")) +
  ggplot2::labs(x = "Masked visits (%)", y = "Pooled RMSE",
                colour = "Method", linetype = "Dataset",
                caption = "Thin lines: independent simulated datasets (1-3). Thick lines: pooled RMSE across matched valid cases. Panel y-scales differ.") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
    strip.text = ggplot2::element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.caption = ggplot2::element_text(size = 8, hjust = 0),
    panel.spacing = grid::unit(0.7, "lines")
  )

dir.create(out_dir, recursive = TRUE)
figure_name <- "Figure_S2_pooled_and_three_datasets"
ggplot2::ggsave(file.path(out_dir, paste0(figure_name, ".pdf")), plot,
                width = 11, height = 5, units = "in")
ggplot2::ggsave(file.path(out_dir, paste0(figure_name, ".png")), plot,
                width = 11, height = 5, units = "in", dpi = 300)
write.csv(individual, file.path(out_dir, "Figure_S2_independent_datasets.csv"), row.names = FALSE)

write.csv(
  summary[, c("scenario", "level", "method", "RMSE")],
  file.path(out_dir, "Figure_S2_pooled_RMSE.csv"),
  row.names = FALSE
)

write.csv(data.frame(file = basename(input_files), md5 = unname(tools::md5sum(input_files))),
          file.path(out_dir, "input_md5.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("Figure S2 saved in:", out_dir, "\n")
cat("PASS:", length(unique(key_dataset(by_dataset))),
    "scenario-dataset-level results and", nrow(pooled), "pooled results.\n")
