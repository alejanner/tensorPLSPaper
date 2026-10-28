# Plot GSE68310 pooled reconstruction RMSE with subject bootstrap intervals.
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 03_plot_gse68310_reconstruction.R BENCHMARK_DIR NEW_FIGURES_DIR"
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0) }
if (length(args) != 2L) stop(usage, call. = FALSE)

run_dir <- normalizePath(path.expand(args[1L]), mustWork = TRUE)
out_dir <- path.expand(args[2L])
if (file.exists(out_dir)) stop("Output directory already exists: ", out_dir, call. = FALSE)

# Keep figures separate from the benchmark results.
parent <- normalizePath(dirname(out_dir), mustWork = TRUE)
out_dir <- file.path(parent, basename(out_dir))
if (identical(out_dir, run_dir) ||
    startsWith(out_dir, paste0(run_dir, .Platform$file.sep))) {
  stop("Choose a figures directory outside the benchmark results.", call. = FALSE)
}

real_file <- file.path(run_dir, "subject_level_errors.csv")
summary_file <- file.path(run_dir, "pooled_RMSE.csv")
required <- c(real_file, summary_file)
if (!all(file.exists(required))) {
  stop("Missing benchmark input: ", paste(required[!file.exists(required)], collapse = ", "), call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.", call. = FALSE)
}

need_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
}
write_result <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE)
}
save_plot <- function(p, filename, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(filename, ".pdf")), p,
                  width = width, height = height, units = "in")
  ggplot2::ggsave(file.path(out_dir, paste0(filename, ".png")), p,
                  width = width, height = height, units = "in", dpi = 300)
}

dir.create(out_dir, recursive = TRUE)

# Figure S1: resample subjects within each held-out fold, retaining repeated masks.
real <- read.csv(real_file, stringsAsFactors = FALSE)
need_columns(real,
  c("subject_id", "fold", "seed", "level", "n", "SSE_mean", "SSE_tucker"),
  "Real-data file")
stopifnot(nrow(real) == 896L, length(unique(real$subject_id)) == 108L,
          length(unique(real$fold)) == 5L,
          identical(sort(unique(real$seed)), 123:125),
          isTRUE(all.equal(sort(unique(real$level)), c(.05, .10, .20, .30, .40))),
          all(is.finite(real$n)), all(real$n > 0),
          all(is.finite(real$SSE_mean)), all(real$SSE_mean >= 0),
          all(is.finite(real$SSE_tucker)), all(real$SSE_tucker >= 0))
subject_folds <- split(real$fold, real$subject_id)
stopifnot(all(vapply(subject_folds, function(x) length(unique(x)) == 1L, logical(1))))

subjects <- aggregate(
  cbind(n, SSE_mean, SSE_tucker) ~ level + fold + subject_id,
  data = real, FUN = sum
)
rmse_real <- function(z) {
  c(Mean = sqrt(sum(z$SSE_mean) / sum(z$n)),
    Tucker3 = sqrt(sum(z$SSE_tucker) / sum(z$n)))
}

set.seed(123L)
B <- 2000L
s1 <- do.call(rbind, lapply(sort(unique(subjects$level)), function(level) {
  z <- subjects[subjects$level == level, , drop = FALSE]
  within_fold <- split(z, z$fold)
  point <- rmse_real(z)
  boot <- replicate(B, {
    sampled <- do.call(rbind, lapply(within_fold, function(x) {
      x[sample.int(nrow(x), size = nrow(x), replace = TRUE), , drop = FALSE]
    }))
    rmse_real(sampled)
  })
  low <- apply(boot, 1, stats::quantile, probs = .025, names = FALSE)
  high <- apply(boot, 1, stats::quantile, probs = .975, names = FALSE)
  data.frame(level = level, method = names(point), pooled_RMSE = unname(point),
    lower_95 = unname(low), upper_95 = unname(high),
    n_subjects_evaluated = nrow(z), bootstrap_reps = B)
}))
rownames(s1) <- NULL

reference_real_file <- file.path(dirname(real_file), "pooled_RMSE.csv")
if (!file.exists(reference_real_file)) {
  stop("Missing original real-data summary: ", reference_real_file)
}
reference_real <- read.csv(reference_real_file, stringsAsFactors = FALSE)
need_columns(reference_real, c("level", "method", "pooled_RMSE"), "Real summary")
for (i in seq_len(nrow(s1))) {
  z <- s1[i, ]
  original_method <- if (z$method == "Mean") "training_mean" else "tucker"
  expected <- reference_real$pooled_RMSE[
    reference_real$method == original_method &
    abs(reference_real$level - z$level) < 1e-10
  ]
  stopifnot(length(expected) == 1L,
    abs(z$pooled_RMSE - expected) < 1e-9)
}
write_result(s1, "Figure_S1_pooled_RMSE_bootstrap.csv")

p1 <- ggplot2::ggplot(s1, ggplot2::aes(
  x = level * 100, y = pooled_RMSE, colour = method, group = method
)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_95, ymax = upper_95),
                         width = 1.3, linewidth = .45) +
  ggplot2::geom_line(linewidth = .85) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_continuous(breaks = c(5, 10, 20, 30, 40)) +
  ggplot2::labs(x = "Masked visits (%)", y = "Pooled RMSE", colour = "Method") +
  ggplot2::theme_classic(base_size = 11)
save_plot(p1, "Figure_S1_pooled_RMSE_bootstrap", 7, 4.5)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("Figure S1 saved in:", out_dir, "\n")
cat("PASS: Pooled RMSE agrees with the saved benchmark results.\n")
