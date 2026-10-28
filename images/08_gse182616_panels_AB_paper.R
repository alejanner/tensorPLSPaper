# Export GSE182616 paper panels A and B from the saved descriptive models.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 08_gse182616_panels_AB_paper.R <TensorPLS_source> <saved_run_dir> <output_dir>",
       call. = FALSE)
}

source_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
run_dir <- normalizePath(path.expand(args[[2L]]), mustWork = TRUE)
out_dir <- path.expand(args[[3L]])

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("pkgload and ggplot2 are required.", call. = FALSE)
}
pkgload::load_all(source_dir, export_all = FALSE, quiet = TRUE)

model_file <- file.path(run_dir, "GSE182616_analysis", "descriptive_models.rds")
selected_file <- file.path(run_dir, "GSE182616_analysis", "selected_probes.csv")
if (!file.exists(model_file) || !file.exists(selected_file)) {
  stop("Missing descriptive_models.rds or selected_probes.csv.", call. = FALSE)
}

saved <- readRDS(model_file)
selected <- read.csv(selected_file, colClasses = "character")$feature
model <- saved$model_selected
if (!inherits(model, "tensorpls_model") || ncol(model$x_scores) < 2L ||
    !identical(model$feature_names, selected) ||
    !identical(names(model$training_class), model$subject_ids)) {
  stop("Saved selected model does not match the selected probes or subject IDs.",
       call. = FALSE)
}
if (length(selected) != 2175L || model$ncomp != 2L || length(model$time_names) != 6L) {
  stop("Selected-feature count, component count or number of times differ from the reference.", call. = FALSE)
}

paper_theme <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(2, 3, 2, 2, "mm")
    )
}

plot_paper_scores <- function(model) {
  # Plot the native scores and classes saved in the fitted model.
  scores <- data.frame(
    Comp1 = model$x_scores[, 1L],
    Comp2 = model$x_scores[, 2L],
    class = unname(model$training_class)
  )
  palette <- stats::setNames(
    c("dodgerblue2", "orange2", "green3", "red3", "purple3", "gold2")[seq_len(nlevels(scores$class))],
    levels(scores$class)
  )
  ggplot2::ggplot(scores, ggplot2::aes(x = Comp1, y = Comp2)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(fill = class), shape = 21,
                        colour = "black", stroke = 0.3, size = 2.1) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Comp 1", y = "Comp 2", fill = NULL) +
    paper_theme()
}

plot_paper_temporal_vip <- function(model, features) {
  # The temporal profile is computed from the selected-feature refit.
  vips <- TensorPLS::tensorpls_vips(model)
  summary <- TensorPLS::temporal_vip_summary(
    vips, features = features, summary = "mean", dispersion = "IQR"
  )
  summary$time_index <- seq_len(nrow(summary))
  orange <- "#E69F00"
  plot <- ggplot2::ggplot(summary, ggplot2::aes(x = time_index, y = importance)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75),
                         fill = orange, alpha = 0.16) +
    ggplot2::geom_line(colour = orange, linewidth = 0.9) +
    ggplot2::geom_point(colour = orange, size = 2.8) +
    ggplot2::scale_x_continuous(breaks = summary$time_index,
                                 labels = summary$time,
                                 expand = ggplot2::expansion(mult = 0.06)) +
    ggplot2::labs(x = "Time point", y = "Temporal VIP") +
    paper_theme()
  list(plot = plot, summary = summary)
}

panels <- c("GSE182616_Scores_Paper.pdf", "GSE182616_Temporal_VIP_Paper.pdf",
            "GSE182616_Temporal_VIP_Summary.csv")
if (any(file.exists(file.path(out_dir, panels)))) {
  stop("Output files already exist. Choose a fresh output directory.", call. = FALSE)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scores_plot <- plot_paper_scores(model)
vip_plot <- plot_paper_temporal_vip(model, selected)

ggplot2::ggsave(file.path(out_dir, panels[[1L]]), scores_plot,
                width = 4.85, height = 3.2, units = "in")
ggplot2::ggsave(file.path(out_dir, panels[[2L]]), vip_plot$plot,
                width = 4.85, height = 3.2, units = "in")
write.csv(vip_plot$summary, file.path(out_dir, panels[[3L]]), row.names = FALSE)

cat("Saved panels A and B to:", normalizePath(out_dir), "\n")
cat("Selected probes:", length(selected), "\n")
print(vip_plot$summary[, c("time", "importance", "q25", "q75")], row.names = FALSE)
