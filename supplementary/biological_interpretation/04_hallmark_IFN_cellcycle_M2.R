# Plot temporal M2 profiles for interferon and cell-cycle Hallmark groups.
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript 04_hallmark_IFN_cellcycle_M2.R",
  "MAPPING_RESULTS ENRICHMENT_RESULTS NEW_OUTPUT_DIR [--preflight]"
)
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0) }
preflight <- "--preflight" %in% args
args <- args[args != "--preflight"]
if (length(args) != 3L) stop(usage, call. = FALSE)

mapping_dir <- normalizePath(path.expand(args[[1L]]), mustWork = TRUE)
enrichment_dir <- normalizePath(path.expand(args[[2L]]), mustWork = TRUE)
output_path <- path.expand(args[[3L]])
if (file.exists(output_path)) stop("Output directory already exists: ", output_path, call. = FALSE)
output_parent <- normalizePath(dirname(output_path), mustWork = TRUE)
out_dir <- file.path(output_parent, basename(output_path))
inside <- function(path, directory) {
  identical(path, directory) || startsWith(path, paste0(directory, .Platform$file.sep))
}
if (inside(out_dir, mapping_dir) || inside(out_dir, enrichment_dir)) {
  stop("Output must be outside both input directories.", call. = FALSE)
}

vip_path <- file.path(mapping_dir, "results/GSE68310/selected_probe_M2_temporal_profiles.csv")
hallmark_dir <- file.path(enrichment_dir, "results/GSE68310/Hallmark")
overlap_path <- file.path(hallmark_dir, "Hallmark_significant_gene_overlaps.csv")
sig_path <- file.path(hallmark_dir, "Hallmark_significant_FDR005.csv")
inputs <- c(vip_path, overlap_path, sig_path)
if (!all(file.exists(inputs))) {
  stop("Missing input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "), call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
if (preflight) { cat("Preflight passed; no outputs created.\n"); quit(status = 0) }

vip <- read.csv(vip_path, colClasses = c(probe_id = "character", gene_id = "character"),
                check.names = FALSE)
overlap <- read.csv(overlap_path, colClasses = "character", check.names = FALSE)
sig <- read.csv(sig_path, colClasses = c(geneID = "character"), check.names = FALSE)

terms <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
times <- c("Day0", "Day2", "Day4", "Day6")
days <- c(0, 2, 4, 6)
stopifnot(all(c("probe_id", "gene_id", "gene_symbol", times) %in% names(vip)),
          all(c("ID", "gene_id") %in% names(overlap)),
          all(c("ID", "geneID", "Count", "p.adjust") %in% names(sig)),
          nrow(vip) == 1132L, !anyDuplicated(vip$probe_id),
          !anyDuplicated(sig$ID), all(terms %in% sig$ID),
          !anyNA(overlap$ID), !anyNA(overlap$gene_id),
          !anyDuplicated(overlap[c("ID", "gene_id")]),
          all(vapply(vip[times], is.numeric, logical(1))),
          all(is.finite(as.matrix(vip[times]))))

chosen <- sig[match(terms, sig$ID), , drop = FALSE]
if (anyNA(chosen$p.adjust) || any(chosen$p.adjust >= 0.05)) {
  stop("Some requested Hallmark terms are not significant", call. = FALSE)
}

# Recover the selected-gene overlaps reported by Hallmark ORA.
sets <- lapply(terms, function(id) unique(overlap$gene_id[overlap$ID == id]))
names(sets) <- terms
for (i in seq_along(terms)) {
  ids_in_result <- strsplit(chosen$geneID[[i]], "/", fixed = TRUE)[[1L]]
  stopifnot(length(sets[[i]]) == chosen$Count[[i]],
            setequal(sets[[i]], ids_in_result),
            all(sets[[i]] %in% vip$gene_id))
}

ifn <- union(sets[[1L]], sets[[2L]])
cell_cycle <- union(sets[[3L]], sets[[4L]])
shared <- intersect(ifn, cell_cycle)
ifn_only <- setdiff(ifn, cell_cycle)
cycle_only <- setdiff(cell_cycle, ifn)

# Assign every selected probe to exactly one group.
group <- rep("Other selected", nrow(vip))
group[vip$gene_id %in% shared] <- "Shared (not plotted)"
group[vip$gene_id %in% ifn_only] <- "Interferon only"
group[vip$gene_id %in% cycle_only] <- "Cell cycle only"
plot_groups <- c("Interferon only", "Cell cycle only", "Other selected")
if (!all(plot_groups %in% group)) stop("At least one comparison group is empty", call. = FALSE)
stopifnot(sum(table(group)) == nrow(vip),
          !any(group[vip$gene_id %in% shared] %in% plot_groups))

# Summarize probe-level M2 values at each time.
summarise <- function(g) {
  i <- group == g
  x <- vip[i, times, drop = FALSE]
  data.frame(
    group = g, time = times, day = days,
    n_probes = sum(i),
    n_genes = length(unique(vip$gene_id[i & !is.na(vip$gene_id) & nzchar(vip$gene_id)])),
    median_M2 = vapply(x, median, numeric(1)),
    mean_M2 = vapply(x, mean, numeric(1)),
    q25_M2 = vapply(x, quantile, numeric(1), probs = 0.25, names = FALSE),
    q75_M2 = vapply(x, quantile, numeric(1), probs = 0.75, names = FALSE)
  )
}
summary <- do.call(rbind, lapply(c(plot_groups, "Shared (not plotted)"), summarise))
rownames(summary) <- NULL

curves <- summary[summary$group %in% plot_groups[1:2], , drop = FALSE]
reference <- summary[summary$group == "Other selected", , drop = FALSE]
p <- ggplot2::ggplot(curves,
                     ggplot2::aes(x = day, y = median_M2, colour = group, fill = group)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_M2, ymax = q75_M2),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_line(data = reference,
                     ggplot2::aes(x = day, y = median_M2,
                                  linetype = "Other selected"),
                     inherit.aes = FALSE, linewidth = 0.9) +
  ggplot2::scale_linetype_manual(
    name = NULL, values = c("Other selected" = "dashed")
  ) +
  ggplot2::geom_point(data = reference,
                      ggplot2::aes(x = day, y = median_M2),
                      inherit.aes = FALSE, size = 1.5) +
  ggplot2::scale_x_continuous(breaks = days, labels = times) +
  ggplot2::labs(x = "Sampling time", y = "Temporal VIP", colour = NULL, fill = NULL) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 11),
    axis.text = ggplot2::element_text(size = 11),
    axis.title = ggplot2::element_text(size = 12),
    axis.line = ggplot2::element_line(linewidth = 0.35),
    axis.ticks = ggplot2::element_line(linewidth = 0.35),
    plot.margin = ggplot2::margin(2, 3, 2, 2, "mm")
  )

# Save results separately from the original analysis.
dir.create(out_dir, recursive = TRUE)
write.csv(chosen[, c("ID", "Count", "p.adjust")],
          file.path(out_dir, "Hallmark_terms.csv"), row.names = FALSE)
all_genes <- union(ifn, cell_cycle)
write.csv(data.frame(gene_id = all_genes,
                     group = ifelse(all_genes %in% shared, "Shared",
                                    ifelse(all_genes %in% ifn_only,
                                           "Interferon only", "Cell cycle only"))),
          file.path(out_dir, "gene_groups.csv"), row.names = FALSE)
write.csv(data.frame(vip[, c("probe_id", "gene_id", "gene_symbol", times), drop = FALSE], group),
          file.path(out_dir, "probe_groups.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "temporal_summary.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(out_dir, "M2_IFN_vs_cellcycle.pdf"), plot = p,
                width = 7, height = 3.2, units = "in")
ggplot2::ggsave(file.path(out_dir, "M2_IFN_vs_cellcycle.png"), plot = p,
                width = 7, height = 3.2, units = "in", dpi = 300)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("Interferon genes:", length(ifn), "| Cell-cycle genes:", length(cell_cycle),
    "| Shared genes:", length(shared), "(", paste(shared, collapse = ", "), ")\n")
cat("Exclusive IFN genes:", length(ifn_only),
    "| Exclusive cell-cycle genes:", length(cycle_only), "\n\n")
print(unique(summary[, c("group", "n_probes", "n_genes")]), row.names = FALSE)
cat("\nMedian M2 by group and time:\n")
print(summary[, c("group", "time", "median_M2")], row.names = FALSE)
cat("\nResults saved in:", out_dir, "\n")
