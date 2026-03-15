#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

make_placeholder_plot <- function(title, message) {
  ggplot() +
    theme_void() +
    ggtitle(title) +
    annotate("text", x = 0.5, y = 0.5, label = message, size = 5) +
    xlim(0, 1) +
    ylim(0, 1)
}

safe_plot <- function(title, expr) {
  plot_expr <- substitute(expr)
  tryCatch(
    eval(plot_expr, envir = parent.frame()),
    error = function(e) {
      message("Plot failed for ", title, ": ", conditionMessage(e))
      make_placeholder_plot(title, conditionMessage(e))
    }
  )
}

option_list <- list(
  make_option(c("--sample-dir"), type = "character", dest = "sample_dir",
              help = "Sample output directory containing *_seurat_qc.rds"),
  make_option(c("--output-file"), type = "character", dest = "output_file",
              default = "", help = "Override qc_overview output path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$sample_dir)) {
  print_help(opt_parser)
  stop("--sample-dir is required")
}

sample_dir <- normalizePath(opt$sample_dir, winslash = "/", mustWork = TRUE)
sample_name <- basename(sample_dir)
rds_candidates <- list.files(sample_dir, pattern = "_seurat_qc\\.rds$", full.names = TRUE)
if (length(rds_candidates) != 1) {
  stop("Expected exactly one *_seurat_qc.rds under ", sample_dir)
}
rds_file <- rds_candidates[1]
output_file <- if (nzchar(opt$output_file)) opt$output_file else file.path(sample_dir, "qc_overview.png")

atac_obj <- readRDS(rds_file)
meta_data <- atac_obj@meta.data
qc_cells <- rownames(meta_data)[meta_data$pass_qc %in% TRUE]

p1 <- safe_plot(
  "QC Metrics",
  VlnPlot(
    atac_obj,
    features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "FRiP",
                 "unique_ratio", "blacklist_fraction", "nucleosome_signal"),
    pt.size = 0,
    ncol = 3
  )
)

p2 <- safe_plot(
  "Pass QC Comparison",
  VlnPlot(
    atac_obj,
    features = c("nCount_ATAC", "TSS.enrichment", "FRiP",
                 "unique_ratio", "blacklist_fraction", "nucleosome_signal"),
    group.by = "pass_qc",
    pt.size = 0,
    ncol = 2
  )
)

p3 <- safe_plot(
  "Count vs TSS",
  FeatureScatter(atac_obj, feature1 = "nCount_ATAC", feature2 = "TSS.enrichment")
)

p4 <- safe_plot(
  "Count vs FRiP",
  FeatureScatter(atac_obj, feature1 = "nCount_ATAC", feature2 = "FRiP")
)

p5 <- safe_plot(
  "TSS vs Nucleosome",
  FeatureScatter(atac_obj, feature1 = "TSS.enrichment", feature2 = "nucleosome_signal")
)

p6 <- safe_plot(
  "TSS Profile",
  TSSPlot(atac_obj, group.by = "pass_qc")
)

p7 <- safe_plot(
  "Fragment Histogram",
  FragmentHistogram(atac_obj, group.by = "pass_qc")
)

p8 <- safe_plot(
  "Doublet Comparison",
  VlnPlot(
    atac_obj,
    features = c("nCount_ATAC", "nFeature_ATAC"),
    group.by = "scDblFinder.class",
    pt.size = 0,
    ncol = 2
  )
)

qc_overview <- wrap_plots(
  list(p1, p2, p3, p4, p5, p6, p7, p8),
  ncol = 2
) + plot_annotation(
  title = paste0(sample_name, " QC Overview"),
  subtitle = paste0(
    meta_data$dataset[1],
    " | initial cells: ", ncol(atac_obj),
    " | pass QC: ", length(qc_cells),
    " | doublets: ", sum(atac_obj$scDblFinder.class == "doublet", na.rm = TRUE)
  )
)

ggsave(
  filename = output_file,
  plot = qc_overview,
  width = 24,
  height = 32,
  units = "in",
  dpi = 150,
  limitsize = FALSE
)

cat("QC overview:", output_file, "\n")
