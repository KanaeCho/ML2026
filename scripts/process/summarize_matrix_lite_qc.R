#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) == 0) {
  stop("Could not determine script path from commandArgs")
}
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

option_list <- list(
  make_option(c("--input-root"), type = "character", default = file.path(project_root, "output")),
  make_option(c("--output-dir"), type = "character", default = file.path(project_root, "output", "1.only_atac", "qc_reports")),
  make_option(c("--tag"), type = "character", default = "baseline")
)
opt <- parse_args(OptionParser(option_list = option_list))

is_bad <- function(x, direction = c("low", "high"), nmads = 3) {
  direction <- match.arg(direction)
  vals <- suppressWarnings(as.numeric(x))
  med <- stats::median(vals, na.rm = TRUE)
  mad_val <- stats::mad(vals, center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(mad_val) || mad_val == 0) {
    return(rep(FALSE, length(vals)))
  }
  if (direction == "low") {
    vals < (med - nmads * mad_val)
  } else {
    vals > (med + nmads * mad_val)
  }
}

load_summary <- function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE)
  out <- stats::setNames(tab$value, tab$metric)
  sample_dir <- dirname(path)
  c(
    gse = out[["gse"]],
    gsm = out[["gsm"]],
    input_cells = out[["input_cells"]],
    pass_qc = out[["pass_qc"]],
    qc_rate = gsub("%", "", out[["qc_rate"]]),
    singlets = out[["singlets"]],
    doublets = out[["doublets"]],
    mean_nCount_ATAC = out[["mean_nCount_ATAC"]],
    median_TSS_enrichment = out[["median_TSS_enrichment"]],
    median_FRiP = out[["median_FRiP"]],
    median_fragments = out[["median_fragments"]],
    median_unique_ratio = out[["median_unique_ratio"]],
    median_blacklist_fraction = out[["median_blacklist_fraction"]],
    median_cima_l4_score = out[["median_cima_l4_score"]],
    cima_unique_l4_labels = out[["cima_unique_l4_labels"]],
    cima_umap_basis = out[["cima_umap_basis"]],
    qc_summary_path = path,
    sample_dir = sample_dir
  )
}

input_root <- normalizePath(opt$`input-root`, winslash = "/", mustWork = TRUE)
output_dir <- opt$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_files <- list.files(input_root, pattern = "^qc_summary\\.csv$", recursive = TRUE, full.names = TRUE)
summary_files <- summary_files[grepl("/GSM[^/]+/qc_summary\\.csv$", summary_files)]
if (length(summary_files) == 0) {
  stop("No qc_summary.csv files found under ", input_root)
}

rows <- lapply(summary_files, load_summary)
df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)

numeric_cols <- c(
  "input_cells", "pass_qc", "qc_rate", "singlets", "doublets",
  "mean_nCount_ATAC", "median_TSS_enrichment", "median_FRiP", "median_fragments",
  "median_unique_ratio", "median_blacklist_fraction", "median_cima_l4_score",
  "cima_unique_l4_labels"
)
for (col in numeric_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
df$doublet_rate <- ifelse(df$input_cells > 0, 100 * df$doublets / df$input_cells, NA_real_)

df$flag_low_qc_rate <- is_bad(df$qc_rate, direction = "low")
df$flag_low_tss <- is_bad(df$median_TSS_enrichment, direction = "low")
df$flag_low_frip <- is_bad(df$median_FRiP, direction = "low")
df$flag_high_blacklist <- is_bad(df$median_blacklist_fraction, direction = "high")
df$flag_high_doublet <- is_bad(df$doublet_rate, direction = "high")
df$flag_low_cima_score <- is_bad(df$median_cima_l4_score, direction = "low")
df$flag_any <- Reduce(`|`, df[, grep("^flag_", colnames(df)), drop = FALSE])

write.csv(df, file.path(output_dir, paste0("matrix_lite_qc_summary_", opt$tag, ".csv")), row.names = FALSE)

flag_cols <- c(
  qc_rate = "flag_low_qc_rate",
  median_TSS_enrichment = "flag_low_tss",
  median_FRiP = "flag_low_frip",
  median_blacklist_fraction = "flag_high_blacklist",
  doublet_rate = "flag_high_doublet",
  median_cima_l4_score = "flag_low_cima_score"
)

plot_metric <- function(metric_col, flag_col, title_text, out_name) {
  plot_df <- df[order(df[[metric_col]], decreasing = FALSE), , drop = FALSE]
  plot_df$sample <- factor(paste(plot_df$gse, plot_df$gsm, sep = "/"), levels = paste(plot_df$gse, plot_df$gsm, sep = "/"))
  plot_df$flag_label <- ifelse(plot_df[[flag_col]], "flagged", "ok")
  p <- ggplot(plot_df, aes(x = sample, y = .data[[metric_col]], fill = flag_label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c(flagged = "#D7191C", ok = "#2C7BB6")) +
    labs(title = title_text, x = NULL, y = metric_col, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      plot.title = element_text(face = "bold")
    )
  ggsave(file.path(output_dir, out_name), p, width = 16, height = 6, units = "in", dpi = 150, limitsize = FALSE)
}

for (metric_col in names(flag_cols)) {
  plot_metric(
    metric_col = metric_col,
    flag_col = flag_cols[[metric_col]],
    title_text = paste("Matrix-lite sample QC:", metric_col),
    out_name = paste0("matrix_lite_", metric_col, "_", opt$tag, ".png")
  )
}

flagged <- df[df$flag_any, c("gse", "gsm", "qc_rate", "median_TSS_enrichment", "median_FRiP", "median_blacklist_fraction", "doublet_rate", "median_cima_l4_score", grep("^flag_", colnames(df), value = TRUE)), drop = FALSE]
write.csv(flagged, file.path(output_dir, paste0("matrix_lite_flagged_samples_", opt$tag, ".csv")), row.names = FALSE)

cat("Summary table:", file.path(output_dir, paste0("matrix_lite_qc_summary_", opt$tag, ".csv")), "\n")
cat("Flagged samples:", file.path(output_dir, paste0("matrix_lite_flagged_samples_", opt$tag, ".csv")), "\n")
