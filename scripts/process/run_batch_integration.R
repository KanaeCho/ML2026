#!/usr/bin/env Rscript

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, args, value = TRUE)
  if (length(matches) == 0) {
    stop("Unable to determine script path from commandArgs()")
  }
  normalizePath(sub(file_arg, "", matches[1]), winslash = "/", mustWork = TRUE)
}

project_root <- normalizePath(
  file.path(dirname(get_script_path()), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

project_r_lib <- if (.Platform$OS.type == "windows") {
  file.path(project_root, ".r-win-library")
} else {
  file.path(project_root, ".r-linux-library")
}
if (dir.exists(project_r_lib)) {
  .libPaths(c(project_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(jsonlite)
})

read_lines_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

write_json <- function(path, payload) {
  writeLines(toJSON(payload, auto_unbox = TRUE, pretty = TRUE), con = path)
}

parse_harmony_vars <- function(value) {
  vars <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  vars[nzchar(vars)]
}

peaks_to_granges <- function(peaks) {
  pieces <- tstrsplit(peaks, "-", fixed = TRUE, keep = 1:3)
  if (length(pieces) != 3) {
    stop("Unexpected peak format in features.tsv.gz")
  }
  GenomicRanges::GRanges(
    seqnames = pieces[[1]],
    ranges = IRanges::IRanges(start = as.integer(pieces[[2]]), end = as.integer(pieces[[3]]))
  )
}

save_plot <- function(path, plot, width = 10, height = 8) {
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

category_palette <- function(values) {
  values <- unique(as.character(values))
  base_colors <- c(
    "#0f766e", "#b45309", "#2563eb", "#be123c", "#7c3aed", "#0891b2", "#65a30d", "#c2410c",
    "#4f46e5", "#15803d", "#dc2626", "#1d4ed8", "#a21caf", "#0369a1", "#ca8a04", "#0f766e",
    "#6d28d9", "#0e7490", "#3f6212", "#9a3412", "#1e40af", "#9f1239", "#5b21b6", "#166534",
    "#0c4a6e", "#854d0e", "#374151", "#111827", "#9333ea", "#0f766e", "#22c55e", "#f97316"
  )
  colors <- rep(base_colors, length.out = length(values))
  names(colors) <- values
  colors
}

make_dim_plot <- function(object, reduction, group_by, title_text, point_size = 0.1) {
  values <- object[[group_by]][, 1]
  palette <- category_palette(values)
  DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = point_size,
    cols = palette,
    raster = FALSE
  ) +
    ggtitle(title_text) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 7)
    )
}

make_feature_plot <- function(object, reduction, features, title_text) {
  FeaturePlot(
    object = object,
    reduction = reduction,
    features = features,
    ncol = 2,
    pt.size = 0.1,
    raster = FALSE,
    cols = c("#f8fafc", "#b45309")
  ) &
    theme_bw(base_size = 12) &
    theme(plot.title = element_text(face = "bold"))
}

compute_row_nnz <- function(counts) {
  tabulate(counts@i + 1L, nbins = nrow(counts))
}

select_features <- function(counts, feature_names, min_peak_cells, max_features) {
  row_nnz <- compute_row_nnz(counts)
  keep <- which(row_nnz >= min_peak_cells)
  if (length(keep) == 0) {
    stop("No peaks left after applying min_peak_cells filter")
  }
  if (length(keep) > max_features) {
    ord <- order(row_nnz[keep], decreasing = TRUE)
    keep <- keep[ord[seq_len(max_features)]]
  }
  keep <- sort(keep)
  list(
    counts = counts[keep, , drop = FALSE],
    features = feature_names[keep],
    row_nnz = row_nnz[keep]
  )
}

load_integration_input <- function(input_dir, min_peak_cells, max_features) {
  matrix_file <- file.path(input_dir, "matrix.mtx")
  features_file <- file.path(input_dir, "features.tsv.gz")
  barcodes_file <- file.path(input_dir, "barcodes.tsv.gz")
  metadata_file <- file.path(input_dir, "merged_metadata.csv")

  if (!file.exists(matrix_file) || !file.exists(features_file) || !file.exists(barcodes_file) || !file.exists(metadata_file)) {
    stop("Missing integration input files under ", input_dir)
  }

  feature_names <- read_lines_gz(features_file)
  barcodes <- read_lines_gz(barcodes_file)
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"global_cell_id" %in% colnames(metadata)) {
    stop("global_cell_id column missing in ", metadata_file)
  }
  metadata <- metadata[match(barcodes, metadata$global_cell_id), , drop = FALSE]
  if (any(is.na(metadata$global_cell_id))) {
    stop("Metadata rows do not match barcodes.tsv.gz ordering")
  }

  counts <- readMM(matrix_file)
  counts <- as(counts, "dgTMatrix")
  selected <- select_features(counts, feature_names, min_peak_cells, max_features)
  rm(counts)
  invisible(gc())

  selected_counts <- as(selected$counts, "CsparseMatrix")
  rownames(selected_counts) <- selected$features
  colnames(selected_counts) <- barcodes

  list(
    counts = selected_counts,
    metadata = metadata,
    total_cells = length(barcodes),
    selected_features = selected$features,
    feature_nnz = selected$row_nnz,
    input_dir = input_dir
  )
}

build_annotation <- function() {
  annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"
  annotations
}

create_object <- function(input_data, project_name, annotations) {
  peak_ranges <- peaks_to_granges(rownames(input_data$counts))
  atac_assay <- CreateChromatinAssay(
    counts = input_data$counts,
    ranges = peak_ranges,
    genome = "hg38"
  )
  Annotation(atac_assay) <- annotations
  obj <- CreateSeuratObject(counts = atac_assay, assay = "ATAC", project = project_name, meta.data = input_data$metadata)
  obj
}

prepare_object <- function(object, npcs) {
  object <- RunTFIDF(object)
  object <- FindTopFeatures(object, min.cutoff = "q0")
  object <- RunSVD(object, n = npcs, reduction.key = "LSI_")
  object
}

run_harmony_integration <- function(object, harmony_vars, dims_use) {
  RunHarmony(
    object = object,
    group.by.vars = harmony_vars,
    reduction.use = "lsi",
    dims.use = dims_use,
    reduction.save = "harmony",
    project.dim = FALSE,
    plot_convergence = FALSE,
    verbose = FALSE
  )
}

compute_neighbor_batch_diversity <- function(embedding, labels, k = 30) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for batch mixing metric")
  }
  nn <- RANN::nn2(data = embedding, query = embedding, k = k + 1)
  neighbors <- nn$nn.idx[, -1, drop = FALSE]
  values <- numeric(nrow(neighbors))
  labels <- as.character(labels)
  unique_batches <- unique(labels)
  max_diversity <- if (length(unique_batches) > 1) 1 - 1 / length(unique_batches) else 1
  for (idx in seq_len(nrow(neighbors))) {
    neighbor_labels <- labels[neighbors[idx, ]]
    probs <- table(neighbor_labels) / length(neighbor_labels)
    simpson <- 1 - sum((as.numeric(probs))^2)
    values[idx] <- if (max_diversity > 0) simpson / max_diversity else 0
  }
  values
}

mixing_metric_safe <- function(object, grouping_var, reduction_name, dims_use, k = 30) {
  embedding <- Embeddings(object[[reduction_name]])[, dims_use, drop = FALSE]
  compute_neighbor_batch_diversity(embedding, object[[grouping_var]][, 1], k = k)
}

save_correlation_outputs <- function(object, output_dir, dims_use, qc_vars) {
  lsi_embeddings <- Embeddings(object[["lsi"]])[, dims_use, drop = FALSE]
  meta <- object@meta.data[, qc_vars, drop = FALSE]
  cor_mat <- matrix(NA_real_, nrow = ncol(lsi_embeddings), ncol = ncol(meta))
  rownames(cor_mat) <- colnames(lsi_embeddings)
  colnames(cor_mat) <- colnames(meta)
  for (i in seq_len(ncol(lsi_embeddings))) {
    for (j in seq_len(ncol(meta))) {
      values <- suppressWarnings(as.numeric(meta[[j]]))
      if (all(is.na(values))) {
        next
      }
      cor_mat[i, j] <- suppressWarnings(cor(lsi_embeddings[, i], values, use = "pairwise.complete.obs"))
    }
  }
  write.csv(as.data.frame(cor_mat), file.path(output_dir, "lsi_qc_correlations.csv"), quote = FALSE)
  png(file.path(output_dir, "lsi_qc_correlation_heatmap.png"), width = 1800, height = 1400, res = 180)
  pheatmap::pheatmap(
    cor_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#1d4ed8", "#f8fafc", "#b91c1c"))(101),
    breaks = seq(-1, 1, length.out = 102),
    main = "LSI dimension correlation with QC metrics"
  )
  dev.off()
}

save_cluster_heatmap <- function(object, output_dir, sample_var) {
  composition <- table(object$seurat_clusters, object[[sample_var]][, 1])
  composition_prop <- prop.table(composition, margin = 2)
  write.csv(as.data.frame.matrix(composition), file.path(output_dir, "cluster_by_sample_counts.csv"), quote = FALSE)
  write.csv(as.data.frame.matrix(composition_prop), file.path(output_dir, "cluster_by_sample_fraction.csv"), quote = FALSE)
  png(file.path(output_dir, "cluster_by_sample_heatmap.png"), width = 2200, height = 1400, res = 180)
  pheatmap::pheatmap(
    composition_prop,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#fff7ed", "#b45309"))(100),
    main = "Cluster by sample composition"
  )
  dev.off()
}

save_marker_summary <- function(object, output_dir, top_n = 10) {
  avg_access <- AverageExpression(object, assays = "ATAC", features = rownames(object), group.by = "seurat_clusters", verbose = FALSE)$ATAC
  nearest <- ClosestFeature(object, regions = rownames(avg_access))
  nearest <- nearest[, c("query_region", "gene_name", "distance")]
  rows <- list()
  for (cluster_name in colnames(avg_access)) {
    ord <- order(avg_access[, cluster_name], decreasing = TRUE)
    peaks <- rownames(avg_access)[ord[seq_len(min(top_n, length(ord)))]]
    cluster_rows <- data.frame(
      cluster = cluster_name,
      peak = peaks,
      average_accessibility = as.numeric(avg_access[peaks, cluster_name]),
      stringsAsFactors = FALSE
    )
    cluster_rows <- merge(cluster_rows, nearest, by.x = "peak", by.y = "query_region", all.x = TRUE, sort = FALSE)
    rows[[cluster_name]] <- cluster_rows
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(output_dir, "cluster_top_accessible_peaks.csv"), row.names = FALSE)
}

save_metadata_outputs <- function(object, output_dir, lsi_dims, harmony_dims) {
  meta <- object@meta.data
  lsi <- Embeddings(object[["lsi"]])[, lsi_dims, drop = FALSE]
  harmony_embed <- Embeddings(object[["harmony"]])[, harmony_dims, drop = FALSE]
  umap_lsi <- Embeddings(object[["umap_lsi"]])
  umap_harmony <- Embeddings(object[["umap_harmony"]])
  for (i in seq_len(ncol(lsi))) {
    meta[[paste0("LSI_", lsi_dims[i])]] <- lsi[, i]
  }
  for (i in seq_len(ncol(harmony_embed))) {
    meta[[paste0("Harmony_", harmony_dims[i])]] <- harmony_embed[, i]
  }
  meta$UMAP_LSI_1 <- umap_lsi[, 1]
  meta$UMAP_LSI_2 <- umap_lsi[, 2]
  meta$UMAP_Harmony_1 <- umap_harmony[, 1]
  meta$UMAP_Harmony_2 <- umap_harmony[, 2]
  fwrite(meta, file.path(output_dir, "integrated_metadata.csv.gz"))
}

write_report <- function(output_dir, summary_rows, mixing_rows, harmony_vars, dims_use) {
  lines <- c(
    "# Batch integration report",
    "",
    paste0("- Harmony variables: `", paste(harmony_vars, collapse = ", "), "`"),
    paste0("- Dimensions used: `", paste(range(dims_use), collapse = ":"), "`"),
    "",
    "## Summary",
    ""
  )
  for (row in summary_rows) {
    lines <- c(lines, paste0("- ", row))
  }
  lines <- c(lines, "", "## Batch mixing", "")
  for (row in mixing_rows) {
    lines <- c(lines, paste0("- ", row))
  }
  writeLines(lines, con = file.path(output_dir, "integration_report.md"))
}

option_list <- list(
  make_option("--input-dir", type = "character", dest = "input_dir"),
  make_option("--output-dir", type = "character", dest = "output_dir"),
  make_option("--label", type = "character", default = "integration_run", dest = "label"),
  make_option("--harmony-vars", type = "character", default = "source_gse,source_gsm", dest = "harmony_vars"),
  make_option("--min-peak-cells", type = "integer", default = 50L, dest = "min_peak_cells"),
  make_option("--max-features", type = "integer", default = 30000L, dest = "max_features"),
  make_option("--npcs", type = "integer", default = 30L, dest = "npcs"),
  make_option("--dims-start", type = "integer", default = 2L, dest = "dims_start"),
  make_option("--dims-end", type = "integer", default = 30L, dest = "dims_end"),
  make_option("--neighbors-k", type = "integer", default = 30L, dest = "neighbors_k"),
  make_option("--resolution", type = "double", default = 0.6, dest = "resolution"),
  make_option("--save-rds", action = "store_true", default = FALSE, dest = "save_rds")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  stop("Both --input-dir and --output-dir are required")
}

input_dir <- normalizePath(opt$input_dir, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

harmony_vars <- parse_harmony_vars(opt$harmony_vars)
if (length(harmony_vars) == 0) {
  stop("At least one harmony variable is required")
}

dims_use <- seq.int(opt$dims_start, opt$dims_end)
qc_vars <- c("nCount_ATAC", "FRiP", "TSS.enrichment", "blacklist_fraction", "nucleosome_signal")

cat("[1/8] Loading integration input...\n")
input_data <- load_integration_input(input_dir, opt$min_peak_cells, opt$max_features)

cat("[2/8] Building Seurat/Signac object...\n")
annotations <- build_annotation()
obj <- create_object(input_data, opt$label, annotations)

missing_harmony_vars <- harmony_vars[!harmony_vars %in% colnames(obj@meta.data)]
if (length(missing_harmony_vars) > 0) {
  stop("Missing harmony metadata columns: ", paste(missing_harmony_vars, collapse = ", "))
}

cat("[3/8] Running TF-IDF and LSI...\n")
obj <- prepare_object(obj, opt$npcs)

cat("[4/8] Building pre-Harmony UMAP...\n")
obj <- RunUMAP(obj, reduction = "lsi", dims = dims_use, reduction.name = "umap_lsi", reduction.key = "UMAPLSI_", verbose = FALSE)

cat("[5/8] Running Harmony...\n")
obj <- run_harmony_integration(obj, harmony_vars, dims_use)
harmony_dims <- seq_len(ncol(Embeddings(obj[["harmony"]])))

cat("[6/8] Running graph, clustering, and Harmony UMAP...\n")
obj <- FindNeighbors(obj, reduction = "harmony", dims = harmony_dims, verbose = FALSE)
obj <- FindClusters(obj, resolution = opt$resolution, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = harmony_dims, reduction.name = "umap_harmony", reduction.key = "UMAPHARM_", verbose = FALSE)

cat("[7/8] Computing integration quality metrics...\n")
mixing_gse <- mixing_metric_safe(obj, "source_gse", "harmony", harmony_dims, opt$neighbors_k)
mixing_gsm <- mixing_metric_safe(obj, "source_gsm", "harmony", harmony_dims, opt$neighbors_k)
obj$mixing_source_gse <- mixing_gse
obj$mixing_source_gsm <- mixing_gsm

mixing_summary <- data.frame(
  metric = c("mixing_source_gse", "mixing_source_gsm"),
  mean = c(mean(mixing_gse, na.rm = TRUE), mean(mixing_gsm, na.rm = TRUE)),
  median = c(median(mixing_gse, na.rm = TRUE), median(mixing_gsm, na.rm = TRUE)),
  q05 = c(quantile(mixing_gse, 0.05, na.rm = TRUE), quantile(mixing_gsm, 0.05, na.rm = TRUE)),
  q95 = c(quantile(mixing_gse, 0.95, na.rm = TRUE), quantile(mixing_gsm, 0.95, na.rm = TRUE))
)
write.csv(mixing_summary, file.path(output_dir, "batch_mixing_metrics.csv"), row.names = FALSE)

cat("[8/8] Saving outputs...\n")
save_plot(
  file.path(output_dir, "pre_harmony_lsi_by_gse.png"),
  make_dim_plot(obj, "umap_lsi", "source_gse", "Pre-Harmony LSI UMAP by GSE"),
  width = 9,
  height = 7
)
save_plot(
  file.path(output_dir, "pre_harmony_lsi_by_gsm.png"),
  make_dim_plot(obj, "umap_lsi", "source_gsm", "Pre-Harmony LSI UMAP by GSM"),
  width = 12,
  height = 9
)
save_plot(
  file.path(output_dir, "post_harmony_umap_by_gse.png"),
  make_dim_plot(obj, "umap_harmony", "source_gse", "Post-Harmony UMAP by GSE"),
  width = 9,
  height = 7
)
save_plot(
  file.path(output_dir, "post_harmony_umap_by_gsm.png"),
  make_dim_plot(obj, "umap_harmony", "source_gsm", "Post-Harmony UMAP by GSM"),
  width = 12,
  height = 9
)
save_plot(
  file.path(output_dir, "post_harmony_umap_by_cluster.png"),
  DimPlot(obj, reduction = "umap_harmony", group.by = "seurat_clusters", pt.size = 0.1, raster = FALSE) +
    ggtitle("Post-Harmony UMAP by cluster") + theme_bw(base_size = 12),
  width = 10,
  height = 8
)
save_plot(
  file.path(output_dir, "post_harmony_umap_by_qc.png"),
  make_feature_plot(obj, "umap_harmony", qc_vars, "Post-Harmony UMAP by QC metrics"),
  width = 12,
  height = 10
)

save_correlation_outputs(obj, output_dir, dims_use, qc_vars)
save_cluster_heatmap(obj, output_dir, "source_gsm")
save_marker_summary(obj, output_dir, top_n = 10)
save_metadata_outputs(obj, output_dir, dims_use, harmony_dims)

summary_payload <- list(
  label = opt$label,
  input_dir = input_dir,
  output_dir = output_dir,
  total_cells = ncol(obj),
  selected_features = nrow(obj),
  harmony_vars = harmony_vars,
  dims_use = dims_use,
  neighbors_k = opt$neighbors_k,
  resolution = opt$resolution,
  clusters = length(unique(obj$seurat_clusters)),
  mixing_source_gse_mean = mean(mixing_gse, na.rm = TRUE),
  mixing_source_gsm_mean = mean(mixing_gsm, na.rm = TRUE)
)
write_json(file.path(output_dir, "integration_summary.json"), summary_payload)

summary_rows <- c(
  paste0("Cells integrated: ", format(ncol(obj), big.mark = ",")),
  paste0("Peaks retained for integration: ", format(nrow(obj), big.mark = ",")),
  paste0("Clusters: ", length(unique(obj$seurat_clusters)))
)
mixing_rows <- c(
  paste0("source_gse mean mixing score: ", sprintf("%.4f", mean(mixing_gse, na.rm = TRUE))),
  paste0("source_gsm mean mixing score: ", sprintf("%.4f", mean(mixing_gsm, na.rm = TRUE)))
)
write_report(output_dir, summary_rows, mixing_rows, harmony_vars, dims_use)

if (isTRUE(opt$save_rds)) {
  saveRDS(obj, file.path(output_dir, paste0(opt$label, "_integration.rds")))
}

cat("Integration outputs:", output_dir, "\n")
