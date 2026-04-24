#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(jsonlite)
  library(RhpcBLASctl)
})

blas_set_num_threads(1)
omp_set_num_threads(1)

parse_args <- function(args) {
  parsed <- list(
    gse = NULL,
    sample_id = NULL,
    output_root = NULL,
    project_root = ".",
    input_type = NULL,
    matrix_path = NULL,
    barcodes_path = NULL,
    features_path = NULL,
    h5_path = NULL,
    archive_path = NULL
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    value <- args[[i + 1]]
    name <- gsub("-", "_", substring(key, 3), fixed = TRUE)
    if (!name %in% names(parsed)) {
      stop("Unknown argument: ", key)
    }
    parsed[[name]] <- value
    i <- i + 2
  }

  required <- c("gse", "sample_id", "output_root", "input_type")
  missing <- required[vapply(required, function(item) {
    is.null(parsed[[item]]) || !nzchar(parsed[[item]])
  }, logical(1))]
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  parsed
}

resolve_data_root <- function(project_root) {
  candidates <- c(
    file.path(project_root, "data"),
    Sys.getenv("ML2026_DATA_ROOT", unset = ""),
    "/mnt/g/ML2026_data"
  )
  for (candidate in candidates) {
    if (nzchar(candidate) && dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Unable to resolve data root")
}

open_maybe_gz <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

read_lines_maybe_gz <- function(path) {
  con <- open_maybe_gz(path)
  on.exit(close(con), add = TRUE)
  readLines(con)
}

read_delim_maybe_gz <- function(path, sep = "\t") {
  con <- open_maybe_gz(path)
  on.exit(close(con), add = TRUE)
  read.delim(con, sep = sep, header = FALSE, stringsAsFactors = FALSE)
}

read_triplet_counts <- function(matrix_path, barcodes_path, features_path) {
  con <- open_maybe_gz(matrix_path)
  on.exit(close(con), add = TRUE)
  counts <- readMM(con)
  counts <- as(counts, "dgCMatrix")

  barcodes <- read_lines_maybe_gz(barcodes_path)
  features <- read_delim_maybe_gz(features_path, sep = "\t")

  feature_ids <- as.character(features[[1]])
  feature_names <- if (ncol(features) >= 2) {
    as.character(features[[2]])
  } else {
    feature_ids
  }

  if (ncol(features) >= 3) {
    feature_type <- as.character(features[[3]])
    if (any(feature_type == "Gene Expression")) {
      keep <- which(feature_type == "Gene Expression")
      counts <- counts[keep, , drop = FALSE]
      feature_ids <- feature_ids[keep]
      feature_names <- feature_names[keep]
    }
  }

  rownames(counts) <- make.unique(feature_names)
  colnames(counts) <- barcodes
  list(
    counts = counts,
    feature_ids = feature_ids,
    feature_names = feature_names
  )
}

read_archive_counts <- function(archive_path) {
  extract_dir <- tempfile(pattern = "ml2026_rna_archive_")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(archive_path, exdir = extract_dir)
  files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)

  matrix_path <- files[grepl("matrix\\.mtx$", files)][1]
  barcode_path <- files[grepl("barcodes\\.tsv$", files)][1]
  features_path <- files[grepl("(features|genes)\\.tsv$", files)][1]

  if (!nzchar(matrix_path) || !nzchar(barcode_path) || !nzchar(features_path)) {
    stop("Archive is missing matrix/barcodes/features triplet: ", archive_path)
  }

  read_triplet_counts(matrix_path, barcode_path, features_path)
}

read_h5_counts <- function(h5_path) {
  loaded <- Read10X_h5(h5_path, use.names = TRUE)
  counts <- if (is.list(loaded)) {
    if ("Gene Expression" %in% names(loaded)) {
      loaded[["Gene Expression"]]
    } else {
      loaded[[1]]
    }
  } else {
    loaded
  }

  counts <- as(counts, "dgCMatrix")
  list(
    counts = counts,
    feature_ids = rownames(counts),
    feature_names = rownames(counts)
  )
}

load_counts <- function(input_type, matrix_path, barcodes_path, features_path, h5_path, archive_path) {
  if (identical(input_type, "triplet")) {
    return(read_triplet_counts(matrix_path, barcodes_path, features_path))
  }
  if (identical(input_type, "h5")) {
    return(read_h5_counts(h5_path))
  }
  if (identical(input_type, "archive")) {
    return(read_archive_counts(archive_path))
  }
  stop("Unsupported RNA input_type: ", input_type)
}

median_mad_lower <- function(values, hard_floor) {
  values <- values[is.finite(values)]
  cutoff <- stats::median(values) - 4 * stats::mad(values)
  max(hard_floor, cutoff)
}

median_mad_upper <- function(values, hard_ceiling) {
  values <- values[is.finite(values)]
  cutoff <- stats::median(values) + 4 * stats::mad(values)
  min(hard_ceiling, cutoff)
}

normalize_rows <- function(mat) {
  norms <- sqrt(rowSums(mat * mat))
  norms[norms == 0] <- 1
  mat / norms
}

build_cima_annotation <- function(qc_obj, reference_dir) {
  feature_model <- read.delim(
    gzfile(file.path(reference_dir, "cima_rna_reference_pca_features.tsv.gz")),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  l1_centroids <- read.delim(
    file.path(reference_dir, "cima_rna_reference_l1_centroids.tsv"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  l2_centroids <- read.delim(
    file.path(reference_dir, "cima_rna_reference_l2_centroids.tsv"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hierarchy <- read.csv(
    file.path(reference_dir, "cima_rna_celltype_hierarchy.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  data_mat <- GetAssayData(qc_obj, assay = "RNA", layer = "data")
  feature_ids <- feature_model$feature_id
  query <- matrix(0, nrow = length(feature_ids), ncol = ncol(data_mat))
  rownames(query) <- feature_ids
  colnames(query) <- colnames(data_mat)

  matched <- intersect(feature_ids, rownames(data_mat))
  query[matched, ] <- as.matrix(data_mat[matched, , drop = FALSE])

  gene_mean <- feature_model$gene_mean
  gene_std <- feature_model$gene_std
  gene_std[gene_std == 0 | is.na(gene_std)] <- 1
  scaled <- sweep(query, 1, gene_mean, "-")
  scaled <- sweep(scaled, 1, gene_std, "/")

  loading_cols <- grep("^pc_dim_", colnames(feature_model), value = TRUE)
  loadings <- as.matrix(feature_model[, loading_cols, drop = FALSE])
  query_embeddings <- t(scaled) %*% loadings
  rownames(query_embeddings) <- colnames(data_mat)

  l1_embedding <- as.matrix(l1_centroids[, loading_cols, drop = FALSE])
  l2_embedding <- as.matrix(l2_centroids[, loading_cols, drop = FALSE])
  query_norm <- normalize_rows(query_embeddings)
  l1_norm <- normalize_rows(l1_embedding)
  l2_norm <- normalize_rows(l2_embedding)

  l1_similarity <- query_norm %*% t(l1_norm)
  l1_order <- t(apply(l1_similarity, 1, order, decreasing = TRUE))
  l1_top <- l1_order[, 1]
  l1_second <- if (ncol(l1_order) >= 2) l1_order[, 2] else l1_top

  l1_labels <- l1_centroids$cell_type_l1[l1_top]
  l1_scores <- l1_similarity[cbind(seq_len(nrow(l1_similarity)), l1_top)]
  l1_margins <- l1_similarity[cbind(seq_len(nrow(l1_similarity)), l1_top)] -
    l1_similarity[cbind(seq_len(nrow(l1_similarity)), l1_second)]

  l2_lookup <- split(hierarchy$cell_type_l2, hierarchy$cell_type_l1)
  l2_labels <- character(nrow(query_embeddings))
  l2_scores <- numeric(nrow(query_embeddings))
  l2_margins <- numeric(nrow(query_embeddings))

  for (i in seq_len(nrow(query_embeddings))) {
    allowed_l2 <- unique(l2_lookup[[l1_labels[i]]])
    keep_idx <- which(l2_centroids$cell_type_l2 %in% allowed_l2)
    if (length(keep_idx) == 0) {
      l2_labels[i] <- NA_character_
      l2_scores[i] <- NA_real_
      l2_margins[i] <- NA_real_
      next
    }
    sims <- drop(query_norm[i, , drop = FALSE] %*% t(l2_norm[keep_idx, , drop = FALSE]))
    ord <- order(sims, decreasing = TRUE)
    top <- ord[1]
    second <- if (length(ord) >= 2) ord[2] else ord[1]
    l2_labels[i] <- l2_centroids$cell_type_l2[keep_idx[top]]
    l2_scores[i] <- sims[top]
    l2_margins[i] <- sims[top] - sims[second]
  }

  low_conf <- (l1_scores < 0.35) | (l1_margins < 0.05)
  masked_l1 <- ifelse(low_conf, "Unknown", l1_labels)

  annotation <- data.frame(
    cell_barcode = rownames(query_embeddings),
    cima_cell_type_l1 = l1_labels,
    cima_cell_type_l2 = l2_labels,
    cima_cell_type_l1_masked = masked_l1,
    cima_l1_low_confidence = low_conf,
    cima_l1_score = l1_scores,
    cima_l1_score_margin = l1_margins,
    cima_l2_score = l2_scores,
    cima_l2_score_margin = l2_margins,
    stringsAsFactors = FALSE
  )

  annotation
}

compute_cluster_purity <- function(labels, clusters) {
  purity <- rep(NA_real_, length(labels))
  valid <- !is.na(labels) & nzchar(labels) & !is.na(clusters) & nzchar(clusters)
  for (cluster_id in unique(as.character(clusters[valid]))) {
    idx <- which(valid & as.character(clusters) == cluster_id)
    counts <- sort(table(labels[idx]), decreasing = TRUE)
    purity[idx] <- as.numeric(counts[1]) / length(idx)
  }
  purity
}

save_umap_plot <- function(df, label_col, out_path, title_text, unknown_gray = FALSE) {
  clean_df <- df[!is.na(df[[label_col]]) & nzchar(df[[label_col]]), , drop = FALSE]
  if (nrow(clean_df) == 0) {
    png(out_path, width = 1200, height = 900, res = 150)
    plot.new()
    title(main = title_text)
    text(0.5, 0.5, "No values available")
    dev.off()
    return(invisible(NULL))
  }

  labels <- sort(unique(as.character(clean_df[[label_col]])))
  palette <- grDevices::hcl.colors(length(labels), "Dynamic")
  names(palette) <- labels
  if (unknown_gray && "Unknown" %in% labels) {
    palette["Unknown"] <- "#7F7F7F"
  }

  clean_df[[label_col]] <- factor(clean_df[[label_col]], levels = labels)
  p <- ggplot(clean_df, aes(x = umap_rna_1, y = umap_rna_2, color = .data[[label_col]])) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_manual(values = palette, drop = FALSE) +
    labs(title = title_text, x = "UMAP_1", y = "UMAP_2", color = NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 8)
    )
  ggsave(filename = out_path, plot = p, width = 10, height = 8, units = "in", dpi = 150)
}

save_qc_overview <- function(metadata, out_path) {
  png(out_path, width = 1600, height = 1200, res = 150)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2))
  hist(metadata$nCount_RNA, breaks = 60, main = "nCount_RNA", xlab = "nCount_RNA", col = "#4C78A8")
  hist(metadata$nFeature_RNA, breaks = 60, main = "nFeature_RNA", xlab = "nFeature_RNA", col = "#F58518")
  hist(metadata$percent.mt, breaks = 60, main = "percent.mt", xlab = "percent.mt", col = "#54A24B")
  hist(metadata$percent.ribo, breaks = 60, main = "percent.ribo", xlab = "percent.ribo", col = "#E45756")
}

write_matrix_outputs <- function(counts, output_dir) {
  matrix_dir <- file.path(output_dir, "matrix")
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_path <- file.path(matrix_dir, "matrix.mtx")
  writeMM(counts, matrix_path)

  barcode_con <- gzfile(file.path(matrix_dir, "barcodes.tsv.gz"), open = "wt")
  writeLines(colnames(counts), barcode_con)
  close(barcode_con)

  feature_con <- gzfile(file.path(matrix_dir, "features.tsv.gz"), open = "wt")
  writeLines(rownames(counts), feature_con)
  close(feature_con)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(args$project_root, winslash = "/", mustWork = TRUE)
data_root <- resolve_data_root(project_root)
reference_dir <- file.path(data_root, "reference", "cima")
output_root <- normalizePath(args$output_root, winslash = "/", mustWork = FALSE)
output_dir <- file.path(output_root, args$gse, args$sample_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_info <- load_counts(
  input_type = args$input_type,
  matrix_path = args$matrix_path,
  barcodes_path = args$barcodes_path,
  features_path = args$features_path,
  h5_path = args$h5_path,
  archive_path = args$archive_path
)

counts <- counts_info$counts
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
counts <- counts[, colSums(counts) > 0, drop = FALSE]

if (ncol(counts) < 50) {
  stop("RNA input has too few columns for single-cell processing: ", ncol(counts))
}

obj <- CreateSeuratObject(counts = counts, project = args$sample_id, min.cells = 0, min.features = 0)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")

sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
sce <- scDblFinder(sce)
obj$scDblFinder.class <- as.character(colData(sce)$scDblFinder.class)
if ("scDblFinder.score" %in% colnames(colData(sce))) {
  obj$scDblFinder.score <- as.numeric(colData(sce)$scDblFinder.score)
}

meta <- obj@meta.data
count_cutoff <- median_mad_lower(meta$nCount_RNA, hard_floor = 500)
feature_cutoff <- median_mad_lower(meta$nFeature_RNA, hard_floor = 300)
mt_cutoff <- median_mad_upper(meta$percent.mt, hard_ceiling = 20)
ribo_cutoff <- median_mad_upper(meta$percent.ribo, hard_ceiling = 60)

meta$below_count_floor <- meta$nCount_RNA < count_cutoff
meta$below_feature_floor <- meta$nFeature_RNA < feature_cutoff
meta$above_mt_floor <- meta$percent.mt > mt_cutoff
meta$above_ribo_floor <- meta$percent.ribo > ribo_cutoff
meta$pass_qc <- !meta$below_count_floor &
  !meta$below_feature_floor &
  !meta$above_mt_floor &
  !meta$above_ribo_floor &
  meta$scDblFinder.class != "doublet"
obj@meta.data <- meta

save_qc_overview(obj@meta.data, file.path(output_dir, "qc_overview.png"))

qc_obj <- subset(obj, subset = pass_qc)
if (ncol(qc_obj) < 50) {
  stop("Too few cells pass RNA QC: ", ncol(qc_obj))
}

qc_obj <- NormalizeData(qc_obj, verbose = FALSE)
qc_obj <- FindVariableFeatures(qc_obj, selection.method = "vst", nfeatures = min(3000, nrow(qc_obj)), verbose = FALSE)
qc_obj <- ScaleData(qc_obj, features = VariableFeatures(qc_obj), verbose = FALSE)
qc_obj <- RunPCA(qc_obj, features = VariableFeatures(qc_obj), npcs = min(50, ncol(qc_obj) - 1), verbose = FALSE)

available_dims <- ncol(Embeddings(qc_obj[["pca"]]))
dims_use <- seq_len(min(30, available_dims))
if (length(dims_use) < 2) {
  stop("Not enough PCA dimensions for RNA UMAP")
}

qc_obj <- FindNeighbors(qc_obj, reduction = "pca", dims = dims_use, verbose = FALSE)
qc_obj <- FindClusters(qc_obj, resolution = 0.5, verbose = FALSE)
qc_obj <- RunUMAP(qc_obj, reduction = "pca", dims = dims_use, reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)

annotation <- build_cima_annotation(qc_obj, reference_dir)
qc_meta <- qc_obj@meta.data
qc_meta$cell_barcode <- rownames(qc_meta)
qc_meta$seurat_clusters <- as.character(qc_obj$seurat_clusters)
umap_embeddings <- Embeddings(qc_obj[["umap"]])
qc_meta$umap_rna_1 <- umap_embeddings[, 1]
qc_meta$umap_rna_2 <- umap_embeddings[, 2]
qc_meta <- merge(qc_meta, annotation, by = "cell_barcode", all.x = TRUE, sort = FALSE)
rownames(qc_meta) <- qc_meta$cell_barcode
qc_meta$cima_l1_cluster_purity <- compute_cluster_purity(qc_meta$cima_cell_type_l1, qc_meta$seurat_clusters)
qc_meta$cima_l2_cluster_purity <- compute_cluster_purity(qc_meta$cima_cell_type_l2, qc_meta$seurat_clusters)

for (column in c(
  "seurat_clusters",
  "umap_rna_1",
  "umap_rna_2",
  "cima_cell_type_l1",
  "cima_cell_type_l2",
  "cima_cell_type_l1_masked",
  "cima_l1_low_confidence",
  "cima_l1_score",
  "cima_l1_score_margin",
  "cima_l2_score",
  "cima_l2_score_margin",
  "cima_l1_cluster_purity",
  "cima_l2_cluster_purity"
)) {
  fill_value <- qc_meta[[column]]
  if (is.logical(fill_value)) {
    obj@meta.data[, column] <- FALSE
  } else if (is.numeric(fill_value)) {
    obj@meta.data[, column] <- NA_real_
  } else {
    obj@meta.data[, column] <- NA_character_
  }
  obj@meta.data[rownames(qc_meta), column] <- fill_value
}

validation_result <- qc_meta[, c(
  "cell_barcode",
  "seurat_clusters",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.ribo",
  "scDblFinder.class",
  "cima_cell_type_l1",
  "cima_cell_type_l2",
  "cima_cell_type_l1_masked",
  "cima_l1_low_confidence",
  "cima_l1_cluster_purity",
  "cima_l2_cluster_purity",
  "cima_l2_score",
  "cima_l2_score_margin",
  "umap_rna_1",
  "umap_rna_2"
), drop = FALSE]

save_umap_plot(validation_result, "seurat_clusters", file.path(output_dir, "umap_rna_clusters.png"),
               paste0(args$sample_id, " UMAP by RNA cluster"))
save_umap_plot(validation_result, "cima_cell_type_l1", file.path(output_dir, "umap_rna_cima_cell_type_l1.png"),
               paste0(args$sample_id, " UMAP by CIMA RNA L1"))
save_umap_plot(validation_result, "cima_cell_type_l2", file.path(output_dir, "umap_rna_cima_cell_type_l2.png"),
               paste0(args$sample_id, " UMAP by CIMA RNA L2"))
save_umap_plot(validation_result, "cima_cell_type_l1_masked", file.path(output_dir, "umap_rna_cima_cell_type_l1_masked.png"),
               paste0(args$sample_id, " UMAP by CIMA RNA L1 (masked)"), unknown_gray = TRUE)

write.csv(cbind(cell_barcode = rownames(obj@meta.data), obj@meta.data),
          file.path(output_dir, "metadata.csv"), row.names = FALSE)
write.csv(validation_result, file.path(output_dir, "validation_result.csv"), row.names = FALSE)
write.csv(qc_meta, file.path(output_dir, "metadata_qc.csv"), row.names = FALSE)

qc_summary <- data.frame(
  sample_id = args$sample_id,
  gse = args$gse,
  input_cells = ncol(obj),
  pass_qc = ncol(qc_obj),
  qc_rate = sprintf("%.2f%%", ncol(qc_obj) / ncol(obj) * 100),
  median_nCount_RNA = median(qc_meta$nCount_RNA),
  median_nFeature_RNA = median(qc_meta$nFeature_RNA),
  median_percent_mt = median(qc_meta$percent.mt),
  median_percent_ribo = median(qc_meta$percent.ribo),
  low_conf_cell_frac_l1 = sprintf("%.2f%%", mean(qc_meta$cima_l1_low_confidence) * 100),
  l1_unique_labels = length(unique(qc_meta$cima_cell_type_l1)),
  l2_unique_labels = length(unique(qc_meta$cima_cell_type_l2)),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)

write_matrix_outputs(GetAssayData(qc_obj, assay = "RNA", layer = "counts"), output_dir)
saveRDS(qc_obj, file.path(output_dir, paste0(args$sample_id, "_seurat_qc.rds")))
