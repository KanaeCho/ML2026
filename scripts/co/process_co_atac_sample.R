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

suppressPackageStartupMessages({
  library(optparse)
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(Matrix)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})

find_single_file <- function(raw_dir, gsm, token, required = TRUE, exts = c("tsv")) {
  ext_pattern <- paste(exts, collapse = "|")
  pattern <- paste0("^", gsm, "_.*", token, ".*\\.(", ext_pattern, ")\\.gz$")
  files <- list.files(raw_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    if (required) {
      stop("No file found for ", gsm, " with token '", token, "' under ", raw_dir)
    }
    return(NULL)
  }
  if (length(files) > 1) {
    stop("Multiple files found for ", gsm, " with token '", token, "': ", paste(basename(files), collapse = ", "))
  }
  files[1]
}

write_lines_gz <- function(lines, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, sep = "\n")
}

write_lines_plain <- function(lines, out_path) {
  con <- file(out_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, useBytes = TRUE)
}

safe_plot <- function(title, expr) {
  plot_expr <- substitute(expr)
  tryCatch(
    eval(plot_expr, envir = parent.frame()),
    error = function(e) {
      ggplot() + theme_void() + ggtitle(title) +
        annotate("text", x = 0.5, y = 0.5, label = conditionMessage(e), size = 5) +
        xlim(0, 1) + ylim(0, 1)
    }
  )
}

ensure_tabix_index <- function(fragment_file) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    return(paste0(fragment_file, ".tbi"))
  }
  tbi_file <- paste0(fragment_file, ".tbi")
  if (!file.exists(tbi_file)) {
    Rsamtools::indexTabix(fragment_file, format = "bed")
  }
  tbi_file
}

infer_barcodes <- function(fragment_file, min_fragments, max_barcodes) {
  fragment_stats <- CountFragments(fragment_file)
  fragment_stats <- fragment_stats[!is.na(fragment_stats$CB) & nzchar(fragment_stats$CB), , drop = FALSE]
  fragment_stats <- fragment_stats[fragment_stats$frequency_count > 0, , drop = FALSE]
  fragment_stats <- fragment_stats[order(fragment_stats$frequency_count, decreasing = TRUE), , drop = FALSE]
  fragment_stats$rank <- seq_len(nrow(fragment_stats))
  keep <- fragment_stats$frequency_count >= min_fragments & fragment_stats$rank <= max_barcodes
  list(
    barcodes = fragment_stats$CB[keep],
    fragment_stats = fragment_stats,
    candidate_threshold = min_fragments,
    candidate_rank_cap = max_barcodes
  )
}

is_low_outlier <- function(dataframe, metric, nmads) {
  values <- dataframe[[metric]]
  center <- median(values, na.rm = TRUE)
  spread <- mad(values, na.rm = TRUE)
  if (is.na(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }
  values < center - nmads * spread
}

is_high_outlier <- function(dataframe, metric, nmads) {
  values <- dataframe[[metric]]
  center <- median(values, na.rm = TRUE)
  spread <- mad(values, na.rm = TRUE)
  if (is.na(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }
  values > center + nmads * spread
}

fail_if_na <- function(values) {
  values[is.na(values)] <- TRUE
  values
}

pass_if_na <- function(values) {
  values[is.na(values)] <- FALSE
  values
}

mix_color <- function(color, target = "#FFFFFF", amount = 0.5) {
  amount <- min(max(amount, 0), 1)
  source_rgb <- grDevices::col2rgb(color) / 255
  target_rgb <- grDevices::col2rgb(target) / 255
  blended <- source_rgb * (1 - amount) + target_rgb * amount
  grDevices::rgb(blended[1], blended[2], blended[3])
}

load_cima_hierarchy <- function(path) {
  hierarchy <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("cell_type_l1", "cell_type_l2", "cell_type_l3", "cell_type_l4")
  missing_cols <- setdiff(required_cols, colnames(hierarchy))
  if (length(missing_cols) > 0) {
    stop("CIMA hierarchy file is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  unique(hierarchy[, required_cols, drop = FALSE])
}

load_cima_reference_feature_model <- function(path) {
  model_df <- read.delim(gzfile(path), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("feature_index", "feature_id", "idf")
  missing_cols <- setdiff(required_cols, colnames(model_df))
  if (length(missing_cols) > 0) {
    stop("CIMA reference feature model is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  dim_cols <- grep("^dim_", colnames(model_df), value = TRUE)
  if (length(dim_cols) < 2) {
    stop("CIMA reference feature model must include at least two LSI dimensions")
  }
  list(
    feature_index = as.integer(model_df$feature_index),
    feature_id = as.character(model_df$feature_id),
    idf = as.numeric(model_df$idf),
    loadings = as.matrix(model_df[, dim_cols, drop = FALSE]),
    dims = dim_cols
  )
}

load_cima_reference_centroids <- function(path, id_col) {
  centroid_df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!id_col %in% colnames(centroid_df)) {
    stop("CIMA centroid table missing id column: ", id_col)
  }
  dim_cols <- grep("^dim_", colnames(centroid_df), value = TRUE)
  if (length(dim_cols) < 2) {
    stop("CIMA centroid table must include at least two LSI dimensions: ", path)
  }
  list(labels = centroid_df[[id_col]], embedding = as.matrix(centroid_df[, dim_cols, drop = FALSE]))
}

normalize_embedding_rows <- function(mat) {
  if (nrow(mat) == 0) {
    return(mat)
  }
  norms <- sqrt(rowSums(mat ^ 2))
  norms[norms == 0] <- 1
  sweep(mat, 1, norms, "/")
}

predict_cima_centroids <- function(query_embeddings, centroids) {
  dims_use <- seq.int(2, min(ncol(query_embeddings), ncol(centroids$embedding)))
  if (length(dims_use) < 2) {
    stop("CIMA centroid prediction requires at least two embedding dimensions")
  }
  query_norm <- normalize_embedding_rows(query_embeddings[, dims_use, drop = FALSE])
  centroid_norm <- normalize_embedding_rows(centroids$embedding[, dims_use, drop = FALSE])
  similarity <- query_norm %*% t(centroid_norm)
  top_index <- max.col(similarity, ties.method = "first")
  top_score <- similarity[cbind(seq_len(nrow(similarity)), top_index)]
  score_margin <- apply(similarity, 1, function(values) {
    ranked <- sort(values, decreasing = TRUE)
    if (length(ranked) < 2) {
      return(NA_real_)
    }
    ranked[1] - ranked[2]
  })
  list(label = centroids$labels[top_index], score = as.numeric(top_score), margin = as.numeric(score_margin))
}

annotate_cima_l1_reference_model <- function(counts, feature_model, l1_centroids) {
  selected_idx <- feature_model$feature_index
  if (max(selected_idx) > nrow(counts)) {
    stop("CIMA reference feature model expects more peaks than provided in query matrix")
  }
  observed_feature_ids <- rownames(counts)[selected_idx]
  if (!identical(observed_feature_ids, feature_model$feature_id)) {
    mismatch_idx <- which(observed_feature_ids != feature_model$feature_id)[1]
    stop(
      "CIMA reference feature identity mismatch at selected feature ",
      mismatch_idx,
      ": observed=", observed_feature_ids[[mismatch_idx]],
      " expected=", feature_model$feature_id[[mismatch_idx]]
    )
  }
  binary_counts <- counts[selected_idx, , drop = FALSE]
  if (length(binary_counts@x) > 0) {
    binary_counts@x[] <- 1
  }
  cell_peak_totals <- Matrix::colSums(binary_counts)
  cell_peak_totals[cell_peak_totals == 0] <- 1
  tf_counts <- binary_counts
  tf_counts@x <- tf_counts@x / rep.int(cell_peak_totals, diff(tf_counts@p))
  tfidf_counts <- Diagonal(x = feature_model$idf) %*% tf_counts
  query_embeddings <- as.matrix(Matrix::crossprod(tfidf_counts, feature_model$loadings))
  rownames(query_embeddings) <- colnames(counts)
  colnames(query_embeddings) <- feature_model$dims
  pred_l1 <- predict_cima_centroids(query_embeddings, l1_centroids)
  list(
    annotations = data.frame(
      cell_barcode = colnames(counts),
      cima_cell_type_l1 = pred_l1$label,
      cima_l1_score = as.numeric(pred_l1$score),
      cima_l1_score_margin = as.numeric(pred_l1$margin),
      stringsAsFactors = FALSE,
      row.names = colnames(counts)
    ),
    embeddings = query_embeddings
  )
}

build_l1_masked_labels <- function(labels, scores, margins, score_threshold = 0.4, margin_threshold = 0.1) {
  masked <- as.character(labels)
  masked[is.na(masked) | !nzchar(masked)] <- "Unknown"
  low_score <- !is.na(scores) & scores < score_threshold
  low_margin <- !is.na(margins) & margins < margin_threshold
  masked[low_score | low_margin] <- "Unknown"
  masked
}

build_cima_l1_palette <- function(hierarchy) {
  base_palette <- c(
    B = "#2C7BB6",
    CD4_T = "#D7191C",
    "CD8_T&unconvensional_T" = "#FDAE61",
    Myeloid = "#1A9641",
    ILC = "#762A83"
  )
  labels <- unique(hierarchy$cell_type_l1)
  fallback <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#FEE08B", "#F46D43", "#A50026")
  missing <- setdiff(labels, names(base_palette))
  if (length(missing) > 0) {
    extra <- fallback[seq_len(length(missing))]
    names(extra) <- missing
    base_palette <- c(base_palette, extra)
  }
  base_palette[labels]
}

cima_l1_display_labels <- function(labels) {
  display <- c(
    B = "B",
    CD4_T = "CD4 T",
    "CD8_T&unconvensional_T" = "CD8/unconv T",
    ILC = "ILC",
    Myeloid = "Myeloid",
    Unknown = "Unknown"
  )
  out <- as.character(labels)
  matched <- out %in% names(display)
  out[matched] <- unname(display[out[matched]])
  out
}

save_cima_l1_umap_plot <- function(plot_df, label_col, palette, out_path, title_text) {
  clean_df <- plot_df[!is.na(plot_df[[label_col]]) & nzchar(plot_df[[label_col]]), , drop = FALSE]
  if (nrow(clean_df) == 0) {
    plot_obj <- ggplot() + theme_void() + ggtitle(title_text) + annotate("text", x = 0.5, y = 0.5, label = paste0("No values available for ", label_col), size = 5) + xlim(0, 1) + ylim(0, 1)
  } else {
    clean_df[[label_col]] <- factor(clean_df[[label_col]], levels = names(palette))
    legend_labels <- stats::setNames(cima_l1_display_labels(names(palette)), names(palette))
    plot_obj <- ggplot(clean_df, aes(UMAP_1, UMAP_2, color = .data[[label_col]])) +
      geom_point(size = 0.25, alpha = 0.85) +
      scale_color_manual(values = palette, labels = legend_labels, drop = FALSE) +
      theme_bw(base_size = 11) +
      labs(title = title_text, x = "UMAP_1", y = "UMAP_2", color = "CIMA L1") +
      guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1), ncol = 1)) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.height = grid::unit(0.38, "cm"),
        legend.key.width = grid::unit(0.38, "cm"),
        legend.spacing.y = grid::unit(0.05, "cm")
      )
  }
  ggsave(out_path, plot_obj, width = 11.2, height = 8, dpi = 150, limitsize = FALSE)
}

select_lsi_dims <- function(atac_obj, lsi_embeddings, candidate_dims, max_qc_cor) {
  qc_metrics <- intersect(
    c("nCount_ATAC", "fragments", "TSS.enrichment", "FRiP", "nucleosome_signal"),
    colnames(atac_obj@meta.data)
  )
  if (length(qc_metrics) == 0) {
    return(list(dims = candidate_dims, audit = data.frame()))
  }

  audit_rows <- list()
  keep_dims <- c()
  for (dim_idx in candidate_dims) {
    dim_values <- lsi_embeddings[, dim_idx]
    cors <- vapply(qc_metrics, function(metric) {
      suppressWarnings(abs(cor(dim_values, atac_obj@meta.data[[metric]], method = "spearman", use = "pairwise.complete.obs")))
    }, numeric(1))
    max_cor <- max(cors, na.rm = TRUE)
    top_metric <- names(cors)[which.max(cors)]
    keep <- is.na(max_cor) || max_cor <= max_qc_cor
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      dim = dim_idx,
      max_abs_spearman_qc_cor = max_cor,
      top_qc_metric = top_metric,
      keep = keep,
      stringsAsFactors = FALSE
    )
    if (keep) {
      keep_dims <- c(keep_dims, dim_idx)
    }
  }
  audit <- do.call(rbind, audit_rows)
  if (length(keep_dims) < 2) {
    keep_dims <- candidate_dims
    audit$keep <- TRUE
    audit$forced_keep_all <- TRUE
  } else {
    audit$forced_keep_all <- FALSE
  }
  list(dims = keep_dims, audit = audit)
}

run_single_sample_umap <- function(atac_obj, min_cutoff, resolution, max_qc_cor) {
  DefaultAssay(atac_obj) <- "ATAC"
  atac_obj <- RunTFIDF(atac_obj)
  atac_obj <- FindTopFeatures(atac_obj, min.cutoff = min_cutoff)
  atac_obj <- RunSVD(atac_obj)
  lsi_embeddings <- Embeddings(atac_obj[["lsi"]])
  candidate_dims <- seq.int(2, min(30, ncol(lsi_embeddings)))
  if (length(candidate_dims) < 2) {
    stop("Not enough LSI dimensions available")
  }
  dim_selection <- select_lsi_dims(atac_obj, lsi_embeddings, candidate_dims, max_qc_cor)
  dims_use <- dim_selection$dims
  atac_obj <- FindNeighbors(atac_obj, reduction = "lsi", dims = dims_use, verbose = FALSE)
  atac_obj <- FindClusters(atac_obj, resolution = resolution, verbose = FALSE)
  atac_obj <- RunUMAP(
    atac_obj,
    reduction = "lsi",
    dims = dims_use,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )
  list(object = atac_obj, dims = paste(dims_use, collapse = ","), lsi_qc_cor = dim_selection$audit)
}

option_list <- list(
  make_option(c("--gse"), type = "character"),
  make_option(c("--gsm"), type = "character"),
  make_option(c("--individual-id"), type = "character", default = ""),
  make_option(c("--nmads"), type = "numeric", default = 4),
  make_option(c("--output-root"), type = "character", default = file.path(project_root, "output", "co", "atac")),
  make_option(c("--min-inferred-fragments"), type = "numeric", default = 1000),
  make_option(c("--max-inferred-barcodes"), type = "integer", default = 20000),
  make_option(c("--min-counts"), type = "numeric", default = 50),
  make_option(c("--min-fragments"), type = "numeric", default = 1000),
  make_option(c("--min-tss"), type = "numeric", default = 0.5),
  make_option(c("--min-frip"), type = "numeric", default = 0.02),
  make_option(c("--max-nucleosome"), type = "numeric", default = 4.0),
  make_option(c("--max-blacklist"), type = "numeric", default = 0.05),
  make_option(c("--feature-cutoff"), type = "character", default = "q5"),
  make_option(c("--cluster-resolution"), type = "numeric", default = 0.5),
  make_option(c("--max-lsi-qc-cor"), type = "numeric", default = 0.35),
  make_option(c("--fragment-file"), type = "character", default = ""),
  make_option(c("--filtered-barcodes"), type = "character", default = ""),
  make_option(c("--sample-label"), type = "character", default = ""),
  make_option(c("--annotate-cima-l1"), action = "store_true", default = TRUE),
  make_option(c("--cima-l1-score-threshold"), type = "numeric", default = 0.4),
  make_option(c("--cima-l1-margin-threshold"), type = "numeric", default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$gse) || is.null(opt$gsm)) {
  stop("Both --gse and --gsm are required")
}

raw_dir <- file.path(project_root, "data", "raw", opt$gse)
reference_dir <- file.path(project_root, "data", "reference")
peak_file <- file.path(reference_dir, "peak.bed")
cima_dir <- file.path(reference_dir, "cima")
cima_hierarchy_file <- file.path(cima_dir, "cima_atac_celltype_hierarchy.csv")
cima_feature_model_file <- file.path(cima_dir, "cima_atac_reference_lsi_features.tsv.gz")
cima_l1_centroid_file <- file.path(cima_dir, "cima_atac_reference_l1_centroids.tsv")
output_dir <- file.path(opt$`output-root`, opt$gse, opt$gsm)
matrix_dir <- file.path(output_dir, "matrix")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

fragment_file <- if (nzchar(opt$`fragment-file`)) {
  normalizePath(opt$`fragment-file`, winslash = "/", mustWork = TRUE)
} else {
  find_single_file(raw_dir, opt$gsm, "fragments", required = TRUE)
}
barcode_file <- if (nzchar(opt$`filtered-barcodes`)) {
  normalizePath(opt$`filtered-barcodes`, winslash = "/", mustWork = TRUE)
} else if (dir.exists(raw_dir)) {
  find_single_file(raw_dir, opt$gsm, "filtered_barcodes", required = FALSE)
} else {
  NULL
}
sample_label <- if (nzchar(opt$`sample-label`)) opt$`sample-label` else opt$gsm

cat(rep("=", 80), "\n", sep = "")
cat("CO scATAC-seq QC\n")
cat(rep("=", 80), "\n", sep = "")
cat("GSE ID:", opt$gse, "\n")
cat("GSM ID:", opt$gsm, "\n")
cat("Fragment file:", fragment_file, "\n")
cat("Output directory:", output_dir, "\n")

ensure_tabix_index(fragment_file)
barcode_source <- "filtered_barcodes"
barcode_candidate_threshold <- NA_real_
barcode_candidate_rank_cap <- NA_real_
barcode_prefilter_stats <- NULL

if (!is.null(barcode_file)) {
  cat("Barcode file:", barcode_file, "\n\n")
  barcodes <- readLines(gzfile(barcode_file))
  barcodes <- unique(barcodes[nzchar(barcodes)])
} else {
  cat("Barcode file: not found, using strict fragment-count barcode inference\n\n")
  inferred <- infer_barcodes(fragment_file, opt$`min-inferred-fragments`, opt$`max-inferred-barcodes`)
  barcodes <- inferred$barcodes
  barcode_source <- "strict_fragments"
  barcode_candidate_threshold <- inferred$candidate_threshold
  barcode_candidate_rank_cap <- inferred$candidate_rank_cap
  barcode_prefilter_stats <- inferred$fragment_stats
}

if (length(barcodes) == 0) {
  stop("No barcodes selected for sample ", opt$gsm)
}
cat("Initial barcodes:", length(barcodes), "(", barcode_source, ")\n\n")

peaks_gr <- rtracklayer::import(peak_file)
seqlevelsStyle(peaks_gr) <- "UCSC"
genome(peaks_gr) <- "hg38"
peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
suppressWarnings(seqlevelsStyle(annotations) <- "UCSC")
genome(annotations) <- "hg38"
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
data("blacklist_hg38_unified", package = "Signac")
suppressWarnings(seqlevelsStyle(blacklist_hg38_unified) <- "UCSC")
genome(blacklist_hg38_unified) <- "hg38"
blacklist_hg38_unified <- keepStandardChromosomes(blacklist_hg38_unified, pruning.mode = "coarse")
common_seq <- intersect(seqlevels(peaks_gr), seqlevels(annotations))
peaks_gr <- keepSeqlevels(peaks_gr, common_seq, pruning.mode = "coarse")
annotations <- keepSeqlevels(annotations, common_seq, pruning.mode = "coarse")
blacklist_hg38_unified <- keepSeqlevels(blacklist_hg38_unified, common_seq, pruning.mode = "coarse")
seqinfo(peaks_gr) <- seqinfo(annotations)[common_seq]
project_peak_ids <- paste0(as.character(seqnames(peaks_gr)), ":", start(peaks_gr) - 1L, "-", end(peaks_gr))

frags <- CreateFragmentObject(path = fragment_file, cells = barcodes)
counts <- FeatureMatrix(fragments = frags, features = peaks_gr, cells = barcodes)
rownames(counts) <- project_peak_ids
atac_assay <- CreateChromatinAssay(counts = counts, fragments = frags, genome = "hg38", sep = c(":", "-"))
slot(atac_assay, "annotation") <- annotations
atac_obj <- CreateSeuratObject(counts = atac_assay, assay = "ATAC", project = sample_label)
atac_obj$sample <- sample_label
atac_obj$dataset <- opt$gse
atac_obj$individual_id <- opt$`individual-id`
atac_obj$barcode_source <- barcode_source

cat("Cells:", ncol(atac_obj), "\n")
cat("Features:", nrow(atac_obj), "\n\n")

cat("Computing QC metrics...\n")
atac_obj <- NucleosomeSignal(atac_obj)
atac_obj <- TSSEnrichment(atac_obj, fast = FALSE)
if (!is.null(barcode_prefilter_stats)) {
  fragment_stats <- barcode_prefilter_stats[barcode_prefilter_stats$CB %in% colnames(atac_obj), , drop = FALSE]
} else {
  fragment_stats <- CountFragments(fragment_file, cells = colnames(atac_obj))
}
rownames(fragment_stats) <- fragment_stats$CB
atac_obj$fragments <- fragment_stats[colnames(atac_obj), "frequency_count"]
atac_obj$total_fragments <- atac_obj$fragments
atac_obj <- FRiP(object = atac_obj, assay = "ATAC", total.fragments = "fragments")
atac_obj$blacklist_fraction <- FractionCountsInRegion(object = atac_obj, assay = "ATAC", regions = blacklist_hg38_unified)

cat("Running doublet detection...\n")
counts_matrix <- GetAssayData(atac_obj[["ATAC"]], layer = "counts")
sce <- SingleCellExperiment(list(counts = counts_matrix))
sce <- scDblFinder(sce, clusters = TRUE, aggregateFeatures = TRUE, nfeatures = 25, processing = "normFeatures")
doublet_info <- as.data.frame(colData(sce)[, c("scDblFinder.class", "scDblFinder.score")])
atac_obj$scDblFinder.class <- doublet_info$scDblFinder.class
atac_obj$scDblFinder.score <- doublet_info$scDblFinder.score

cat("Applying co-ATAC QC filters...\n")
meta_data <- atac_obj@meta.data
meta_data$fails_count_floor <- fail_if_na(meta_data$nCount_ATAC < opt$`min-counts`)
meta_data$fails_fragment_floor <- fail_if_na(meta_data$fragments < opt$`min-fragments`)
meta_data$fails_tss_floor <- fail_if_na(meta_data$TSS.enrichment < opt$`min-tss`)
meta_data$fails_frip_floor <- fail_if_na(meta_data$FRiP < opt$`min-frip`)
meta_data$fails_nucleosome_ceiling <- pass_if_na(meta_data$nucleosome_signal > opt$`max-nucleosome`)
meta_data$fails_blacklist_ceiling <- pass_if_na(meta_data$blacklist_fraction > opt$`max-blacklist`)
meta_data$fails_count_mad_low <- pass_if_na(is_low_outlier(meta_data, "nCount_ATAC", opt$nmads))
meta_data$fails_tss_mad_low <- pass_if_na(is_low_outlier(meta_data, "TSS.enrichment", opt$nmads))
meta_data$fails_frip_mad_low <- pass_if_na(is_low_outlier(meta_data, "FRiP", opt$nmads))
meta_data$fails_count_mad_high <- pass_if_na(is_high_outlier(meta_data, "nCount_ATAC", opt$nmads))
meta_data$fails_doublet <- fail_if_na(meta_data$scDblFinder.class == "doublet")
meta_data$pass_qc <- !(
  meta_data$fails_count_floor |
    meta_data$fails_fragment_floor |
    meta_data$fails_tss_floor |
    meta_data$fails_frip_floor |
    meta_data$fails_nucleosome_ceiling |
    meta_data$fails_blacklist_ceiling |
    meta_data$fails_count_mad_low |
    meta_data$fails_tss_mad_low |
    meta_data$fails_frip_mad_low |
    meta_data$fails_count_mad_high |
    meta_data$fails_doublet
)
atac_obj@meta.data <- meta_data
fail_cols <- grep("^fails_", colnames(meta_data), value = TRUE)
for (fail_col in fail_cols) {
  cat(fail_col, ":", sum(meta_data[[fail_col]] %in% TRUE, na.rm = TRUE), "\n")
}
qc_cells <- rownames(meta_data)[meta_data$pass_qc %in% TRUE]
cat("Cells before QC:", nrow(meta_data), "\n")
cat("Cells after QC:", length(qc_cells), "\n\n")
if (length(qc_cells) < 3) {
  stop("Too few QC-passing cells for co-ATAC UMAP: ", length(qc_cells))
}
qc_counts <- GetAssayData(atac_obj[["ATAC"]], layer = "counts")[, qc_cells, drop = FALSE]
rownames(qc_counts) <- project_peak_ids

cima_l1_status <- "disabled"
cima_l1_detail <- ""
if (isTRUE(opt$`annotate-cima-l1`)) {
  cat("Running CIMA ATAC L1 annotation...\n")
  tryCatch({
    cima_hierarchy <- load_cima_hierarchy(cima_hierarchy_file)
    cima_feature_model <- load_cima_reference_feature_model(cima_feature_model_file)
    cima_l1_centroids <- load_cima_reference_centroids(cima_l1_centroid_file, "cell_type_l1")
    cima_result <- annotate_cima_l1_reference_model(qc_counts, cima_feature_model, cima_l1_centroids)
    annotation_df <- cima_result$annotations
    annotation_df$cima_cell_type_l1_masked <- build_l1_masked_labels(
      annotation_df$cima_cell_type_l1,
      annotation_df$cima_l1_score,
      annotation_df$cima_l1_score_margin,
      score_threshold = opt$`cima-l1-score-threshold`,
      margin_threshold = opt$`cima-l1-margin-threshold`
    )
    annotation_df$cima_l1_low_confidence <- annotation_df$cima_cell_type_l1_masked == "Unknown"
    annotation_cols <- c(
      "cima_cell_type_l1",
      "cima_cell_type_l1_masked",
      "cima_l1_low_confidence",
      "cima_l1_score",
      "cima_l1_score_margin"
    )
    atac_obj@meta.data[annotation_df$cell_barcode, annotation_cols] <- annotation_df[, annotation_cols]
    cima_l1_status <- "ok"
    cima_l1_detail <- "reference_projected_lsi_l1"
  }, error = function(e) {
    cima_l1_status <<- "failed"
    cima_l1_detail <<- conditionMessage(e)
    warning("CIMA ATAC L1 annotation failed: ", conditionMessage(e))
  })
}

cat("Running query-native LSI/UMAP...\n")
atac_qc_obj <- subset(atac_obj, cells = qc_cells)
umap_result <- run_single_sample_umap(atac_qc_obj, opt$`feature-cutoff`, opt$`cluster-resolution`, opt$`max-lsi-qc-cor`)
atac_qc_obj <- umap_result$object
atac_obj@meta.data[rownames(atac_qc_obj@meta.data), "seurat_clusters"] <- as.character(atac_qc_obj$seurat_clusters)
umap_embeddings <- Embeddings(atac_qc_obj[["umap"]])
atac_obj@meta.data[rownames(umap_embeddings), c("umap_atac_1", "umap_atac_2")] <- umap_embeddings[, 1:2, drop = FALSE]
if (cima_l1_status == "ok") {
  l1_cols <- c("cima_cell_type_l1", "cima_cell_type_l1_masked", "cima_l1_low_confidence", "cima_l1_score", "cima_l1_score_margin")
  atac_qc_obj@meta.data[, l1_cols] <- atac_obj@meta.data[rownames(atac_qc_obj@meta.data), l1_cols, drop = FALSE]
}

cluster_df <- data.frame(
  UMAP_1 = umap_embeddings[, 1],
  UMAP_2 = umap_embeddings[, 2],
  cluster = as.character(atac_qc_obj$seurat_clusters),
  nCount_ATAC = atac_qc_obj$nCount_ATAC,
  TSS.enrichment = atac_qc_obj$TSS.enrichment,
  FRiP = atac_qc_obj$FRiP,
  fragments = atac_qc_obj$fragments,
  stringsAsFactors = FALSE
)

umap_cluster <- ggplot(cluster_df, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.25, alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(title = paste0(sample_label, " co-ATAC query-native UMAP by cluster"), color = "Cluster")
ggsave(file.path(output_dir, "umap_atac_clusters.png"), umap_cluster, width = 10, height = 8, dpi = 150)

qc_feature_umap <- (FeaturePlot(atac_qc_obj, features = "nCount_ATAC", reduction = "umap") |
  FeaturePlot(atac_qc_obj, features = "TSS.enrichment", reduction = "umap")) /
  (FeaturePlot(atac_qc_obj, features = "FRiP", reduction = "umap") |
    FeaturePlot(atac_qc_obj, features = "fragments", reduction = "umap"))
ggsave(file.path(output_dir, "umap_atac_qc_features.png"), qc_feature_umap, width = 14, height = 11, dpi = 150)

if (cima_l1_status == "ok") {
  cima_hierarchy <- load_cima_hierarchy(cima_hierarchy_file)
  l1_palette <- build_cima_l1_palette(cima_hierarchy)
  l1_masked_palette <- c(l1_palette, Unknown = "#7F7F7F")
  l1_df <- data.frame(
    UMAP_1 = umap_embeddings[, 1],
    UMAP_2 = umap_embeddings[, 2],
    cell_barcode = rownames(umap_embeddings),
    atac_qc_obj@meta.data[, c("cima_cell_type_l1", "cima_cell_type_l1_masked"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  save_cima_l1_umap_plot(l1_df, "cima_cell_type_l1", l1_palette, file.path(output_dir, "umap_atac_cima_l1.png"), paste0(sample_label, " query-native UMAP by CIMA ATAC L1"))
  save_cima_l1_umap_plot(l1_df, "cima_cell_type_l1_masked", l1_masked_palette, file.path(output_dir, "umap_atac_cima_l1_masked.png"), paste0(sample_label, " query-native UMAP by CIMA ATAC L1 (masked)"))
}

p1 <- safe_plot("QC Metrics", VlnPlot(atac_obj, features = c("nCount_ATAC", "fragments", "TSS.enrichment", "FRiP", "blacklist_fraction", "nucleosome_signal"), group.by = "pass_qc", pt.size = 0, ncol = 3))
p2 <- safe_plot("Count vs TSS", FeatureScatter(atac_obj, feature1 = "nCount_ATAC", feature2 = "TSS.enrichment", group.by = "pass_qc"))
p3 <- safe_plot("Count vs FRiP", FeatureScatter(atac_obj, feature1 = "nCount_ATAC", feature2 = "FRiP", group.by = "pass_qc"))
p4 <- safe_plot("Fragment Histogram", FragmentHistogram(atac_obj, group.by = "pass_qc"))
qc_overview <- wrap_plots(list(p1, p2, p3, p4), ncol = 2) +
  plot_annotation(title = paste0(sample_label, " co-ATAC QC Overview"))
ggsave(file.path(output_dir, "qc_overview.png"), qc_overview, width = 20, height = 16, dpi = 150, limitsize = FALSE)

if (!is.null(barcode_prefilter_stats)) {
  barcode_rank <- ggplot(barcode_prefilter_stats, aes(rank, frequency_count)) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = barcode_candidate_threshold, linetype = "dashed", color = "#d55e00") +
    geom_vline(xintercept = barcode_candidate_rank_cap, linetype = "dotdash", color = "#0072b2") +
    scale_x_log10() + scale_y_log10() + theme_bw(base_size = 11) +
    labs(title = "Strict barcode rank", x = "Barcode rank", y = "Fragment count")
  ggsave(file.path(output_dir, "barcode_rank.png"), barcode_rank, width = 8, height = 6, dpi = 150)
}

if (nrow(umap_result$lsi_qc_cor) > 0) {
  write.csv(umap_result$lsi_qc_cor, file.path(output_dir, "lsi_qc_correlation.csv"), row.names = FALSE)
}

meta_data <- atac_obj@meta.data
metadata_all <- cbind(cell_barcode = rownames(meta_data), meta_data)
metadata_qc <- metadata_all[metadata_all$pass_qc %in% TRUE, , drop = FALSE]

summary_stats <- data.frame(
  metric = c(
    "gse", "gsm", "individual_id", "fragment_file", "barcode_source",
    "barcode_candidate_threshold", "barcode_candidate_rank_cap",
    "input_cells", "pass_qc", "qc_rate", "singlets", "doublets",
    "mean_nCount_ATAC", "median_TSS_enrichment", "median_FRiP", "median_fragments",
    "median_blacklist_fraction", "median_nucleosome_signal", "query_cluster_count",
    "cima_l1_status", "cima_l1_detail", "cima_l1_low_confidence_count", "cima_l1_low_confidence_frac",
    "median_cima_l1_score", "median_cima_l1_score_margin",
    "min_counts", "min_fragments", "min_tss", "min_frip", "max_nucleosome", "max_blacklist",
    "feature_cutoff", "max_lsi_qc_cor", "lsi_dims"
  ),
  value = c(
    opt$gse, sample_label, opt$`individual-id`, fragment_file, barcode_source,
    if (is.na(barcode_candidate_threshold)) "" else sprintf("%.0f", barcode_candidate_threshold),
    if (is.na(barcode_candidate_rank_cap)) "" else sprintf("%.0f", barcode_candidate_rank_cap),
    ncol(atac_obj), nrow(metadata_qc), sprintf("%.2f%%", nrow(metadata_qc) / ncol(atac_obj) * 100),
    sum(atac_obj$scDblFinder.class == "singlet", na.rm = TRUE),
    sum(atac_obj$scDblFinder.class == "doublet", na.rm = TRUE),
    sprintf("%.0f", mean(metadata_qc$nCount_ATAC, na.rm = TRUE)),
    sprintf("%.2f", median(metadata_qc$TSS.enrichment, na.rm = TRUE)),
    sprintf("%.4f", median(metadata_qc$FRiP, na.rm = TRUE)),
    sprintf("%.0f", median(metadata_qc$fragments, na.rm = TRUE)),
    sprintf("%.4f", median(metadata_qc$blacklist_fraction, na.rm = TRUE)),
    sprintf("%.2f", median(metadata_qc$nucleosome_signal, na.rm = TRUE)),
    length(unique(metadata_qc$seurat_clusters[!is.na(metadata_qc$seurat_clusters)])),
    cima_l1_status, cima_l1_detail,
    if ("cima_l1_low_confidence" %in% colnames(metadata_qc)) sum(metadata_qc$cima_l1_low_confidence %in% TRUE, na.rm = TRUE) else "",
    if ("cima_l1_low_confidence" %in% colnames(metadata_qc)) sprintf("%.2f%%", mean(metadata_qc$cima_l1_low_confidence %in% TRUE, na.rm = TRUE) * 100) else "",
    if ("cima_l1_score" %in% colnames(metadata_qc)) sprintf("%.4f", median(metadata_qc$cima_l1_score, na.rm = TRUE)) else "",
    if ("cima_l1_score_margin" %in% colnames(metadata_qc)) sprintf("%.4f", median(metadata_qc$cima_l1_score_margin, na.rm = TRUE)) else "",
    opt$`min-counts`, opt$`min-fragments`, opt$`min-tss`, opt$`min-frip`, opt$`max-nucleosome`, opt$`max-blacklist`,
    opt$`feature-cutoff`, opt$`max-lsi-qc-cor`, umap_result$dims
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_stats, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)
write.csv(metadata_all, file.path(output_dir, "metadata.csv"), row.names = FALSE)
write.csv(metadata_qc, file.path(output_dir, "metadata_qc.csv"), row.names = FALSE)
write.csv(metadata_qc[, intersect(c("cell_barcode", "seurat_clusters", "nCount_ATAC", "fragments", "TSS.enrichment", "FRiP", "blacklist_fraction", "nucleosome_signal", "scDblFinder.class", "cima_cell_type_l1", "cima_cell_type_l1_masked", "cima_l1_low_confidence", "cima_l1_score", "cima_l1_score_margin", "umap_atac_1", "umap_atac_2"), colnames(metadata_qc)), drop = FALSE], file.path(output_dir, "validation_result.csv"), row.names = FALSE)
writeMM(qc_counts, file.path(matrix_dir, "matrix.mtx"))
write_lines_gz(rownames(qc_counts), file.path(matrix_dir, "features.tsv.gz"))
write_lines_gz(colnames(qc_counts), file.path(matrix_dir, "barcodes.tsv.gz"))
saveRDS(atac_obj, file.path(output_dir, paste0(opt$gsm, "_co_atac_qc.rds")))

cat("Summary:", file.path(output_dir, "qc_summary.csv"), "\n")
cat("UMAP clusters:", file.path(output_dir, "umap_atac_clusters.png"), "\n")
cat("Metadata:", file.path(output_dir, "metadata.csv"), "\n")
