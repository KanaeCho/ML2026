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

  if (length(files) == 1) {
    return(files[1])
  }

  filtered <- files[grepl("filtered_barcodes", basename(files), ignore.case = TRUE)]
  if (token == "barcodes" && length(filtered) == 1) {
    return(filtered[1])
  }

  stop(
    "Multiple candidate files found for ", gsm, " with token '", token, "': ",
    paste(basename(files), collapse = ", ")
  )
}

normalize_metric_colnames <- function(values) {
  normalized <- tolower(gsub("[^A-Za-z0-9]+", "_", values))
  normalized <- gsub("_+", "_", normalized)
  gsub("^_|_$", "", normalized)
}

coerce_logical_flag <- function(values) {
  lower <- tolower(trimws(as.character(values)))
  lower[lower %in% c("", "na", "nan", "none", "null")] <- NA_character_
  ifelse(lower %in% c("1", "true", "t", "yes", "y"), TRUE,
         ifelse(lower %in% c("0", "false", "f", "no", "n"), FALSE, NA))
}

load_singlecell_barcodes <- function(singlecell_file) {
  metrics <- read.csv(gzfile(singlecell_file), stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(metrics) == 0) {
    stop("singlecell file is empty: ", singlecell_file)
  }

  original_names <- colnames(metrics)
  normalized_names <- normalize_metric_colnames(original_names)
  colnames(metrics) <- normalized_names

  barcode_col <- intersect(c("barcode", "cb"), normalized_names)
  if (length(barcode_col) == 0) {
    stop("Unable to find barcode column in singlecell file: ", singlecell_file)
  }
  barcode_col <- barcode_col[1]
  metrics$barcode <- as.character(metrics[[barcode_col]])

  selector <- NULL
  selector_source <- NULL

  flag_cols <- intersect(
    c("is_cell_barcode", "is_cell", "cell_barcode", "cell_called"),
    normalized_names
  )
  if (length(flag_cols) > 0) {
    selector <- coerce_logical_flag(metrics[[flag_cols[1]]])
    selector_source <- flag_cols[1]
  }

  if (is.null(selector) || all(is.na(selector))) {
    cell_id_cols <- intersect(c("cell_id", "cellid"), normalized_names)
    if (length(cell_id_cols) > 0) {
      cell_ids <- trimws(as.character(metrics[[cell_id_cols[1]]]))
      selector <- !(is.na(cell_ids) | cell_ids == "" | tolower(cell_ids) %in% c("none", "na", "nan"))
      selector_source <- cell_id_cols[1]
    }
  }

  if (is.null(selector) || all(is.na(selector))) {
    count_cols <- intersect(c("passed_filters", "peak_region_fragments", "on_target_fragments"), normalized_names)
    if (length(count_cols) > 0) {
      selector <- suppressWarnings(as.numeric(metrics[[count_cols[1]]]) > 0)
      selector_source <- count_cols[1]
    }
  }

  if (is.null(selector) || all(is.na(selector))) {
    stop("Unable to infer cell barcode flag from singlecell file: ", singlecell_file)
  }

  selector[is.na(selector)] <- FALSE
  barcodes <- unique(metrics$barcode[selector & nzchar(metrics$barcode)])
  metrics_selected <- metrics[metrics$barcode %in% barcodes, , drop = FALSE]

  list(
    barcodes = barcodes,
    metrics = metrics_selected,
    selector_source = selector_source
  )
}

is_outlier <- function(dataframe, metric, nmads) {
  if (!metric %in% colnames(dataframe)) {
    stop("Metric not found in dataframe: ", metric)
  }

  values <- dataframe[[metric]]
  center <- median(values, na.rm = TRUE)
  spread <- mad(values, na.rm = TRUE)

  if (is.na(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }

  (values < center - nmads * spread) | (center + nmads * spread < values)
}

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

make_line_df <- function(name, value) {
  if (is.na(value) || !is.finite(value) || value <= 0) {
    return(data.frame(name = character(), value = numeric()))
  }
  data.frame(name = name, value = value, stringsAsFactors = FALSE)
}

combine_line_df <- function(...) {
  items <- Filter(function(df) nrow(df) > 0, list(...))
  if (length(items) == 0) {
    return(data.frame(name = character(), value = numeric()))
  }
  do.call(rbind, items)
}

infer_barcode_threshold <- function(counts) {
  counts <- sort(counts[counts > 0], decreasing = TRUE)

  if (length(counts) == 0) {
    stop("No positive fragment counts available for barcode inference")
  }

  if (length(counts) < 100) {
    return(tail(counts, n = 1))
  }

  ranks <- seq_along(counts)
  x <- log10(ranks)
  y <- log10(counts)
  fit <- smooth.spline(x = x, y = y, spar = 0.6)
  d1 <- predict(fit, x = x, deriv = 1)$y
  d2 <- predict(fit, x = x, deriv = 2)$y
  curvature <- abs(d2) / (1 + d1^2)^(3 / 2)
  knee_idx <- which.max(curvature)
  counts[knee_idx]
}

infer_barcode_inflection <- function(counts, search_end = NULL) {
  counts <- sort(counts[counts > 0], decreasing = TRUE)

  if (length(counts) == 0) {
    stop("No positive fragment counts available for barcode inflection inference")
  }

  if (length(counts) < 100) {
    return(list(threshold = tail(counts, n = 1), rank = length(counts)))
  }

  ranks <- seq_along(counts)
  x <- log10(ranks)
  y <- log10(counts)
  fit <- smooth.spline(x = x, y = y, spar = 0.6)
  d2 <- predict(fit, x = x, deriv = 2)$y

  # Focus on the early portion of the rank curve where true-cell signal should live.
  if (is.null(search_end) || !is.finite(search_end)) {
    search_end <- min(length(counts), max(1000, floor(length(counts) * 0.1)))
  }
  search_end <- min(length(counts), max(1000, as.integer(search_end)))
  search_idx <- seq_len(search_end)
  inflection_idx <- search_idx[which.min(d2[search_idx])]

  list(threshold = counts[inflection_idx], rank = inflection_idx)
}

infer_barcode_floor_auto <- function(counts) {
  positive_counts <- counts[counts > 0]
  if (length(positive_counts) < 100) {
    return(NA_real_)
  }

  log_counts <- log10(positive_counts + 1)
  dens <- density(log_counts, n = 2048, adjust = 1)
  maxima <- which(diff(sign(diff(dens$y))) == -2) + 1
  minima <- which(diff(sign(diff(dens$y))) == 2) + 1

  if (length(maxima) < 2 || length(minima) == 0) {
    return(NA_real_)
  }

  best_score <- -Inf
  best_valley <- NA_integer_

  for (i in seq_len(length(maxima) - 1)) {
    for (j in seq.int(i + 1, length(maxima))) {
      left_peak <- maxima[i]
      right_peak <- maxima[j]
      if (left_peak > right_peak) {
        tmp <- left_peak
        left_peak <- right_peak
        right_peak <- tmp
      }

      valley_candidates <- minima[minima > left_peak & minima < right_peak]
      if (length(valley_candidates) == 0) {
        next
      }

      valley_idx <- valley_candidates[which.min(dens$y[valley_candidates])]
      separation <- dens$x[right_peak] - dens$x[left_peak]
      if (separation < 0.15) {
        next
      }

      valley_depth <- min(dens$y[left_peak], dens$y[right_peak]) - dens$y[valley_idx]
      if (valley_depth <= 0) {
        next
      }

      score <- separation * valley_depth
      if (score > best_score) {
        best_score <- score
        best_valley <- valley_idx
      }
    }
  }

  if (is.na(best_valley)) {
    return(NA_real_)
  }

  max(1, ceiling((10 ^ dens$x[best_valley]) - 1))
}

infer_barcodes <- function(fragment_file) {
  fragment_stats <- CountFragments(fragment_file)
  fragment_stats <- fragment_stats[!is.na(fragment_stats$CB) & nzchar(fragment_stats$CB), , drop = FALSE]
  fragment_stats <- fragment_stats[fragment_stats$frequency_count > 0, , drop = FALSE]
  sorted_counts <- sort(fragment_stats$frequency_count, decreasing = TRUE)

  knee_threshold <- infer_barcode_threshold(sorted_counts)
  floor_auto <- infer_barcode_floor_auto(sorted_counts)
  if (is.na(floor_auto) || !is.finite(floor_auto)) {
    floor_auto <- max(10, ceiling(knee_threshold * 0.1))
  }
  floor_rank <- max(which(sorted_counts >= floor_auto))
  inflection_search_end <- min(floor_rank, max(2000, ceiling(floor_rank * 0.25)))
  inflection <- infer_barcode_inflection(sorted_counts, search_end = inflection_search_end)
  candidate_threshold <- max(knee_threshold, floor_auto, inflection$threshold)

  rank_df <- data.frame(
    CB = fragment_stats$CB,
    frequency_count = fragment_stats$frequency_count,
    stringsAsFactors = FALSE
  )
  rank_df <- rank_df[order(rank_df$frequency_count, decreasing = TRUE), , drop = FALSE]
  rank_df$rank <- seq_len(nrow(rank_df))
  candidate_rank_cap <- max(1, inflection$rank)

  keep <- fragment_stats$frequency_count >= candidate_threshold &
    rank_df$rank[match(fragment_stats$CB, rank_df$CB)] <= candidate_rank_cap

  list(
    barcodes = fragment_stats$CB[keep],
    knee_threshold = knee_threshold,
    inflection_threshold = inflection$threshold,
    candidate_rank_cap = candidate_rank_cap,
    inflection_search_end = inflection_search_end,
    floor_auto = floor_auto,
    candidate_threshold = candidate_threshold,
    total_candidates = nrow(fragment_stats),
    fragment_stats = rank_df
  )
}

write_lines_gz <- function(lines, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n")
}

make_peak_ids <- function(gr) {
  paste0(as.character(seqnames(gr)), ":", start(gr) - 1L, "-", end(gr))
}

mix_color <- function(color, target = "#FFFFFF", amount = 0.5) {
  amount <- min(max(amount, 0), 1)
  source_rgb <- grDevices::col2rgb(color) / 255
  target_rgb <- grDevices::col2rgb(target) / 255
  blended <- source_rgb * (1 - amount) + target_rgb * amount
  grDevices::rgb(blended[1], blended[2], blended[3])
}

make_shade_palette <- function(base_color, labels) {
  labels <- unique(labels)
  n <- length(labels)
  if (n == 0) {
    return(setNames(character(), character()))
  }
  if (n == 1) {
    return(stats::setNames(base_color, labels))
  }
  shades <- grDevices::colorRampPalette(
    c(
      mix_color(base_color, "#FFFFFF", 0.55),
      base_color,
      mix_color(base_color, "#000000", 0.18)
    )
  )(n)
  stats::setNames(shades, labels)
}

load_cima_hierarchy <- function(path) {
  hierarchy <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("cell_type_l1", "cell_type_l2", "cell_type_l3", "cell_type_l4")
  missing_cols <- setdiff(required_cols, colnames(hierarchy))
  if (length(missing_cols) > 0) {
    stop("CIMA hierarchy file is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  hierarchy <- unique(hierarchy[, required_cols, drop = FALSE])
  if (anyDuplicated(hierarchy$cell_type_l4)) {
    stop("CIMA hierarchy has non-unique cell_type_l4 labels")
  }
  hierarchy[order(hierarchy$cell_type_l1, hierarchy$cell_type_l2, hierarchy$cell_type_l3, hierarchy$cell_type_l4), , drop = FALSE]
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

  list(
    labels = centroid_df[[id_col]],
    embedding = as.matrix(centroid_df[, dim_cols, drop = FALSE])
  )
}

build_cima_palettes <- function(hierarchy) {
  base_palette_defaults <- c(
    B = "#2C7BB6",
    CD4_T = "#D7191C",
    "CD8_T&unconvensional_T" = "#FDAE61",
    Myeloid = "#1A9641",
    ILC = "#762A83"
  )
  fallback_colors <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#FEE08B", "#F46D43", "#A50026")

  l1_labels <- unique(hierarchy$cell_type_l1)
  missing_l1 <- setdiff(l1_labels, names(base_palette_defaults))
  if (length(missing_l1) > 0) {
    extra_colors <- fallback_colors[seq_len(length(missing_l1))]
    names(extra_colors) <- missing_l1
    base_palette_defaults <- c(base_palette_defaults, extra_colors)
  }
  l1_palette <- base_palette_defaults[l1_labels]

  build_level_palette <- function(level_col) {
    palette <- c()
    for (l1 in l1_labels) {
      labels <- sort(unique(hierarchy[hierarchy$cell_type_l1 == l1, level_col]))
      palette <- c(palette, make_shade_palette(l1_palette[[l1]], labels))
    }
    palette
  }

  list(
    cima_cell_type_l1 = l1_palette,
    cima_cell_type_l2 = build_level_palette("cell_type_l2"),
    cima_cell_type_l3 = build_level_palette("cell_type_l3"),
    cima_cell_type_l4 = build_level_palette("cell_type_l4")
  )
}

normalize_embedding_rows <- function(mat) {
  if (nrow(mat) == 0) {
    return(mat)
  }
  norms <- sqrt(rowSums(mat ^ 2))
  norms[norms == 0] <- 1
  sweep(mat, 1, norms, "/")
}

predict_cima_centroids <- function(query_embeddings, centroids, allowed_labels = NULL) {
  labels <- centroids$labels
  centroid_embeddings <- centroids$embedding
  dims_use <- seq.int(2, min(ncol(query_embeddings), ncol(centroid_embeddings)))
  if (length(dims_use) < 2) {
    stop("CIMA centroid prediction requires at least two embedding dimensions")
  }

  query_norm <- normalize_embedding_rows(query_embeddings[, dims_use, drop = FALSE])
  centroid_norm <- normalize_embedding_rows(centroid_embeddings[, dims_use, drop = FALSE])

  pick_from_similarity <- function(similarity, candidate_labels) {
    order_idx <- order(similarity, decreasing = TRUE)
    top_idx <- order_idx[1]
    top_score <- similarity[top_idx]
    margin <- if (length(order_idx) > 1) {
      top_score - similarity[order_idx[2]]
    } else {
      NA_real_
    }
    list(label = candidate_labels[top_idx], score = top_score, margin = margin)
  }

  if (is.null(allowed_labels)) {
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
    return(list(
      label = labels[top_index],
      score = as.numeric(top_score),
      margin = as.numeric(score_margin)
    ))
  }

  label_out <- character(nrow(query_norm))
  score_out <- numeric(nrow(query_norm))
  margin_out <- numeric(nrow(query_norm))
  margin_out[] <- NA_real_

  label_pos <- stats::setNames(seq_along(labels), labels)
  full_similarity <- query_norm %*% t(centroid_norm)
  for (i in seq_len(nrow(query_norm))) {
    candidates <- allowed_labels[[i]]
    if (length(candidates) == 0) {
      candidates <- labels
    }
    pos <- unname(label_pos[candidates])
    pos <- pos[!is.na(pos)]
    if (length(pos) == 0) {
      pos <- seq_along(labels)
    }
    picked <- pick_from_similarity(full_similarity[i, pos], labels[pos])
    label_out[i] <- picked$label
    score_out[i] <- picked$score
    margin_out[i] <- picked$margin
  }

  list(label = label_out, score = score_out, margin = margin_out)
}

annotate_cima_reference_model <- function(counts, hierarchy, feature_model, centroids_by_level) {
  if (ncol(counts) == 0) {
    return(list(
      annotations = data.frame(
        cell_barcode = character(),
        cima_cell_type_l1 = character(),
        cima_cell_type_l2 = character(),
        cima_cell_type_l3 = character(),
        cima_cell_type_l4 = character(),
        cima_l4_score = numeric(),
        cima_l4_score_margin = numeric(),
        stringsAsFactors = FALSE
      ),
      embeddings = matrix(numeric(), nrow = 0, ncol = 0)
    ))
  }

  selected_idx <- feature_model$feature_index
  if (max(selected_idx) > nrow(counts)) {
    stop("CIMA reference feature model expects more peaks than provided in query matrix")
  }
  expected_feature_ids <- feature_model$feature_id
  observed_feature_ids <- rownames(counts)[selected_idx]
  if (!identical(observed_feature_ids, expected_feature_ids)) {
    mismatch_idx <- which(observed_feature_ids != expected_feature_ids)[1]
    stop(
      "CIMA reference feature identity mismatch at selected feature ",
      mismatch_idx,
      ": observed=",
      observed_feature_ids[[mismatch_idx]],
      " expected=",
      expected_feature_ids[[mismatch_idx]]
    )
  }
  selected_counts <- counts[selected_idx, , drop = FALSE]

  binary_counts <- selected_counts
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

  l2_by_l1 <- lapply(split(hierarchy$cell_type_l2, hierarchy$cell_type_l1), unique)
  l3_by_l2 <- lapply(
    split(hierarchy$cell_type_l3, paste(hierarchy$cell_type_l1, hierarchy$cell_type_l2, sep = "||")),
    unique
  )
  l4_by_l3 <- lapply(
    split(hierarchy$cell_type_l4, paste(hierarchy$cell_type_l1, hierarchy$cell_type_l2, hierarchy$cell_type_l3, sep = "||")),
    unique
  )

  pred_l1 <- predict_cima_centroids(query_embeddings, centroids_by_level$l1)
  allowed_l2 <- lapply(pred_l1$label, function(label) l2_by_l1[[label]])
  pred_l2 <- predict_cima_centroids(query_embeddings, centroids_by_level$l2, allowed_labels = allowed_l2)

  allowed_l3 <- lapply(seq_along(pred_l1$label), function(i) {
    key <- paste(pred_l1$label[[i]], pred_l2$label[[i]], sep = "||")
    l3_by_l2[[key]]
  })
  pred_l3 <- predict_cima_centroids(query_embeddings, centroids_by_level$l3, allowed_labels = allowed_l3)

  allowed_l4 <- lapply(seq_along(pred_l1$label), function(i) {
    key <- paste(pred_l1$label[[i]], pred_l2$label[[i]], pred_l3$label[[i]], sep = "||")
    l4_by_l3[[key]]
  })
  pred_l4 <- predict_cima_centroids(query_embeddings, centroids_by_level$l4, allowed_labels = allowed_l4)

  list(
    annotations = data.frame(
      cell_barcode = colnames(binary_counts),
      cima_cell_type_l1 = pred_l1$label,
      cima_cell_type_l2 = pred_l2$label,
      cima_cell_type_l3 = pred_l3$label,
      cima_cell_type_l4 = pred_l4$label,
      cima_l4_score = as.numeric(pred_l4$score),
      cima_l4_score_margin = as.numeric(pred_l4$margin),
      stringsAsFactors = FALSE,
      row.names = colnames(binary_counts)
    ),
    embeddings = query_embeddings
  )
}

run_cima_reference_umap <- function(atac_obj, query_embeddings) {
  available_dims <- ncol(query_embeddings)
  dims_use <- seq.int(2, min(30, available_dims))
  if (length(dims_use) < 2) {
    stop("Not enough reference-projected dimensions available for CIMA UMAP")
  }

  reduction_embeddings <- query_embeddings
  colnames(reduction_embeddings) <- paste0("CIMAREFLSI_", seq_len(ncol(reduction_embeddings)))
  atac_obj[["cima_ref_lsi"]] <- CreateDimReducObject(
    embeddings = reduction_embeddings,
    key = "CIMAREFLSI_",
    assay = DefaultAssay(atac_obj)
  )

  atac_obj <- RunUMAP(
    atac_obj,
    reduction = "cima_ref_lsi",
    dims = dims_use,
    reduction.name = "cima_ref_umap",
    reduction.key = "CIMAUMAP_",
    verbose = FALSE
  )

  list(object = atac_obj, dims = dims_use)
}

run_single_sample_umap <- function(atac_obj) {
  DefaultAssay(atac_obj) <- "ATAC"
  atac_obj <- RunTFIDF(atac_obj)
  atac_obj <- FindTopFeatures(atac_obj, min.cutoff = "q0")
  atac_obj <- RunSVD(atac_obj)

  lsi_embeddings <- Embeddings(atac_obj[["lsi"]])
  available_dims <- ncol(lsi_embeddings)
  dims_use <- seq.int(2, min(30, available_dims))
  if (length(dims_use) < 2) {
    stop("Not enough LSI dimensions available for per-sample UMAP")
  }

  atac_obj <- RunUMAP(
    atac_obj,
    reduction = "lsi",
    dims = dims_use,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )

  list(object = atac_obj, dims = dims_use)
}

save_cima_umap_plot <- function(plot_df, label_col, palette, out_path, title_text) {
  clean_df <- plot_df[!is.na(plot_df[[label_col]]) & nzchar(plot_df[[label_col]]), , drop = FALSE]
  plot_obj <- if (nrow(clean_df) == 0) {
    make_placeholder_plot(title_text, paste0("No values available for ", label_col))
  } else {
    clean_df[[label_col]] <- factor(clean_df[[label_col]], levels = names(palette))
    ggplot(clean_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[label_col]])) +
      geom_point(size = 0.25, alpha = 0.85) +
      scale_color_manual(values = palette, drop = FALSE) +
      labs(title = title_text, x = "UMAP_1", y = "UMAP_2", color = NULL) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.key.height = grid::unit(0.35, "cm"),
        legend.text = element_text(size = 8)
      )
  }

  ggsave(
    filename = out_path,
    plot = plot_obj,
    width = 10,
    height = 8,
    units = "in",
    dpi = 150,
    limitsize = FALSE
  )
}

ensure_tabix_index <- function(fragment_file) {
  is_bgzf_file <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    header <- readBin(con, what = "raw", n = 18)
    if (length(header) < 18) {
      return(FALSE)
    }
    gzip_magic <- identical(as.integer(header[1:2]), c(31L, 139L))
    has_fextra <- bitwAnd(as.integer(header[4]), 4L) != 0L
    has_bgzf_tag <- identical(as.integer(header[13:16]), c(66L, 67L, 2L, 0L))
    gzip_magic && has_fextra && has_bgzf_tag
  }

  ensure_bgzf <- function(path) {
    if (is_bgzf_file(path)) {
      return(path)
    }

    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop(
        "Fragment file is not bgzip-compressed and Rsamtools is unavailable for conversion: ",
        path
      )
    }

    tmp_bgzf <- paste0(path, ".bgzip_tmp")
    if (file.exists(tmp_bgzf)) {
      file.remove(tmp_bgzf)
    }

    cat("Fragment file is regular gzip, converting to bgzip in place:", basename(path), "\n")
    Rsamtools::bgzip(path, dest = tmp_bgzf, overwrite = TRUE)

    if (!file.exists(tmp_bgzf)) {
      stop("bgzip conversion did not produce: ", tmp_bgzf)
    }

    if (file.exists(path) && !file.remove(path)) {
      stop("Failed to remove original gzip fragment file before replacing with bgzip: ", path)
    }
    if (!file.rename(tmp_bgzf, path)) {
      stop("Failed to replace fragment file with bgzip version: ", path)
    }
    path
  }

  fragment_file <- ensure_bgzf(fragment_file)
  tbi_file <- paste0(fragment_file, ".tbi")
  if (file.exists(tbi_file)) {
    return(tbi_file)
  }

  cat("Tabix index missing, creating:", basename(tbi_file), "\n")

  tabix_bin <- Sys.which("tabix")
  if (nzchar(tabix_bin)) {
    output <- system2(
      tabix_bin,
      args = c("-p", "bed", fragment_file),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0
    }
    if (status != 0) {
      stop(
        "Failed to create tabix index with system tabix: ",
        paste(output, collapse = "\n")
      )
    }
  } else if (requireNamespace("Rsamtools", quietly = TRUE)) {
    Rsamtools::indexTabix(fragment_file, format = "bed")
  } else {
    stop(
      "Missing tabix index and no indexer available. Install tabix or Rsamtools."
    )
  }

  if (!file.exists(tbi_file)) {
    stop("Tabix index creation did not produce: ", tbi_file)
  }

  tbi_file
}

option_list <- list(
  make_option(c("--gse"), type = "character", help = "GSE accession"),
  make_option(c("--gsm"), type = "character", help = "GSM accession"),
  make_option(c("--nmads"), type = "numeric", default = 4, help = "MAD multiplier"),
  make_option(c("--output-profile"), type = "character", default = "full", help = "Output profile: full or matrix-lite"),
  make_option(c("--output-root"), type = "character", default = file.path(project_root, "output"), help = "Output root directory")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gse) || is.null(opt$gsm)) {
  print_help(opt_parser)
  stop("Both --gse and --gsm are required")
}

valid_output_profiles <- c("full", "matrix-lite", "validation-lite")
if (!opt$`output-profile` %in% valid_output_profiles) {
  stop("Unsupported --output-profile: ", opt$`output-profile`, ". Choose one of: ", paste(valid_output_profiles, collapse = ", "))
}
is_lite_output <- opt$`output-profile` %in% c("matrix-lite", "validation-lite")

raw_dir <- file.path(project_root, "data", "raw", opt$gse)
peak_file <- file.path(project_root, "data", "reference", "peak.bed")
cima_dir <- file.path(project_root, "data", "reference", "cima")
cima_reference_file <- file.path(cima_dir, "CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad")
cima_hierarchy_file <- file.path(cima_dir, "cima_atac_celltype_hierarchy.csv")
cima_feature_model_file <- file.path(cima_dir, "cima_atac_reference_lsi_features.tsv.gz")
cima_l1_centroid_file <- file.path(cima_dir, "cima_atac_reference_l1_centroids.tsv")
cima_l2_centroid_file <- file.path(cima_dir, "cima_atac_reference_l2_centroids.tsv")
cima_l3_centroid_file <- file.path(cima_dir, "cima_atac_reference_l3_centroids.tsv")
cima_l4_centroid_file <- file.path(cima_dir, "cima_atac_reference_l4_centroids.tsv")
output_dir <- file.path(opt$`output-root`, opt$gse, opt$gsm)
matrix_dir <- file.path(output_dir, "matrix")

if (!dir.exists(raw_dir)) {
  stop("Raw directory not found: ", raw_dir)
}
if (!file.exists(peak_file)) {
  stop("Peak file not found: ", peak_file)
}
if (!file.exists(cima_hierarchy_file)) {
  stop("CIMA ATAC hierarchy file not found: ", cima_hierarchy_file)
}
required_cima_model_files <- c(
  cima_feature_model_file,
  cima_l1_centroid_file,
  cima_l2_centroid_file,
  cima_l3_centroid_file,
  cima_l4_centroid_file
)
missing_cima_model_files <- required_cima_model_files[!file.exists(required_cima_model_files)]
if (length(missing_cima_model_files) > 0) {
  stop("Missing CIMA reference model files: ", paste(missing_cima_model_files, collapse = ", "))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_gz <- function(dataframe, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(dataframe, con, row.names = FALSE)
}

write_lines_plain <- function(lines, out_path) {
  con <- file(out_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, useBytes = TRUE)
}

make_barcode_rank_plot <- function(fragment_stats, knee_threshold, inflection_threshold, floor_auto, candidate_threshold, candidate_rank_cap) {
  plot_df <- fragment_stats
  threshold_df <- combine_line_df(
    make_line_df("T_knee", knee_threshold),
    make_line_df("T_inflection", inflection_threshold),
    make_line_df("T_floor_auto", floor_auto),
    make_line_df("T_candidate", candidate_threshold)
  )
  rank_df <- combine_line_df(
    make_line_df("N_candidate_auto", candidate_rank_cap)
  )

  ggplot(plot_df, aes(x = rank, y = frequency_count)) +
    geom_line(linewidth = 0.4, color = "#2f4f4f") +
    geom_hline(
      data = threshold_df,
      aes(yintercept = value, color = name),
      linewidth = 0.7,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_vline(
      data = rank_df,
      aes(xintercept = value),
      linewidth = 0.8,
      linetype = "dotdash",
      color = "#5e3c99",
      inherit.aes = FALSE
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(
      values = c(
        T_knee = "#d55e00",
        T_inflection = "#cc79a7",
        T_floor_auto = "#0072b2",
        T_candidate = "#009e73"
      )
    ) +
    labs(
      title = "Barcode Rank",
      x = "Barcode rank",
      y = "Fragment count",
      color = NULL
    ) +
    theme_bw(base_size = 11)
}

make_barcode_density_plot <- function(fragment_stats, floor_auto, candidate_threshold) {
  plot_df <- data.frame(log_count = log10(fragment_stats$frequency_count + 1))
  line_df <- combine_line_df(
    make_line_df("T_floor_auto", log10(floor_auto + 1)),
    make_line_df("T_candidate", log10(candidate_threshold + 1))
  )

  ggplot(plot_df, aes(x = log_count)) +
    geom_density(fill = "#9ecae1", alpha = 0.6, color = "#08519c", linewidth = 0.6) +
    geom_vline(
      data = line_df,
      aes(xintercept = value, color = name),
      linewidth = 0.8,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        T_floor_auto = "#0072b2",
        T_candidate = "#009e73"
      )
    ) +
    labs(
      title = "Barcode Count Density",
      x = "log10(fragment_count + 1)",
      y = "Density",
      color = NULL
    ) +
    theme_bw(base_size = 11)
}

fragment_file <- find_single_file(raw_dir, opt$gsm, "fragments", required = TRUE)
barcode_file <- find_single_file(raw_dir, opt$gsm, "filtered_barcodes", required = FALSE)
singlecell_file <- find_single_file(raw_dir, opt$gsm, "singlecell", required = FALSE, exts = c("csv", "tsv"))

cat(rep("=", 80), "\n", sep = "")
cat("scATAC-seq Single Sample QC\n")
cat(rep("=", 80), "\n", sep = "")
cat("Project root:", project_root, "\n")
cat("GSE ID:", opt$gse, "\n")
cat("GSM ID:", opt$gsm, "\n")
cat("Fragment file:", fragment_file, "\n")
cat("Peak file:", peak_file, "\n")
cat("Output directory:", output_dir, "\n")

barcode_source <- "filtered_barcodes"
barcode_knee_threshold <- NA_real_
barcode_inflection_threshold <- NA_real_
barcode_floor_auto <- NA_real_
barcode_candidate_threshold <- NA_real_
barcode_candidate_rank_cap <- NA_real_
barcode_prefilter_stats <- NULL
singlecell_metrics <- NULL
barcode_rank_plot <- make_placeholder_plot("Barcode Rank", "Using provided filtered_barcodes")
barcode_density_plot <- make_placeholder_plot("Barcode Count Density", "Using provided filtered_barcodes")

if (!is.null(barcode_file)) {
  cat("Barcode file:", barcode_file, "\n\n")
  barcodes <- readLines(gzfile(barcode_file))
  barcodes <- unique(barcodes[nzchar(barcodes)])
} else if (!is.null(singlecell_file)) {
  cat("Barcode file: not found\n")
  cat("singlecell file:", singlecell_file, "\n\n")
  singlecell <- load_singlecell_barcodes(singlecell_file)
  barcodes <- singlecell$barcodes
  singlecell_metrics <- singlecell$metrics
  barcode_source <- "singlecell_csv"
  barcode_rank_plot <- make_placeholder_plot("Barcode Rank", paste0("Using ", basename(singlecell_file)))
  barcode_density_plot <- make_placeholder_plot("Barcode Count Density", paste0("Using ", singlecell$selector_source))
} else {
  cat("Barcode file: not found, inferring initial barcode set from fragment counts\n\n")
  inferred <- infer_barcodes(fragment_file)
  barcodes <- inferred$barcodes
  barcode_source <- "prefiltered_from_fragments"
  barcode_knee_threshold <- inferred$knee_threshold
  barcode_inflection_threshold <- inferred$inflection_threshold
  barcode_floor_auto <- inferred$floor_auto
  barcode_candidate_threshold <- inferred$candidate_threshold
  barcode_candidate_rank_cap <- inferred$candidate_rank_cap
  barcode_prefilter_stats <- inferred$fragment_stats
  barcode_rank_plot <- make_barcode_rank_plot(
    barcode_prefilter_stats,
    barcode_knee_threshold,
    barcode_inflection_threshold,
    barcode_floor_auto,
    barcode_candidate_threshold,
    barcode_candidate_rank_cap
  )
  barcode_density_plot <- make_barcode_density_plot(
    barcode_prefilter_stats,
    barcode_floor_auto,
    barcode_candidate_threshold
  )
}

if (length(barcodes) == 0) {
  stop("No barcodes selected for sample ", opt$gsm)
}

cat("Initial barcodes:", length(barcodes), "(", barcode_source, ")\n\n")

cat("[1/7] Loading peaks and annotations...\n")
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
project_peak_ids <- make_peak_ids(peaks_gr)

cima_hierarchy <- load_cima_hierarchy(cima_hierarchy_file)
cima_feature_model <- load_cima_reference_feature_model(cima_feature_model_file)
cima_centroids <- list(
  l1 = load_cima_reference_centroids(cima_l1_centroid_file, "cell_type_l1"),
  l2 = load_cima_reference_centroids(cima_l2_centroid_file, "cell_type_l2"),
  l3 = load_cima_reference_centroids(cima_l3_centroid_file, "cell_type_l3"),
  l4 = load_cima_reference_centroids(cima_l4_centroid_file, "cell_type_l4")
)
if (max(cima_feature_model$feature_index) > length(project_peak_ids)) {
  stop(
    "CIMA reference feature model exceeds loaded project peak space: max index=",
    max(cima_feature_model$feature_index),
    " peaks=",
    length(project_peak_ids)
  )
}

cat("Peaks:", length(peaks_gr), "\n")
cat("Common chromosomes:", length(common_seq), "\n\n")

cat("[2/7] Creating fragment object...\n")
tbi_file <- ensure_tabix_index(fragment_file)

frags <- CreateFragmentObject(path = fragment_file, cells = barcodes)
cat("Fragment object ready\n\n")

cat("[3/7] Building peak-by-cell matrix...\n")
counts <- FeatureMatrix(
  fragments = frags,
  features = peaks_gr,
  cells = barcodes
)
if (nrow(counts) != length(project_peak_ids)) {
  stop(
    "FeatureMatrix row count does not match project peak reference order: matrix rows=",
    nrow(counts),
    " expected peaks=",
    length(project_peak_ids)
  )
}
rownames(counts) <- project_peak_ids
cat("Matrix shape:", nrow(counts), "peaks x", ncol(counts), "cells\n\n")

cat("[4/7] Creating Seurat object...\n")
atac_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = frags,
  genome = "hg38",
  sep = c(":", "-")
)
slot(atac_assay, "annotation") <- annotations

atac_obj <- CreateSeuratObject(
  counts = atac_assay,
  assay = "ATAC",
  project = opt$gsm
)

atac_obj$sample <- opt$gsm
atac_obj$dataset <- opt$gse
atac_obj$barcode_source <- barcode_source

cat("Cells:", ncol(atac_obj), "\n")
cat("Features:", nrow(atac_obj), "\n\n")

cat("[5/7] Computing QC metrics...\n")
atac_obj <- NucleosomeSignal(atac_obj)
atac_obj <- TSSEnrichment(atac_obj, fast = FALSE)

if (!is.null(barcode_prefilter_stats)) {
  fragment_stats <- barcode_prefilter_stats[barcode_prefilter_stats$CB %in% colnames(atac_obj), , drop = FALSE]
} else {
  fragment_stats <- tryCatch(
    CountFragments(fragment_file, cells = colnames(atac_obj)),
    error = function(e) {
      warning(
        "CountFragments failed for ",
        basename(fragment_file),
        ": ",
        conditionMessage(e),
        ". Falling back to nCount_ATAC for fragment totals."
      )
      fallback <- data.frame(
        CB = colnames(atac_obj),
        frequency_count = as.numeric(atac_obj$nCount_ATAC),
        stringsAsFactors = FALSE
      )
      fallback
    }
  )
}
rownames(fragment_stats) <- fragment_stats$CB
atac_obj$fragments <- fragment_stats[colnames(atac_obj), "frequency_count"]
atac_obj$total_fragments <- atac_obj$fragments
if ("reads_count" %in% colnames(fragment_stats)) {
  atac_obj$reads_count <- fragment_stats[colnames(atac_obj), "reads_count"]
  atac_obj$unique_ratio <- ifelse(
    atac_obj$reads_count > 0,
    atac_obj$fragments / atac_obj$reads_count,
    NA_real_
  )
}

if (!is.null(singlecell_metrics) && nrow(singlecell_metrics) > 0) {
  rownames(singlecell_metrics) <- singlecell_metrics$barcode
  common_metrics <- intersect(rownames(singlecell_metrics), colnames(atac_obj))
  singlecell_cols <- setdiff(colnames(singlecell_metrics), "barcode")
  for (metric_name in singlecell_cols) {
    meta_col <- paste0("singlecell_", metric_name)
    atac_obj@meta.data[, meta_col] <- NA
    atac_obj@meta.data[common_metrics, meta_col] <- singlecell_metrics[common_metrics, metric_name]
  }
}
atac_obj <- FRiP(object = atac_obj, assay = "ATAC", total.fragments = "fragments")
atac_obj$blacklist_fraction <- FractionCountsInRegion(
  object = atac_obj,
  assay = "ATAC",
  regions = blacklist_hg38_unified
)
cat("QC metrics ready\n\n")

cat("[6/7] Running doublet detection...\n")
counts_matrix <- GetAssayData(atac_obj[["ATAC"]], layer = "counts")
sce <- SingleCellExperiment(list(counts = counts_matrix))
sce <- scDblFinder(
  sce,
  clusters = TRUE,
  aggregateFeatures = TRUE,
  nfeatures = 25,
  processing = "normFeatures"
)
doublet_info <- as.data.frame(colData(sce)[, c("scDblFinder.class", "scDblFinder.score")])
atac_obj$scDblFinder.class <- doublet_info$scDblFinder.class
atac_obj$scDblFinder.score <- doublet_info$scDblFinder.score
cat("Singlets:", sum(doublet_info$scDblFinder.class == "singlet", na.rm = TRUE), "\n")
cat("Doublets:", sum(doublet_info$scDblFinder.class == "doublet", na.rm = TRUE), "\n\n")

cat("[7/7] Applying QC filters...\n")
meta_data <- atac_obj@meta.data
meta_data <- subset(meta_data, nCount_ATAC > 0)

outlier_count <- is_outlier(meta_data, "nCount_ATAC", opt$nmads)
outlier_tss <- is_outlier(meta_data, "TSS.enrichment", opt$nmads)
outlier_frip <- if ("FRiP" %in% colnames(meta_data)) {
  is_outlier(meta_data, "FRiP", opt$nmads)
} else {
  rep(FALSE, nrow(meta_data))
}

mad_outlier <- outlier_count | outlier_tss | outlier_frip
is_doublet <- meta_data$scDblFinder.class == "doublet"
meta_data$outlier <- mad_outlier | is_doublet
meta_data$pass_qc <- !meta_data$outlier
atac_obj@meta.data <- meta_data

qc_cells <- rownames(meta_data)[meta_data$pass_qc %in% TRUE]
qc_counts <- GetAssayData(atac_obj[["ATAC"]], layer = "counts")[, qc_cells, drop = FALSE]
rownames(qc_counts) <- project_peak_ids

meta_data$cima_cell_type_l1 <- NA_character_
meta_data$cima_cell_type_l2 <- NA_character_
meta_data$cima_cell_type_l3 <- NA_character_
meta_data$cima_cell_type_l4 <- NA_character_
meta_data$cima_l4_score <- NA_real_
meta_data$cima_l4_score_margin <- NA_real_
meta_data$umap_atac_1 <- NA_real_
meta_data$umap_atac_2 <- NA_real_
meta_data$cima_ref_umap_1 <- NA_real_
meta_data$cima_ref_umap_2 <- NA_real_
atac_obj@meta.data <- meta_data

cat("Cells before QC:", nrow(meta_data), "\n")
cat("Cells after QC:", length(qc_cells), "\n")
cat("QC rate:", sprintf("%.2f%%", length(qc_cells) / nrow(meta_data) * 100), "\n\n")

umap_l1_file <- file.path(output_dir, "umap_cima_cell_type_l1.png")
umap_l2_file <- file.path(output_dir, "umap_cima_cell_type_l2.png")
umap_l3_file <- file.path(output_dir, "umap_cima_cell_type_l3.png")
umap_l4_file <- file.path(output_dir, "umap_cima_cell_type_l4.png")

cat("[8/10] Assigning CIMA labels...\n")
annotation_result <- annotate_cima_reference_model(qc_counts, cima_hierarchy, cima_feature_model, cima_centroids)
annotation_df <- annotation_result$annotations
query_cima_embeddings <- annotation_result$embeddings
annotation_cols <- c(
  "cima_cell_type_l1",
  "cima_cell_type_l2",
  "cima_cell_type_l3",
  "cima_cell_type_l4",
  "cima_l4_score",
  "cima_l4_score_margin"
)
atac_obj@meta.data[annotation_df$cell_barcode, annotation_cols] <- annotation_df[, annotation_cols]
cat("Assigned L4 labels:", length(unique(annotation_df$cima_cell_type_l4)), "\n\n")

cat("[9/10] Running per-sample UMAP...\n")
umap_plot_df <- NULL
cima_palettes <- build_cima_palettes(cima_hierarchy)
cima_umap_basis <- ""

if (length(qc_cells) >= 3) {
  atac_qc_obj <- subset(atac_obj, cells = qc_cells)
  umap_result <- tryCatch(
    run_single_sample_umap(atac_qc_obj),
    error = function(e) {
      message("UMAP failed for ", opt$gsm, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(umap_result)) {
    atac_qc_obj <- umap_result$object
    umap_embeddings <- Embeddings(atac_qc_obj[["umap"]])
    atac_obj@meta.data[rownames(umap_embeddings), c("umap_atac_1", "umap_atac_2")] <- umap_embeddings[, 1:2, drop = FALSE]
    umap_plot_df <- data.frame(
      cell_barcode = rownames(umap_embeddings),
      UMAP_1 = umap_embeddings[, 1],
      UMAP_2 = umap_embeddings[, 2],
      atac_qc_obj@meta.data[, c("cima_cell_type_l1", "cima_cell_type_l2", "cima_cell_type_l3", "cima_cell_type_l4"), drop = FALSE],
      stringsAsFactors = FALSE,
      row.names = rownames(umap_embeddings)
    )
  }

  cima_ref_result <- tryCatch(
    run_cima_reference_umap(atac_qc_obj, query_cima_embeddings[qc_cells, , drop = FALSE]),
    error = function(e) {
      message("CIMA reference-space UMAP failed for ", opt$gsm, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(cima_ref_result)) {
    atac_cima_obj <- cima_ref_result$object
    cima_ref_embeddings <- Embeddings(atac_cima_obj[["cima_ref_umap"]])
    atac_obj@meta.data[rownames(cima_ref_embeddings), c("cima_ref_umap_1", "cima_ref_umap_2")] <- cima_ref_embeddings[, 1:2, drop = FALSE]
    umap_plot_df <- data.frame(
      cell_barcode = rownames(cima_ref_embeddings),
      UMAP_1 = cima_ref_embeddings[, 1],
      UMAP_2 = cima_ref_embeddings[, 2],
      atac_cima_obj@meta.data[, c("cima_cell_type_l1", "cima_cell_type_l2", "cima_cell_type_l3", "cima_cell_type_l4"), drop = FALSE],
      stringsAsFactors = FALSE,
      row.names = rownames(cima_ref_embeddings)
    )
    cima_umap_basis <- "reference_projected_lsi"
  } else if (!is.null(umap_plot_df)) {
    cima_umap_basis <- "query_native_lsi"
  }
}

if (is.null(umap_plot_df)) {
  placeholder_df <- data.frame(
    UMAP_1 = numeric(),
    UMAP_2 = numeric(),
    cima_cell_type_l1 = character(),
    cima_cell_type_l2 = character(),
    cima_cell_type_l3 = character(),
    cima_cell_type_l4 = character(),
    stringsAsFactors = FALSE
  )
  umap_plot_df <- placeholder_df
  cima_umap_basis <- ""
}

save_cima_umap_plot(umap_plot_df, "cima_cell_type_l1", cima_palettes$cima_cell_type_l1, umap_l1_file, paste0(opt$gsm, " UMAP by CIMA L1"))
save_cima_umap_plot(umap_plot_df, "cima_cell_type_l2", cima_palettes$cima_cell_type_l2, umap_l2_file, paste0(opt$gsm, " UMAP by CIMA L2"))
save_cima_umap_plot(umap_plot_df, "cima_cell_type_l3", cima_palettes$cima_cell_type_l3, umap_l3_file, paste0(opt$gsm, " UMAP by CIMA L3"))
save_cima_umap_plot(umap_plot_df, "cima_cell_type_l4", cima_palettes$cima_cell_type_l4, umap_l4_file, paste0(opt$gsm, " UMAP by CIMA L4"))

cat("[10/10] Building QC overview figure...\n")
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
  title = paste0(opt$gsm, " QC Overview"),
  subtitle = paste0(
    opt$gse,
    " | initial cells: ", ncol(atac_obj),
    " | pass QC: ", length(qc_cells),
    " | doublets: ", sum(atac_obj$scDblFinder.class == "doublet", na.rm = TRUE)
  )
)

qc_plot_file <- file.path(output_dir, "qc_overview.png")
ggsave(
  filename = qc_plot_file,
  plot = qc_overview,
  width = 24,
  height = 32,
  units = "in",
  dpi = 150,
  limitsize = FALSE
)

meta_data <- atac_obj@meta.data
metadata_all <- cbind(cell_barcode = rownames(meta_data), meta_data)
metadata_qc <- metadata_all[metadata_all$pass_qc %in% TRUE, , drop = FALSE]
median_unique_ratio <- if ("unique_ratio" %in% colnames(atac_obj@meta.data)) {
  sprintf("%.2f", median(atac_obj$unique_ratio[atac_obj$pass_qc], na.rm = TRUE))
} else {
  ""
}
median_cima_l4_score <- if (all(is.na(metadata_qc$cima_l4_score))) {
  ""
} else {
  sprintf("%.4f", median(metadata_qc$cima_l4_score, na.rm = TRUE))
}
unique_cima_l4 <- length(unique(metadata_qc$cima_cell_type_l4[!is.na(metadata_qc$cima_cell_type_l4)]))

summary_stats <- data.frame(
  metric = c(
    "gse",
    "gsm",
    "fragment_file",
    "barcode_source",
    "barcode_knee_threshold",
    "barcode_inflection_threshold",
    "barcode_floor_auto",
    "barcode_candidate_threshold",
    "barcode_candidate_rank_cap",
    "input_cells",
    "pass_qc",
    "qc_rate",
    "singlets",
    "doublets",
    "mean_nCount_ATAC",
    "median_TSS_enrichment",
    "median_FRiP",
    "median_fragments",
    "median_unique_ratio",
    "median_blacklist_fraction",
    "cima_reference_atac_h5ad",
    "cima_reference_model_features",
    "cima_annotation_levels",
    "cima_umap_basis",
    "cima_unique_l4_labels",
    "median_cima_l4_score"
  ),
  value = c(
    opt$gse,
    opt$gsm,
    fragment_file,
    barcode_source,
    if (is.na(barcode_knee_threshold)) "" else sprintf("%.0f", barcode_knee_threshold),
    if (is.na(barcode_inflection_threshold)) "" else sprintf("%.0f", barcode_inflection_threshold),
    if (is.na(barcode_floor_auto)) "" else sprintf("%.0f", barcode_floor_auto),
    if (is.na(barcode_candidate_threshold)) "" else sprintf("%.0f", barcode_candidate_threshold),
    if (is.na(barcode_candidate_rank_cap)) "" else sprintf("%.0f", barcode_candidate_rank_cap),
    ncol(atac_obj),
    length(qc_cells),
    sprintf("%.2f%%", length(qc_cells) / ncol(atac_obj) * 100),
    sum(atac_obj$scDblFinder.class == "singlet", na.rm = TRUE),
    sum(atac_obj$scDblFinder.class == "doublet", na.rm = TRUE),
    sprintf("%.0f", mean(atac_obj$nCount_ATAC[atac_obj$pass_qc], na.rm = TRUE)),
    sprintf("%.2f", median(atac_obj$TSS.enrichment[atac_obj$pass_qc], na.rm = TRUE)),
    sprintf("%.4f", median(atac_obj$FRiP[atac_obj$pass_qc], na.rm = TRUE)),
    sprintf("%.0f", median(atac_obj$fragments[atac_obj$pass_qc], na.rm = TRUE)),
    median_unique_ratio,
    sprintf("%.4f", median(atac_obj$blacklist_fraction[atac_obj$pass_qc], na.rm = TRUE)),
    cima_reference_file,
    cima_feature_model_file,
    "l1,l2,l3,l4",
    cima_umap_basis,
    unique_cima_l4,
    median_cima_l4_score
  ),
  stringsAsFactors = FALSE
)

feature_file <- if (is_lite_output) file.path(matrix_dir, "features.tsv") else file.path(matrix_dir, "features.tsv.gz")
barcode_out_file <- if (is_lite_output) file.path(matrix_dir, "barcodes.tsv") else file.path(matrix_dir, "barcodes.tsv.gz")
matrix_file <- file.path(matrix_dir, "matrix.mtx")
metadata_file <- file.path(output_dir, "metadata.csv")
metadata_qc_file <- file.path(output_dir, "metadata_qc.csv")
summary_file <- file.path(output_dir, "qc_summary.csv")
rds_file <- file.path(output_dir, paste0(opt$gsm, "_seurat_qc.rds"))
validation_file <- file.path(output_dir, "validation_result.csv")

validation_keep <- intersect(
  c(
    "cell_barcode",
    "seurat_clusters",
    "nCount_ATAC",
    "nFeature_ATAC",
    "TSS.enrichment",
    "nucleosome_signal",
    "FRiP",
    "blacklist_fraction",
    "scDblFinder.class",
    "cima_cell_type_l1",
    "cima_cell_type_l2",
    "cima_cell_type_l3",
    "cima_cell_type_l4",
    "cima_l4_score",
    "cima_l4_score_margin",
    "umap_atac_1",
    "umap_atac_2",
    "cima_ref_umap_1",
    "cima_ref_umap_2"
  ),
  colnames(metadata_qc)
)
validation_result <- metadata_qc[, validation_keep, drop = FALSE]

write.csv(summary_stats, summary_file, row.names = FALSE)
write.csv(validation_result, validation_file, row.names = FALSE)
writeMM(qc_counts, matrix_file)
if (is_lite_output) {
  write_lines_plain(rownames(qc_counts), feature_file)
  write_lines_plain(colnames(qc_counts), barcode_out_file)
} else {
  write_lines_gz(rownames(qc_counts), feature_file)
  write_lines_gz(colnames(qc_counts), barcode_out_file)
}

if (!is_lite_output) {
  write.csv(metadata_all, metadata_file, row.names = FALSE)
  write.csv(metadata_qc, metadata_qc_file, row.names = FALSE)
  saveRDS(atac_obj, rds_file)
} else {
  suppressWarnings(file.remove(qc_plot_file, umap_l2_file, umap_l3_file, umap_l4_file))
}

cat("UMAP L1:", umap_l1_file, "\n")
cat("Summary:", summary_file, "\n")
cat("Validation result:", validation_file, "\n")
cat("QC matrix:", matrix_file, "\n")
if (!is_lite_output) {
  cat("QC overview:", qc_plot_file, "\n")
  cat("UMAP L2:", umap_l2_file, "\n")
  cat("UMAP L3:", umap_l3_file, "\n")
  cat("UMAP L4:", umap_l4_file, "\n")
  cat("Metadata:", metadata_file, "\n")
  cat("QC metadata:", metadata_qc_file, "\n")
  cat("RDS:", rds_file, "\n")
}
