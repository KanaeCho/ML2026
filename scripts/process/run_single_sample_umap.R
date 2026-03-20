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
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})

write_csv_auto <- function(dt, path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    fwrite(dt, path)
    return(invisible(path))
  }

  tmp_path <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_path), add = TRUE)
  fwrite(dt, tmp_path)

  in_con <- file(tmp_path, open = "rb")
  out_con <- gzfile(path, open = "wb")
  on.exit(close(in_con), add = TRUE)
  on.exit(close(out_con), add = TRUE)

  repeat {
    chunk <- readBin(in_con, what = raw(), n = 1024 * 1024)
    if (length(chunk) == 0) {
      break
    }
    writeBin(chunk, out_con)
  }

  invisible(path)
}

parse_int_list <- function(start_value, end_value) {
  if (is.na(start_value) || is.na(end_value) || end_value < start_value) {
    stop("Invalid dimension range")
  }
  seq.int(start_value, end_value)
}

format_dims <- function(values) {
  if (length(values) == 0) {
    return("none")
  }
  paste(values, collapse = ",")
}

category_palette <- function(values) {
  values <- unique(as.character(values))
  base_colors <- c(
    "#0f766e", "#b45309", "#2563eb", "#be123c", "#7c3aed", "#0891b2", "#65a30d", "#c2410c",
    "#4f46e5", "#15803d", "#dc2626", "#1d4ed8", "#a21caf", "#0369a1", "#ca8a04", "#374151",
    "#0e7490", "#92400e", "#166534", "#9f1239", "#7f1d1d", "#1f2937", "#854d0e", "#0c4a6e"
  )
  cols <- rep(base_colors, length.out = length(values))
  names(cols) <- values
  cols
}

normalize_string <- function(values) {
  values <- tolower(trimws(as.character(values)))
  values[is.na(values)] <- ""
  gsub("[^a-z0-9]+", "", values)
}

normalize_barcode <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values)] <- ""
  output <- character(length(values))

  for (idx in seq_along(values)) {
    value <- values[idx]
    if (!nzchar(value)) {
      output[idx] <- ""
      next
    }

    upper_value <- toupper(value)
    barcode_match <- regexpr("[ACGTN]+-[0-9]+", upper_value, perl = TRUE)
    if (barcode_match[1] > 0) {
      output[idx] <- regmatches(upper_value, barcode_match)
      next
    }

    pieces <- unlist(strsplit(upper_value, "[#|:;,_/\\\\]+", perl = TRUE), use.names = FALSE)
    pieces <- pieces[nzchar(pieces)]
    piece_match <- grep("^[ACGTN]+-[0-9]+$", pieces, perl = TRUE, value = TRUE)
    if (length(piece_match) > 0) {
      output[idx] <- piece_match[length(piece_match)]
      next
    }

    output[idx] <- gsub("[^A-Z0-9-]+", "", upper_value)
  }

  output
}

choose_existing_column <- function(df, user_value, candidates, required = FALSE) {
  if (!is.null(user_value) && nzchar(user_value)) {
    if (!user_value %in% colnames(df)) {
      stop("Column not found: ", user_value)
    }
    return(user_value)
  }

  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) > 0) {
    return(hits[1])
  }

  if (isTRUE(required)) {
    stop("Unable to resolve required column from candidates: ", paste(candidates, collapse = ", "))
  }
  NULL
}

discover_sample_aliases <- function(project_root, gse, gsm) {
  aliases <- gsm
  raw_dir <- file.path(project_root, "data", "raw", gse)
  if (!dir.exists(raw_dir)) {
    return(unique(normalize_string(aliases)))
  }

  sample_files <- list.files(raw_dir, pattern = paste0("^", gsm, "_"), full.names = FALSE)
  if (length(sample_files) == 0) {
    return(unique(normalize_string(aliases)))
  }

  cleaned <- sub(paste0("^", gsm, "_"), "", sample_files)
  cleaned <- sub("_(fragments|filtered_barcodes|raw_barcodes|barcodes|singlecell).*$", "", cleaned)
  aliases <- c(aliases, cleaned)
  aliases <- c(aliases, unlist(strsplit(cleaned, "[-_]+"), use.names = FALSE))
  aliases <- aliases[nzchar(aliases)]
  unique(normalize_string(aliases))
}

make_dim_plot <- function(object, group_by, title_text, point_size = 0.16) {
  values <- as.character(object[[group_by]][, 1])
  values[is.na(values) | !nzchar(values)] <- "Unlabeled"
  object[[group_by]] <- values
  palette <- category_palette(values)
  DimPlot(
    object = object,
    reduction = "umap_atac",
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
      legend.text = element_text(size = 8)
    )
}

make_feature_plot <- function(object, features, title_text) {
  FeaturePlot(
    object = object,
    reduction = "umap_atac",
    features = features,
    ncol = 2,
    pt.size = 0.16,
    raster = FALSE,
    cols = c("#f8fafc", "#b45309")
  ) &
    ggtitle(title_text) &
    theme_bw(base_size = 12) &
    theme(plot.title = element_text(face = "bold"))
}

save_plot <- function(path, plot, width = 10, height = 8) {
  ggsave(path, plot = plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

resolve_sample_dir <- function(project_root, gse, gsm, sample_dir) {
  if (!is.null(sample_dir) && nzchar(sample_dir)) {
    return(normalizePath(sample_dir, winslash = "/", mustWork = TRUE))
  }
  if (is.null(gse) || !nzchar(gse) || is.null(gsm) || !nzchar(gsm)) {
    stop("Provide either --sample-dir or both --gse and --gsm")
  }
  normalizePath(file.path(project_root, "output", gse, gsm), winslash = "/", mustWork = TRUE)
}

resolve_rds_path <- function(sample_dir, gsm) {
  preferred <- file.path(sample_dir, paste0(gsm, "_seurat_qc.rds"))
  if (file.exists(preferred)) {
    return(preferred)
  }
  rds_files <- list.files(sample_dir, pattern = "_seurat_qc\\.rds$", full.names = TRUE)
  if (length(rds_files) == 1) {
    return(normalizePath(rds_files[1], winslash = "/", mustWork = TRUE))
  }
  stop("Unable to resolve single-sample QC RDS under ", sample_dir)
}

resolve_fragment_path <- function(project_root, gse, gsm) {
  if (is.null(gse) || !nzchar(gse)) {
    return(NULL)
  }
  raw_dir <- file.path(project_root, "data", "raw", gse)
  if (!dir.exists(raw_dir)) {
    return(NULL)
  }
  matches <- list.files(
    raw_dir,
    pattern = paste0("^", gsm, ".*fragments.*tsv\\.gz$"),
    full.names = TRUE
  )
  if (length(matches) == 0) {
    return(NULL)
  }
  normalizePath(matches[1], winslash = "/", mustWork = TRUE)
}

repair_fragment_path <- function(object, project_root, gse, gsm) {
  current_frags <- Fragments(object)
  if (length(current_frags) == 0) {
    return(object)
  }

  current_path <- current_frags[[1]]@path
  if (file.exists(current_path)) {
    return(object)
  }

  replacement_path <- resolve_fragment_path(project_root, gse, gsm)
  if (is.null(replacement_path)) {
    warning("Unable to repair fragment path for ", gsm, "; GeneActivity may fail")
    return(object)
  }

  current_frags[[1]]@path <- replacement_path
  slot(object[[DefaultAssay(object)]], "fragments") <- current_frags
  object
}

prepare_query_object <- function(object, npcs, dims_use, resolution) {
  DefaultAssay(object) <- "ATAC"
  object <- RunTFIDF(object)
  object <- FindTopFeatures(object, min.cutoff = "q0")
  object <- RunSVD(object, n = npcs, reduction.key = "LSI_")

  available_dims <- seq_len(ncol(Embeddings(object[["lsi"]])))
  dims_use <- intersect(dims_use, available_dims)
  if (length(dims_use) < 2) {
    stop("Need at least two LSI dimensions for single-sample UMAP")
  }

  object <- FindNeighbors(object, reduction = "lsi", dims = dims_use, verbose = FALSE)
  object <- FindClusters(object, resolution = resolution, verbose = FALSE)
  object <- RunUMAP(
    object,
    reduction = "lsi",
    dims = dims_use,
    reduction.name = "umap_atac",
    reduction.key = "UMAPATAC_",
    verbose = FALSE
  )

  list(object = object, dims_use = dims_use)
}

collapse_annotation_rows <- function(annotation_df, barcode_col, label_col, subtype_col = NULL) {
  ann <- as.data.table(annotation_df)
  ann[, barcode_norm := normalize_barcode(get(barcode_col))]
  ann[, label_value := as.character(get(label_col))]
  if (!is.null(subtype_col) && subtype_col %in% colnames(ann)) {
    ann[, subtype_value := as.character(get(subtype_col))]
  } else {
    ann[, subtype_value := NA_character_]
  }

  ann <- ann[nzchar(barcode_norm) & !is.na(label_value) & nzchar(label_value)]
  if (nrow(ann) == 0) {
    stop("No non-empty label rows found in annotation metadata")
  }

  ann[, subtype_value := fifelse(is.na(subtype_value), "", subtype_value)]
  ann[, subtype_rank := fifelse(nzchar(subtype_value), 1L, 0L)]
  setorder(ann, barcode_norm, -subtype_rank, label_value, subtype_value)
  ann <- ann[, .SD[1], by = barcode_norm]
  ann[, subtype_rank := NULL]
  ann[]
}

apply_direct_labels <- function(object, annotation_path, gsm, gse, barcode_col, sample_col, label_col, subtype_col, project_root) {
  annotation <- fread(annotation_path)
  barcode_col <- choose_existing_column(
    annotation,
    barcode_col,
    c("cell_barcode", "barcode", "barcodes", "obs_name", "cell_id", "global_cell_id"),
    required = TRUE
  )
  sample_col <- choose_existing_column(
    annotation,
    sample_col,
    c("gsm", "source_gsm", "sample_gsm", "sample", "sample_id", "orig.ident", "pbmc_sample_id")
  )
  if (!label_col %in% colnames(annotation)) {
    stop("Label column not found in annotation metadata: ", label_col)
  }
  if (!is.null(subtype_col) && nzchar(subtype_col) && !subtype_col %in% colnames(annotation)) {
    subtype_col <- NULL
  }

  if (!is.null(sample_col)) {
    sample_aliases <- discover_sample_aliases(project_root, gse, gsm)
    sample_norm <- normalize_string(annotation[[sample_col]])
    keep <- sample_norm %in% sample_aliases
    if (any(keep)) {
      annotation <- annotation[keep, , drop = FALSE]
    }
  }

  collapsed <- collapse_annotation_rows(annotation, barcode_col, label_col, subtype_col)
  query_barcodes <- normalize_barcode(colnames(object))
  matched_idx <- match(query_barcodes, collapsed$barcode_norm)

  object$geo_celltype <- NA_character_
  matched <- !is.na(matched_idx)
  object$geo_celltype[matched] <- collapsed$label_value[matched_idx[matched]]

  if (!is.null(subtype_col)) {
    object$geo_celltype_subtype <- NA_character_
    object$geo_celltype_subtype[matched] <- collapsed$subtype_value[matched_idx[matched]]
    object$geo_celltype_subtype[!nzchar(object$geo_celltype_subtype)] <- NA_character_
  }

  list(
    object = object,
    matched_cells = sum(matched),
    total_cells = length(query_barcodes),
    barcode_column = barcode_col,
    sample_column = sample_col
  )
}

downsample_reference <- function(reference, group_col, max_cells_per_group, random_seed) {
  if (is.null(max_cells_per_group) || is.na(max_cells_per_group) || max_cells_per_group <= 0) {
    return(reference)
  }
  if (!group_col %in% colnames(reference@meta.data)) {
    return(reference)
  }

  grouping <- as.character(reference@meta.data[[group_col]])
  grouping[is.na(grouping) | !nzchar(grouping)] <- "Unlabeled"
  if (all(table(grouping) <= max_cells_per_group)) {
    return(reference)
  }

  set.seed(random_seed)
  cell_index <- seq_along(grouping)
  sampled_index <- unlist(
    lapply(split(cell_index, grouping), function(idx) {
      if (length(idx) <= max_cells_per_group) {
        idx
      } else {
        sample(idx, size = max_cells_per_group)
      }
    }),
    use.names = FALSE
  )
  sampled_index <- sort(sampled_index)
  subset(reference, cells = colnames(reference)[sampled_index])
}

run_rna_transfer <- function(
  object,
  reference_path,
  reference_assay,
  label_col,
  subtype_col,
  dims_use,
  max_cells_per_label,
  random_seed
) {
  reference <- readRDS(reference_path)
  if (!inherits(reference, "Seurat")) {
    stop("Reference RDS must contain a Seurat object")
  }
  if (!reference_assay %in% Assays(reference)) {
    stop("Reference assay not found: ", reference_assay)
  }
  if (!label_col %in% colnames(reference@meta.data)) {
    stop("Reference label column not found: ", label_col)
  }
  if (!is.null(subtype_col) && nzchar(subtype_col) && !subtype_col %in% colnames(reference@meta.data)) {
    subtype_col <- NULL
  }

  reference <- downsample_reference(reference, label_col, max_cells_per_label, random_seed)
  reference <- DietSeurat(reference, assays = reference_assay, dimreducs = NULL, graphs = NULL, misc = TRUE)
  DefaultAssay(reference) <- reference_assay
  if (identical(reference_assay, "RNA")) {
    reference <- NormalizeData(reference, verbose = FALSE)
  }
  if (length(VariableFeatures(reference)) == 0) {
    reference <- FindVariableFeatures(reference, nfeatures = min(5000, nrow(reference)), verbose = FALSE)
  }

  gene_activities <- GeneActivity(object)
  gene_activities <- gene_activities[intersect(rownames(gene_activities), rownames(reference[[reference_assay]])), , drop = FALSE]
  gene_activities <- gene_activities[Matrix::rowSums(gene_activities) > 0, , drop = FALSE]
  object[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
  DefaultAssay(object) <- "ACTIVITY"
  object <- NormalizeData(object, verbose = FALSE)

  transfer_features <- intersect(VariableFeatures(reference), rownames(object[["ACTIVITY"]]))
  if (length(transfer_features) < 200) {
    transfer_features <- intersect(rownames(reference[[reference_assay]]), rownames(object[["ACTIVITY"]]))
    transfer_features <- head(transfer_features, n = min(5000, length(transfer_features)))
  }
  if (length(transfer_features) < 50) {
    stop("Too few shared genes for RNA -> ATAC label transfer")
  }
  object <- ScaleData(object, features = transfer_features, verbose = FALSE)

  anchors <- FindTransferAnchors(
    reference = reference,
    query = object,
    reference.assay = reference_assay,
    query.assay = "ACTIVITY",
    reduction = "cca",
    features = transfer_features,
    dims = 1:30,
    verbose = FALSE
  )

  pred_main <- TransferData(
    anchorset = anchors,
    refdata = reference[[label_col]][, 1],
    weight.reduction = object[["lsi"]],
    dims = dims_use,
    verbose = FALSE
  )
  object$transferred_celltype <- pred_main$predicted.id
  main_score_cols <- grep("^prediction.score", colnames(pred_main), value = TRUE)
  object$transferred_celltype_score <- if (length(main_score_cols) > 0) {
    apply(pred_main[, main_score_cols, drop = FALSE], 1, max, na.rm = TRUE)
  } else {
    NA_real_
  }

  if (!is.null(subtype_col)) {
    pred_sub <- TransferData(
      anchorset = anchors,
      refdata = reference[[subtype_col]][, 1],
      weight.reduction = object[["lsi"]],
      dims = dims_use,
      verbose = FALSE
    )
    object$transferred_celltype_subtype <- pred_sub$predicted.id
    sub_score_cols <- grep("^prediction.score", colnames(pred_sub), value = TRUE)
    object$transferred_celltype_subtype_score <- if (length(sub_score_cols) > 0) {
      apply(pred_sub[, sub_score_cols, drop = FALSE], 1, max, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  DefaultAssay(object) <- "ATAC"
  list(object = object, features_used = length(transfer_features), reference_cells_used = ncol(reference))
}

save_metadata <- function(object, output_path) {
  meta <- as.data.table(object@meta.data)
  meta[, cell_barcode := colnames(object)]

  umap_embed <- Embeddings(object[["umap_atac"]])
  meta[, UMAP_ATAC_1 := umap_embed[, 1]]
  meta[, UMAP_ATAC_2 := umap_embed[, 2]]

  lsi_embed <- Embeddings(object[["lsi"]])
  max_dims <- min(10, ncol(lsi_embed))
  for (idx in seq_len(max_dims)) {
    meta[, (paste0("LSI_", idx)) := lsi_embed[, idx]]
  }

  write_csv_auto(meta, output_path)
}

option_list <- list(
  make_option("--gse", type = "character", default = NULL, dest = "gse"),
  make_option("--gsm", type = "character", default = NULL, dest = "gsm"),
  make_option("--sample-dir", type = "character", default = NULL, dest = "sample_dir"),
  make_option("--output-dir", type = "character", default = NULL, dest = "output_dir"),
  make_option("--annotation-csv", type = "character", default = NULL, dest = "annotation_csv"),
  make_option("--barcode-col", type = "character", default = NULL, dest = "barcode_col"),
  make_option("--sample-col", type = "character", default = NULL, dest = "sample_col"),
  make_option("--label-col", type = "character", default = "celltype", dest = "label_col"),
  make_option("--subtype-col", type = "character", default = "celltype_subtype", dest = "subtype_col"),
  make_option("--reference-rds", type = "character", default = NULL, dest = "reference_rds"),
  make_option("--reference-assay", type = "character", default = "RNA", dest = "reference_assay"),
  make_option("--reference-label-col", type = "character", default = "predicted.celltype.l2", dest = "reference_label_col"),
  make_option("--reference-subtype-col", type = "character", default = "predicted.celltype.l3", dest = "reference_subtype_col"),
  make_option("--reference-max-cells-per-label", type = "integer", default = 300L, dest = "reference_max_cells_per_label"),
  make_option("--reference-random-seed", type = "integer", default = 20260320L, dest = "reference_random_seed"),
  make_option("--npcs", type = "integer", default = 30L, dest = "npcs"),
  make_option("--dims-start", type = "integer", default = 2L, dest = "dims_start"),
  make_option("--dims-end", type = "integer", default = 30L, dest = "dims_end"),
  make_option("--resolution", type = "double", default = 0.5, dest = "resolution"),
  make_option("--save-rds", action = "store_true", default = FALSE, dest = "save_rds")
)

opt <- parse_args(OptionParser(option_list = option_list))
sample_dir <- resolve_sample_dir(project_root, opt$gse, opt$gsm, opt$sample_dir)

if (is.null(opt$gsm) || !nzchar(opt$gsm)) {
  opt$gsm <- basename(sample_dir)
}
if (is.null(opt$output_dir) || !nzchar(opt$output_dir)) {
  opt$output_dir <- sample_dir
}
output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- resolve_rds_path(sample_dir, opt$gsm)
query_obj <- readRDS(rds_path)
if (!inherits(query_obj, "Seurat")) {
  stop("Single-sample QC RDS must contain a Seurat object")
}
if (!"pass_qc" %in% colnames(query_obj@meta.data)) {
  stop("Single-sample QC object is missing pass_qc metadata")
}

query_obj <- subset(query_obj, subset = pass_qc %in% TRUE)
if (ncol(query_obj) < 50) {
  stop("Too few QC-passed cells for single-sample UMAP: ", ncol(query_obj))
}
query_obj <- repair_fragment_path(query_obj, project_root, opt$gse, opt$gsm)

dims_requested <- parse_int_list(opt$dims_start, opt$dims_end)
prepared <- prepare_query_object(query_obj, opt$npcs, dims_requested, opt$resolution)
query_obj <- prepared$object
dims_use <- prepared$dims_use

qc_features <- c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "FRiP", "nucleosome_signal", "blacklist_fraction", "scDblFinder.score")
qc_features <- qc_features[qc_features %in% colnames(query_obj@meta.data)]

cluster_plot_path <- file.path(output_dir, "single_sample_umap_by_cluster.png")
qc_plot_path <- file.path(output_dir, "single_sample_umap_by_qc.png")
metadata_path <- file.path(output_dir, "single_sample_umap_metadata.csv.gz")
report_path <- file.path(output_dir, "single_sample_umap_report.md")

save_plot(
  cluster_plot_path,
  make_dim_plot(query_obj, "seurat_clusters", paste0(opt$gsm, " single-sample ATAC UMAP by cluster")),
  width = 10,
  height = 8
)
if (length(qc_features) > 0) {
  save_plot(
    qc_plot_path,
    make_feature_plot(query_obj, qc_features, paste0(opt$gsm, " single-sample ATAC UMAP by QC metrics")),
    width = 12,
    height = 10
  )
}

report_lines <- c(
  "# Single-sample UMAP report",
  "",
  paste0("- Sample directory: `", sample_dir, "`"),
  paste0("- QC RDS: `", rds_path, "`"),
  paste0("- Cells after QC: ", format(ncol(query_obj), big.mark = ",")),
  paste0("- Peaks: ", format(nrow(query_obj), big.mark = ",")),
  paste0("- Clusters: ", length(unique(query_obj$seurat_clusters))),
  paste0("- LSI dims used: `", format_dims(dims_use), "`")
)

if (!is.null(opt$annotation_csv) && nzchar(opt$annotation_csv)) {
  annotation_path <- normalizePath(opt$annotation_csv, winslash = "/", mustWork = TRUE)
  gse_value <- if (is.null(opt$gse) || !nzchar(opt$gse)) "" else opt$gse
  direct_result <- apply_direct_labels(
    query_obj,
    annotation_path,
    opt$gsm,
    gse_value,
    opt$barcode_col,
    opt$sample_col,
    opt$label_col,
    opt$subtype_col,
    project_root
  )
  query_obj <- direct_result$object

  geo_plot_path <- file.path(output_dir, "single_sample_umap_by_geo_celltype.png")
  save_plot(
    geo_plot_path,
    make_dim_plot(query_obj, "geo_celltype", paste0(opt$gsm, " single-sample UMAP by GEO celltype")),
    width = 11,
    height = 8
  )

  report_lines <- c(
    report_lines,
    paste0("- GEO annotation metadata: `", annotation_path, "`"),
    paste0("- GEO label mapping matched cells: ", format(direct_result$matched_cells, big.mark = ","), " / ", format(direct_result$total_cells, big.mark = ",")),
    paste0("- GEO barcode column: `", direct_result$barcode_column, "`")
  )
  if (!is.null(direct_result$sample_column)) {
    report_lines <- c(report_lines, paste0("- GEO sample column: `", direct_result$sample_column, "`"))
  }

  if ("geo_celltype_subtype" %in% colnames(query_obj@meta.data) && any(!is.na(query_obj$geo_celltype_subtype))) {
    geo_subtype_plot_path <- file.path(output_dir, "single_sample_umap_by_geo_celltype_subtype.png")
    save_plot(
      geo_subtype_plot_path,
      make_dim_plot(query_obj, "geo_celltype_subtype", paste0(opt$gsm, " single-sample UMAP by GEO subtype")),
      width = 13,
      height = 10
    )
  }
}

if (!is.null(opt$reference_rds) && nzchar(opt$reference_rds)) {
  reference_path <- normalizePath(opt$reference_rds, winslash = "/", mustWork = TRUE)
  transfer_result <- run_rna_transfer(
    query_obj,
    reference_path,
    opt$reference_assay,
    opt$reference_label_col,
    opt$reference_subtype_col,
    dims_use,
    opt$reference_max_cells_per_label,
    opt$reference_random_seed
  )
  query_obj <- transfer_result$object

  transfer_plot_path <- file.path(output_dir, "single_sample_umap_by_transferred_celltype.png")
  save_plot(
    transfer_plot_path,
    make_dim_plot(query_obj, "transferred_celltype", paste0(opt$gsm, " single-sample UMAP by transferred celltype")),
    width = 11,
    height = 8
  )

  report_lines <- c(
    report_lines,
    paste0("- RNA reference RDS: `", reference_path, "`"),
      paste0("- RNA reference assay: `", opt$reference_assay, "`"),
      paste0("- RNA transfer label column: `", opt$reference_label_col, "`"),
      paste0("- RNA reference sampled cells: ", format(transfer_result$reference_cells_used, big.mark = ",")),
      paste0("- RNA transfer shared genes used: ", transfer_result$features_used)
    )

  if ("transferred_celltype_score" %in% colnames(query_obj@meta.data)) {
    report_lines <- c(
      report_lines,
      paste0(
        "- Median transferred-celltype score: ",
        sprintf("%.4f", median(query_obj$transferred_celltype_score, na.rm = TRUE))
      )
    )
  }

  if ("transferred_celltype_subtype" %in% colnames(query_obj@meta.data) && any(!is.na(query_obj$transferred_celltype_subtype))) {
    transfer_subtype_plot_path <- file.path(output_dir, "single_sample_umap_by_transferred_celltype_subtype.png")
    save_plot(
      transfer_subtype_plot_path,
      make_dim_plot(query_obj, "transferred_celltype_subtype", paste0(opt$gsm, " single-sample UMAP by transferred subtype")),
      width = 13,
      height = 10
    )
  }
}

save_metadata(query_obj, metadata_path)

if (isTRUE(opt$save_rds)) {
  saveRDS(query_obj, file.path(output_dir, paste0(opt$gsm, "_single_sample_umap.rds")))
}

output_files <- c(
  "- `single_sample_umap_by_cluster.png`",
  if (length(qc_features) > 0) "- `single_sample_umap_by_qc.png`",
  if ("geo_celltype" %in% colnames(query_obj@meta.data)) "- `single_sample_umap_by_geo_celltype.png`",
  if ("geo_celltype_subtype" %in% colnames(query_obj@meta.data) && any(!is.na(query_obj$geo_celltype_subtype))) "- `single_sample_umap_by_geo_celltype_subtype.png`",
  if ("transferred_celltype" %in% colnames(query_obj@meta.data)) "- `single_sample_umap_by_transferred_celltype.png`",
  if ("transferred_celltype_subtype" %in% colnames(query_obj@meta.data) && any(!is.na(query_obj$transferred_celltype_subtype))) "- `single_sample_umap_by_transferred_celltype_subtype.png`",
  "- `single_sample_umap_metadata.csv.gz`",
  "- `single_sample_umap_report.md`"
)
output_files <- output_files[!is.na(output_files)]

report_lines <- c(
  report_lines,
  "",
  "## Output files",
  "",
  output_files
)
report_lines <- report_lines[!is.na(report_lines)]
writeLines(report_lines, con = report_path)

cat("Single-sample UMAP outputs:", output_dir, "\n")
