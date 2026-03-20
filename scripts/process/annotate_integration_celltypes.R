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
  library(data.table)
  library(ggplot2)
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

default_celltype_map <- data.frame(
  seurat_clusters = as.character(0:21),
  celltype = c(
    "Unknown / non-PBMC-like",
    "Unknown / non-PBMC-like",
    "B cell lineage",
    "T / NK",
    "Myeloid / innate-like",
    "T / NK",
    "T / NK",
    "Unknown / non-PBMC-like",
    "Unknown / non-PBMC-like",
    "Unknown / non-PBMC-like",
    "Myeloid / innate-like",
    "Unknown / non-PBMC-like",
    "T / NK",
    "Myeloid / innate-like",
    "T / NK",
    "Myeloid / innate-like",
    "Megakaryocyte / platelet",
    "T / NK",
    "T / NK",
    "B cell lineage",
    "Dendritic cell",
    "Unknown / non-PBMC-like"
  ),
  celltype_subtype = c(
    "Unknown_EDAR_ISM1",
    "Unknown_ITGB4_NECTIN4",
    "B_CD79A_BLK",
    "CD8_T_CRTAM",
    "Innate_APOBEC3A_RORC",
    "Cytotoxic_PRF1_KLF2",
    "T_NK_BCL11B_ADGRG1",
    "Unknown_DAB1_ANXA6",
    "Unknown_ITPR3_BIRC3",
    "Unknown_KRT74_AKNA",
    "Monocyte_like_CALD1_PRG4",
    "Unknown_LPCAT1_PTPRN2",
    "Activated_T_NK_CD8A_TNFSF9",
    "Myeloid_like_PRG4_PFKFB3",
    "Activated_T_NK_IL2RB_TNFRSF4",
    "Myeloid_like_SULF2_ADTRP",
    "Platelet_GP5",
    "T_NK_BCL11B_ADGRG1_2",
    "Cytotoxic_UNC13D",
    "B_IGLL5_CECR3",
    "cDC1_like_XCR1",
    "Unknown_PIK3AP1_GAD1"
  ),
  evidence_markers = c(
    "EDAR; ISM1; PROX1-AS1",
    "ITGB4; NECTIN4; TIFAB",
    "CD79A; BLK; PRDM4",
    "CD8A; CD8B; CRTAM",
    "APOBEC3A; RORC; FAM129A",
    "PRF1; KLF2; CCL18",
    "BCL11B; ADGRG1; F2R",
    "DAB1; ANXA6; ADARB1",
    "ITPR3; BIRC3; ARHGEF4",
    "KRT74; AKNA; LZTS1",
    "CALD1; PRG4; ISM1",
    "LPCAT1; PTPRN2; CDH23",
    "CD8A; TNFSF9; PRKCQ",
    "PRG4; PFKFB3; ZNF787",
    "IL2RB; TNFRSF4; IFNL1",
    "SULF2; ADTRP; METTL22",
    "GP5; CECR3; MYBL2",
    "BCL11B; ADGRG1; MARCH4",
    "UNC13D; MAPK14; NEDD9",
    "IGLL5; CECR3; TNFRSF1B",
    "XCR1; BLK; SPON1",
    "PIK3AP1; GAD1; ARAP1"
  ),
  annotation_confidence = c(
    "low", "low", "high", "high", "medium", "medium", "medium", "low", "low", "low", "medium",
    "low", "high", "medium", "high", "medium", "high", "medium", "high", "high", "medium", "low"
  ),
  stringsAsFactors = FALSE
)

normalize_cluster_id <- function(values) {
  sub("^g", "", as.character(values))
}

build_marker_lookup <- function(analysis_dir) {
  marker_path <- file.path(analysis_dir, "cluster_top_accessible_peaks.csv")
  if (!file.exists(marker_path)) {
    return(data.table(seurat_clusters = character(), evidence_markers = character()))
  }

  markers <- fread(marker_path)
  if (!all(c("cluster", "gene_name") %in% colnames(markers))) {
    return(data.table(seurat_clusters = character(), evidence_markers = character()))
  }

  markers[, seurat_clusters := normalize_cluster_id(cluster)]
  markers <- markers[!is.na(gene_name) & nzchar(gene_name)]
  marker_summary <- markers[, .(
    evidence_markers = paste(unique(gene_name)[seq_len(min(3L, uniqueN(gene_name)))], collapse = "; ")
  ), by = seurat_clusters]
  marker_summary[]
}

derive_confidence <- function(celltype_fraction, subtype_fraction) {
  if (is.na(celltype_fraction) || is.na(subtype_fraction)) {
    return("low")
  }
  if (celltype_fraction >= 0.9 && subtype_fraction >= 0.75) {
    return("high")
  }
  if (celltype_fraction >= 0.7 && subtype_fraction >= 0.5) {
    return("medium")
  }
  "low"
}

derive_celltype_map_from_reference <- function(meta, reference_meta, analysis_dir) {
  required_reference_cols <- c("global_cell_id", "celltype", "celltype_subtype")
  missing_reference_cols <- required_reference_cols[!required_reference_cols %in% colnames(reference_meta)]
  if (length(missing_reference_cols) > 0) {
    stop("Missing required reference metadata columns: ", paste(missing_reference_cols, collapse = ", "))
  }
  if (!"global_cell_id" %in% colnames(meta)) {
    stop("Current integrated metadata must contain global_cell_id for reference-driven annotation")
  }

  marker_lookup <- build_marker_lookup(analysis_dir)
  reference_subset <- as.data.table(reference_meta)[, .(
    global_cell_id,
    reference_celltype = celltype,
    reference_subtype = celltype_subtype
  )]

  overlap <- merge(
    as.data.table(meta)[, .(global_cell_id, seurat_clusters)],
    reference_subset,
    by = "global_cell_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(overlap$reference_celltype))) {
    missing_n <- sum(is.na(overlap$reference_celltype))
    stop("Reference metadata missing annotations for ", missing_n, " cells in current analysis")
  }

  overlap_summary <- overlap[, .N, by = .(seurat_clusters, reference_celltype, reference_subtype)]
  overlap_summary[, cluster_order := as.integer(seurat_clusters)]
  setorder(overlap_summary, cluster_order, -N, reference_celltype, reference_subtype)
  overlap_summary[, cluster_order := NULL]

  celltype_winners <- overlap[, .N, by = .(seurat_clusters, reference_celltype)]
  celltype_winners[, cluster_order := as.integer(seurat_clusters)]
  setorder(celltype_winners, cluster_order, -N, reference_celltype)
  celltype_winners[, cluster_order := NULL]
  celltype_winners[, total_cells := sum(N), by = seurat_clusters]
  celltype_winners <- celltype_winners[, .SD[1], by = seurat_clusters]
  setnames(celltype_winners, c("reference_celltype", "N"), c("celltype", "celltype_votes"))
  celltype_winners[, celltype_fraction := celltype_votes / total_cells]

  subtype_winners <- overlap[, .N, by = .(seurat_clusters, reference_subtype)]
  subtype_winners[, cluster_order := as.integer(seurat_clusters)]
  setorder(subtype_winners, cluster_order, -N, reference_subtype)
  subtype_winners[, cluster_order := NULL]
  subtype_winners[, total_cells := sum(N), by = seurat_clusters]
  subtype_winners <- subtype_winners[, .SD[1], by = seurat_clusters]
  setnames(subtype_winners, c("reference_subtype", "N"), c("celltype_subtype", "subtype_votes"))
  subtype_winners[, subtype_fraction := subtype_votes / total_cells]

  celltype_map <- merge(
    celltype_winners[, .(seurat_clusters, celltype, total_cells, celltype_votes, celltype_fraction)],
    subtype_winners[, .(seurat_clusters, celltype_subtype, subtype_votes, subtype_fraction)],
    by = "seurat_clusters",
    all = TRUE,
    sort = FALSE
  )
  celltype_map <- merge(celltype_map, marker_lookup, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
  celltype_map[is.na(evidence_markers) | !nzchar(evidence_markers), evidence_markers := "reference-transfer"]
  celltype_map[, annotation_confidence := mapply(derive_confidence, celltype_fraction, subtype_fraction)]
  celltype_map[, cluster_order := as.integer(seurat_clusters)]
  setorder(celltype_map, cluster_order)
  celltype_map[, cluster_order := NULL]

  fwrite(overlap_summary, file.path(analysis_dir, "cluster_annotation_reference_overlap.csv"))
  celltype_map[]
}

category_palette <- function(values) {
  values <- unique(as.character(values))
  base_colors <- c(
    "#0f766e", "#b45309", "#2563eb", "#be123c", "#7c3aed", "#0891b2", "#65a30d", "#c2410c",
    "#4f46e5", "#15803d", "#dc2626", "#1d4ed8", "#a21caf", "#0369a1", "#ca8a04", "#0f766e",
    "#6d28d9", "#0e7490", "#3f6212", "#9a3412", "#1e40af", "#9f1239", "#5b21b6", "#166534"
  )
  cols <- rep(base_colors, length.out = length(values))
  names(cols) <- values
  cols
}

fixed_celltype_palette <- c(
  "B cell lineage" = "#be123c",
  "T / NK" = "#b45309",
  "Myeloid / innate-like" = "#0f766e",
  "Dendritic cell" = "#7c3aed",
  "Megakaryocyte / platelet" = "#db2777",
  "Unknown / unresolved" = "#2563eb",
  "Unknown / non-PBMC-like" = "#2563eb"
)

resolve_plot_palette <- function(values, color_by) {
  values <- unique(as.character(values))
  if (!identical(color_by, "celltype")) {
    return(list(levels = values, colors = category_palette(values)))
  }

  ordered_values <- c(
    names(fixed_celltype_palette)[names(fixed_celltype_palette) %in% values],
    sort(setdiff(values, names(fixed_celltype_palette)))
  )
  ordered_values <- unique(ordered_values)

  palette <- fixed_celltype_palette[ordered_values[ordered_values %in% names(fixed_celltype_palette)]]
  extra_values <- setdiff(ordered_values, names(fixed_celltype_palette))
  if (length(extra_values) > 0) {
    extra_palette <- category_palette(extra_values)
    palette <- c(palette, extra_palette[extra_values])
  }

  list(levels = ordered_values, colors = palette)
}

save_umap_plot <- function(df, color_by, out_path, title_text, width = 10, height = 8, point_size = 0.14) {
  palette_spec <- resolve_plot_palette(df[[color_by]], color_by)
  if (identical(color_by, "celltype")) {
    df[[color_by]] <- factor(as.character(df[[color_by]]), levels = palette_spec$levels)
  }

  p <- ggplot(df, aes(x = UMAP_Harmony_1, y = UMAP_Harmony_2, color = .data[[color_by]])) +
    geom_point(size = point_size, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = palette_spec$colors, breaks = palette_spec$levels, drop = TRUE) +
    labs(title = title_text, x = "UMAPHARM_1", y = "UMAPHARM_2", color = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.45, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.8, alpha = 1)))
  ggsave(out_path, plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

option_list <- list(
  make_option("--analysis-dir", type = "character", dest = "analysis_dir"),
  make_option("--reference-metadata", type = "character", default = NULL, dest = "reference_metadata"),
  make_option("--annotation-map", type = "character", default = NULL, dest = "annotation_map"),
  make_option("--output-suffix", type = "character", default = "", dest = "output_suffix")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$analysis_dir)) {
  stop("--analysis-dir is required")
}

analysis_dir <- normalizePath(opt$analysis_dir, winslash = "/", mustWork = TRUE)
metadata_path <- file.path(analysis_dir, "integrated_metadata.csv.gz")
if (!file.exists(metadata_path)) {
  stop("Missing integrated_metadata.csv.gz under ", analysis_dir)
}

meta <- fread(metadata_path)
required_cols <- c("seurat_clusters", "UMAP_Harmony_1", "UMAP_Harmony_2")
missing_cols <- required_cols[!required_cols %in% colnames(meta)]
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta[, seurat_clusters := as.character(seurat_clusters)]
if (!is.null(opt$annotation_map)) {
  annotation_map_path <- normalizePath(opt$annotation_map, winslash = "/", mustWork = TRUE)
  celltype_map <- fread(annotation_map_path)
} else if (is.null(opt$reference_metadata)) {
  celltype_map <- copy(default_celltype_map)
} else {
  reference_metadata_path <- normalizePath(opt$reference_metadata, winslash = "/", mustWork = TRUE)
  reference_meta <- fread(reference_metadata_path)
  reference_meta[, seurat_clusters := NULL]
  celltype_map <- derive_celltype_map_from_reference(meta, reference_meta, analysis_dir)
}

required_map_cols <- c("seurat_clusters", "celltype", "celltype_subtype", "evidence_markers", "annotation_confidence")
missing_map_cols <- required_map_cols[!required_map_cols %in% colnames(celltype_map)]
if (length(missing_map_cols) > 0) {
  stop("Missing required annotation-map columns: ", paste(missing_map_cols, collapse = ", "))
}
celltype_map[, seurat_clusters := as.character(seurat_clusters)]

output_suffix <- opt$output_suffix
if (!nzchar(output_suffix)) {
  output_suffix <- ""
}

map_out <- file.path(analysis_dir, paste0("cluster_celltype_annotation_map", output_suffix, ".csv"))
annotated_out <- file.path(analysis_dir, paste0("integrated_metadata_celltyped", output_suffix, ".csv.gz"))
umap_celltype_out <- file.path(analysis_dir, paste0("post_harmony_umap_by_celltype", output_suffix, ".png"))
umap_subtype_out <- file.path(analysis_dir, paste0("post_harmony_umap_by_celltype_subtype", output_suffix, ".png"))
notes_out <- file.path(analysis_dir, paste0("celltype_annotation_notes", output_suffix, ".md"))

annotated <- merge(meta, celltype_map, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
if (any(is.na(annotated$celltype))) {
  missing_clusters <- sort(unique(annotated$seurat_clusters[is.na(annotated$celltype)]))
  stop("Unmapped clusters detected: ", paste(missing_clusters, collapse = ", "))
}

write_csv_auto(celltype_map, map_out)
write_csv_auto(annotated, annotated_out)

save_umap_plot(
  annotated,
  "celltype",
  umap_celltype_out,
  "Post-Harmony UMAP by broad cell type",
  width = 10,
  height = 8
)
save_umap_plot(
  annotated,
  "celltype_subtype",
  umap_subtype_out,
  "Post-Harmony UMAP by cell type subtype",
  width = 14,
  height = 10,
  point_size = 0.18
)

summary_lines <- c(
  "# Celltype annotation notes",
  "",
  if (!is.null(opt$annotation_map)) {
    paste0("- Cluster labels are read from external annotation map `", normalizePath(opt$annotation_map, winslash = "/", mustWork = TRUE), "`.")
  } else if (is.null(opt$reference_metadata)) {
    "- Source datasets `GSE190992` and `GSE283744` are both PBMC-based, so unexpected non-immune marker clusters are conservatively labeled `Unknown / non-PBMC-like`."
  } else {
    paste0("- Cluster labels are transferred from reference metadata `", normalizePath(opt$reference_metadata, winslash = "/", mustWork = TRUE), "` using shared `global_cell_id` and majority vote within each new cluster.")
  },
  "- Broad labels are intended for visualization only; subtype labels are marker-based shorthand rather than definitive manual curation.",
  "- Evidence markers are taken from the current analysis `cluster_top_accessible_peaks.csv` nearest genes when available.",
  if (!is.null(opt$reference_metadata)) {
    "- `cluster_annotation_reference_overlap.csv` records the old-annotation composition of each new cluster for traceability."
  }
)
summary_lines <- summary_lines[!is.na(summary_lines)]
writeLines(summary_lines, con = notes_out)

cat("Celltype annotation outputs:", analysis_dir, "\n")
