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

celltype_map <- data.frame(
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

save_umap_plot <- function(df, color_by, out_path, title_text, width = 10, height = 8) {
  palette <- category_palette(df[[color_by]])
  p <- ggplot(df, aes(x = UMAP_Harmony_1, y = UMAP_Harmony_2, color = .data[[color_by]])) +
    geom_point(size = 0.08, alpha = 0.8, stroke = 0) +
    scale_color_manual(values = palette) +
    labs(title = title_text, x = "UMAPHARM_1", y = "UMAPHARM_2", color = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.35, "cm")
    )
  ggsave(out_path, plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

option_list <- list(
  make_option("--analysis-dir", type = "character", dest = "analysis_dir")
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
annotated <- merge(meta, celltype_map, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
if (any(is.na(annotated$celltype))) {
  missing_clusters <- sort(unique(annotated$seurat_clusters[is.na(annotated$celltype)]))
  stop("Unmapped clusters detected: ", paste(missing_clusters, collapse = ", "))
}

fwrite(celltype_map, file.path(analysis_dir, "cluster_celltype_annotation_map.csv"))
fwrite(annotated, file.path(analysis_dir, "integrated_metadata_celltyped.csv.gz"))

save_umap_plot(
  annotated,
  "celltype",
  file.path(analysis_dir, "post_harmony_umap_by_celltype.png"),
  "Post-Harmony UMAP by broad cell type",
  width = 10,
  height = 8
)
save_umap_plot(
  annotated,
  "celltype_subtype",
  file.path(analysis_dir, "post_harmony_umap_by_celltype_subtype.png"),
  "Post-Harmony UMAP by cell type subtype",
  width = 14,
  height = 10
)

summary_lines <- c(
  "# Celltype annotation notes",
  "",
  "- Source datasets `GSE190992` and `GSE283744` are both PBMC-based, so unexpected non-immune marker clusters are conservatively labeled `Unknown / non-PBMC-like`.",
  "- Broad labels are intended for visualization only; subtype labels are marker-based shorthand rather than definitive manual curation.",
  "- Evidence markers are taken from `cluster_top_accessible_peaks.csv` nearest genes."
)
writeLines(summary_lines, con = file.path(analysis_dir, "celltype_annotation_notes.md"))

cat("Celltype annotation outputs:", analysis_dir, "\n")
