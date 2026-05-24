#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(data.table)
})

option_list <- list(
  make_option(c("--fragment-file"), type = "character", help = "Input fragments.tsv.gz file"),
  make_option(c("--sample-id"), type = "character", help = "Sample identifier"),
  make_option(c("--output-dir"), type = "character", help = "Output directory"),
  make_option(c("--min-fragments"), type = "integer", default = 200L, help = "Minimum fragments per barcode [default %default]"),
  make_option(c("--max-barcodes"), type = "integer", default = 50000L, help = "Maximum selected barcodes after ranking [default %default]"),
  make_option(c("--min-tss"), type = "double", default = 0, help = "Minimum ArchR TSSEnrichment when available [default %default]"),
  make_option(c("--rank-by"), type = "character", default = "fragments", help = "Barcode ranking method: fragments or tss_then_fragments [default %default]"),
  make_option(c("--force"), action = "store_true", default = FALSE, help = "Overwrite existing outputs")
)

opt <- parse_args(OptionParser(option_list = option_list))

get_opt <- function(name) {
  opt[[name]]
}

fragment_file <- get_opt("fragment-file")
sample_id <- get_opt("sample-id")
output_dir_arg <- get_opt("output-dir")
min_fragments <- get_opt("min-fragments")
max_barcodes <- get_opt("max-barcodes")
min_tss <- get_opt("min-tss")
rank_by <- get_opt("rank-by")
if (!rank_by %in% c("fragments", "tss_then_fragments")) {
  stop("Unsupported --rank-by: ", rank_by, call. = FALSE)
}

required <- list(
  "fragment-file" = fragment_file,
  "sample-id" = sample_id,
  "output-dir" = output_dir_arg
)
for (field in names(required)) {
  if (is.null(required[[field]]) || !nzchar(required[[field]])) {
    stop("Missing required option --", field, call. = FALSE)
  }
}

if (!file.exists(fragment_file)) {
  stop("Fragment file not found: ", fragment_file, call. = FALSE)
}

if (!requireNamespace("ArchR", quietly = TRUE)) {
  stop(
    "ArchR is required for longevity ATAC barcode preprocessing but is not installed. ",
    "Install ArchR in the R environment, then rerun this command.",
    call. = FALSE
  )
}

output_dir <- normalizePath(output_dir_arg, winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

barcode_path <- file.path(output_dir, "filtered_barcodes.tsv.gz")
qc_path <- file.path(output_dir, "barcode_qc.csv.gz")
summary_path <- file.path(output_dir, "summary.json")
archr_dir <- file.path(output_dir, "archr")
arrow_path <- file.path(archr_dir, paste0(sample_id, ".arrow"))

existing <- c(barcode_path, qc_path, summary_path)
if (!opt$force && any(file.exists(existing))) {
  stop("Barcode preprocessing outputs already exist. Use --force to overwrite: ", output_dir, call. = FALSE)
}

dir.create(archr_dir, recursive = TRUE, showWarnings = FALSE)

message("[archr-barcode-preprocess] sample=", sample_id)
message("fragment_file=", fragment_file)
message("output_dir=", output_dir)

suppressPackageStartupMessages({
  library(ArchR)
})

ArchR::addArchRGenome("hg38")
ArchR::addArchRThreads(threads = max(1L, parallel::detectCores(logical = TRUE) - 2L))

existing_arrow <- c(arrow_path, paste0(arrow_path, ".arrow"))
existing_arrow <- existing_arrow[file.exists(existing_arrow)]
if (length(existing_arrow) > 0) {
  arrow_files <- existing_arrow[1]
} else {
  arrow_files <- ArchR::createArrowFiles(
    inputFiles = setNames(fragment_file, sample_id),
    sampleNames = sample_id,
    outputNames = arrow_path,
    minTSS = min_tss,
    minFrags = min_fragments,
    addTileMat = TRUE,
    addGeneScoreMat = FALSE,
    force = TRUE
  )
}

project <- ArchR::ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = archr_dir,
  copyArrows = FALSE
)

cell_qc <- as.data.frame(ArchR::getCellColData(project))
cell_qc$cell_name <- rownames(cell_qc)

fragment_col <- intersect(c("nFrags", "nfrags", "fragments"), colnames(cell_qc))
if (length(fragment_col) == 0) {
  stop("ArchR cell QC did not contain an nFrags column", call. = FALSE)
}
fragment_col <- fragment_col[1]

tss_col <- intersect(c("TSSEnrichment", "tss_enrichment", "TSS"), colnames(cell_qc))
tss_col <- if (length(tss_col) > 0) tss_col[1] else NA_character_

barcode <- sub("^.*#", "", cell_qc$cell_name)
cell_qc$barcode <- barcode
cell_qc$n_fragments_for_filter <- suppressWarnings(as.numeric(cell_qc[[fragment_col]]))
if (!is.na(tss_col)) {
  cell_qc$tss_for_filter <- suppressWarnings(as.numeric(cell_qc[[tss_col]]))
} else {
  cell_qc$tss_for_filter <- NA_real_
}

selected <- !is.na(cell_qc$n_fragments_for_filter) & cell_qc$n_fragments_for_filter >= min_fragments
if (!is.na(tss_col)) {
  selected <- selected & !is.na(cell_qc$tss_for_filter) & cell_qc$tss_for_filter >= min_tss
}

ranked <- cell_qc[selected & nzchar(cell_qc$barcode), , drop = FALSE]
if (rank_by == "tss_then_fragments" && !is.na(tss_col)) {
  ranked <- ranked[order(ranked$tss_for_filter, ranked$n_fragments_for_filter, decreasing = TRUE), , drop = FALSE]
} else {
  ranked <- ranked[order(ranked$n_fragments_for_filter, decreasing = TRUE), , drop = FALSE]
}
if (!is.na(max_barcodes) && is.finite(max_barcodes) && max_barcodes > 0 && nrow(ranked) > max_barcodes) {
  ranked <- ranked[seq_len(max_barcodes), , drop = FALSE]
}

cell_qc$selected <- cell_qc$barcode %in% ranked$barcode

qc_con <- gzfile(qc_path, open = "wt")
utils::write.csv(cell_qc, qc_con, row.names = FALSE, quote = TRUE)
close(qc_con)
con <- gzfile(barcode_path, open = "wt")
writeLines(ranked$barcode, con = con, sep = "\n")
close(con)

summary <- list(
  sample_id = sample_id,
  fragment_file = normalizePath(fragment_file, winslash = "/", mustWork = TRUE),
  output_dir = output_dir,
  method = "ArchR_createArrowFiles_plus_project_filter",
  min_fragments = min_fragments,
  max_barcodes = max_barcodes,
  min_tss = min_tss,
  rank_by = rank_by,
  archr_version = as.character(utils::packageVersion("ArchR")),
  n_archr_cells = nrow(cell_qc),
  n_selected_barcodes = length(ranked$barcode),
  fragment_column = fragment_col,
  tss_column = ifelse(is.na(tss_col), "", tss_col),
  barcode_file = barcode_path,
  barcode_qc_file = qc_path,
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), summary_path)

message("selected_barcodes=", length(ranked$barcode))
message("barcode_file=", barcode_path)
