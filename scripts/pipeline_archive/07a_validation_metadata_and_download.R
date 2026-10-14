# CRC-ImmuneMap | Step 07A
# GSE144735 (KUL3 cohort): retrieve processed annotation/raw UMI, inspect suitability
# This is a PREPARATION / FEASIBILITY step, not validation scoring or a hypothesis test.
# No renormalization, subclustering, or inference is performed here.

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
if (!dir.exists(project_dir)) stop("Project directory missing. Edit project_dir at top of script.")
setwd(project_dir)

options(timeout = max(getOption("timeout"), 1200))
raw_dir <- file.path("data", "raw", "GSE144735")
out_dir <- file.path("results", "tables", "validation")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144735/suppl/"
files <- c(
  annotation = "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz",
  raw_counts = "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
)

# First prepare the small annotation and raw UMI (~54 MB compressed).
# The natural-log TPM matrix is unnecessary for the next stage.
get_gse_file <- function(filename) {
  dest <- file.path(raw_dir, filename)
  if (file.exists(dest) && file.size(dest) > 0L) {
    cat("Already downloaded:", dest, "\n")
    return(dest)
  }
  tmp <- paste0(dest, ".partial")
  if (file.exists(tmp)) unlink(tmp)
  url <- paste0(base_url, filename)
  cat("Downloading:", url, "\n")
  tryCatch(
    download.file(url, destfile = tmp, mode = "wb", method = "libcurl", quiet = FALSE),
    error = function(e) stop(
      "Download failed: ", conditionMessage(e), "\n",
      "Try downloading from the GSE144735 GEO supplementary-files page, ",
      "place the .gz file in: ", normalizePath(raw_dir, winslash = "/", mustWork = FALSE),
      " then rerun."
    )
  )
  if (!file.exists(tmp) || file.size(tmp) <= 0L) stop("Download missing/empty: ", tmp)
  if (!file.rename(tmp, dest)) stop("Unable to rename downloaded file: ", tmp)
  dest
}

paths <- vapply(files, get_gse_file, character(1))

# Basic readable-gzip check and text preview; does not load the count matrix in memory.
preview_gz <- function(path, n = 2L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  lines <- readLines(con, n = n, warn = FALSE)
  if (length(lines) < n) stop("Unexpected/truncated content in: ", path)
  lines
}

annotation_preview <- preview_gz(paths[["annotation"]])
count_preview <- preview_gz(paths[["raw_counts"]])
cat("\n================ GSE144735 FILE CHECK ================\n")
cat("Annotation file bytes:", file.size(paths[["annotation"]]), "\n")
cat("Count file bytes:", file.size(paths[["raw_counts"]]), "\n")
cat("Count matrix header entries:",
    length(strsplit(count_preview[1], "\t", fixed = TRUE)[[1]]), "\n")
cat("Count matrix first data-row entries:",
    length(strsplit(count_preview[2], "\t", fixed = TRUE)[[1]]), "\n")
cat("First 6 count header fields:\n")
print(head(strsplit(count_preview[1], "\t", fixed = TRUE)[[1]], 6))
cat("First 6 count first-row fields:\n")
print(head(strsplit(count_preview[2], "\t", fixed = TRUE)[[1]], 6))

# Read author-provided metadata (small file only).
annotation <- read.delim(
  paths[["annotation"]], sep = "\t", header = TRUE,
  check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""
)
if (nrow(annotation) < 1L) stop("Empty annotation; inspect delimiter and format.")
cat("\n================ ANNOTATION CHECK ================\n")
cat("Annotation dimensions:", nrow(annotation), "rows x", ncol(annotation), "columns\n")
print(colnames(annotation))
cat("First 4 metadata rows:\n")
print(utils::head(annotation, 4), row.names = FALSE)

# Save schema; never infer if supplied column names differ.
schema <- data.frame(
  column_name = colnames(annotation),
  R_class = vapply(annotation, function(x) paste(class(x), collapse = "/"), character(1)),
  n_missing = vapply(annotation, function(x) sum(is.na(x) | trimws(as.character(x)) == ""), integer(1)),
  n_unique = vapply(annotation, function(x) length(unique(x)), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(schema, file.path(out_dir, "validation_annotation_schema.csv"), row.names = FALSE)

pick_col <- function(nms, possibilities) {
  first <- possibilities[possibilities %in% nms]
  if (length(first)) first[1] else NA_character_
}
col_cell <- pick_col(names(annotation), c("Index", "index", "Cell", "cell", "cell_id", "Cell_ID", "Barcode", "barcode"))
col_sample <- pick_col(names(annotation), c("Sample", "sample", "sample_id", "SampleID"))
col_class <- pick_col(names(annotation), c("Class", "class", "condition", "Condition", "Tissue", "tissue"))
col_type <- pick_col(names(annotation), c("Cell_type", "CellType", "Cell_type_1", "Cell_Type", "cell_type", "Celltype", "CellType1", "celltype", "Cell type"))

selected_columns <- data.frame(
  field = c("cell_id", "sample", "class", "authors_cell_type"),
  selected_column = c(col_cell, col_sample, col_class, col_type),
  stringsAsFactors = FALSE
)
write.csv(selected_columns, file.path(out_dir, "validation_selected_columns.csv"), row.names = FALSE)
print(selected_columns, row.names = FALSE)

if (!is.na(col_cell)) {
  ids <- as.character(annotation[[col_cell]])
  cat("Cell ID duplicates:", anyDuplicated(ids), "; missing:",
      sum(is.na(ids) | ids == ""), "\n")
}
if (is.na(col_sample) || is.na(col_class)) {
  cat("\nWARNING: Cannot confidently identify sample and/or tissue columns.\n")
  cat("Inspect validation_annotation_schema.csv and the metadata before further analysis.\n")
}

# Sample and patient audit. Strictly map KULxx-T/B/N only when identifiers confirm it.
if (!is.na(col_sample)) {
  sample_raw <- trimws(as.character(annotation[[col_sample]]))
  sample_id <- sample_raw
  if (!is.na(col_class)) {
    tissue_raw <- trimws(as.character(annotation[[col_class]]))
  } else tissue_raw <- rep(NA_character_, length(sample_id))
  # Patient id can only be assigned when the KUL## pattern is present in sample ID.
  patient_id <- ifelse(grepl("KUL[0-9]+", sample_id, ignore.case = TRUE),
                       toupper(sub(".*(KUL[0-9]+).*", "\\1", sample_id, ignore.case = TRUE)),
                       NA_character_)
  # Class is displayed as-is; condition is only mapped for unambiguous values.
  canonical_tissue <- rep(NA_character_, length(tissue_raw))
  canonical_tissue[tolower(tissue_raw) %in% c("normal", "n", "normal mucosa")] <- "Normal"
  canonical_tissue[tolower(tissue_raw) %in% c("tumor", "tumour", "t", "tumor core", "tumour core")] <- "Tumor_Core"
  canonical_tissue[tolower(tissue_raw) %in% c("border", "b", "tumor border", "tumour border")] <- "Tumor_Border"
  # Use suffix if Class has no recognized mapping (and only when it is unambiguous).
  suffix <- toupper(sub("^.*-([TBN])$", "\\1", sample_id))
  suffix[!grepl("-[TBN]$", sample_id, ignore.case = TRUE)] <- NA_character_
  canonical_tissue[which(is.na(canonical_tissue) & !is.na(suffix) & suffix == "T")] <- "Tumor_Core"
  canonical_tissue[which(is.na(canonical_tissue) & !is.na(suffix) & suffix == "B")] <- "Tumor_Border"
  canonical_tissue[which(is.na(canonical_tissue) & !is.na(suffix) & suffix == "N")] <- "Normal"
  # Conflict flag: metadata Class and suffix disagree (requires human review).
  suffix_condition <- ifelse(suffix == "T", "Tumor_Core", ifelse(suffix == "B", "Tumor_Border", ifelse(suffix == "N", "Normal", NA_character_)))
  class_condition <- rep(NA_character_, length(tissue_raw))
  class_condition[tolower(tissue_raw) %in% c("normal", "n", "normal mucosa")] <- "Normal"
  class_condition[tolower(tissue_raw) %in% c("tumor", "tumour", "t", "tumor core", "tumour core")] <- "Tumor_Core"
  class_condition[tolower(tissue_raw) %in% c("border", "b", "tumor border", "tumour border")] <- "Tumor_Border"
  conflicts <- !is.na(class_condition) & !is.na(suffix_condition) & class_condition != suffix_condition
  cat("Tissue/suffix mapping conflicts:", sum(conflicts), "\n")
  if (any(conflicts)) canonical_tissue[conflicts] <- NA_character_

  cell_type_raw <- if (!is.na(col_type)) as.character(annotation[[col_type]]) else rep(NA_character_, nrow(annotation))
  review <- data.frame(sample_id = sample_id, patient_id = patient_id,
                       class_original = tissue_raw, tissue_region = canonical_tissue,
                       authors_cell_type = cell_type_raw,
                       mapping_conflict = conflicts, stringsAsFactors = FALSE)
  # Never export single-cell raw data; use grouped counts instead.
  sample_keys <- unique(review[, c("sample_id", "patient_id", "class_original", "tissue_region", "mapping_conflict")])
  sample_keys$n_cells <- vapply(seq_len(nrow(sample_keys)), function(i) {
    sum(review$sample_id == sample_keys$sample_id[i] &
          !is.na(review$sample_id), na.rm = TRUE)
  }, integer(1))
  sample_keys <- sample_keys[order(sample_keys$sample_id), , drop = FALSE]
  write.csv(sample_keys, file.path(out_dir, "validation_sample_mapping_review.csv"), row.names = FALSE)
  cat("\n================ PATIENT & TISSUE SUMMARY ================\n")
  print(sample_keys, row.names = FALSE)
  cat("Distinct inferred patients:", length(unique(stats::na.omit(patient_id))), "\n")
  cat("Unmapped patient IDs:", sum(is.na(patient_id)), "cells\n")
  cat("Unmapped/ambiguous tissue regions:", sum(is.na(canonical_tissue)), "cells\n")
  if (!is.na(col_type)) {
    label_counts <- sort(table(cell_type_raw, useNA = "ifany"), decreasing = TRUE)
    labels <- data.frame(authors_cell_type = names(label_counts), n_cells = as.integer(label_counts))
    write.csv(labels, file.path(out_dir, "validation_authors_cell_type_counts.csv"), row.names = FALSE)
    cat("\n================ AUTHOR CELL TYPE LABELS ================\n")
    print(labels, row.names = FALSE)
  } else {
    cat("\nWARNING: No unambiguous author-cell-type column selected. Check schema.\n")
  }
}

# Freeze discovery signature exactly as-is; does not adapt to validation outcomes.
discovery_manifest <- file.path("results", "tables", "programs", "candidate_signature_manifest_v1.csv")
if (!file.exists(discovery_manifest)) {
  stop("Missing discovery signature manifest: ", discovery_manifest,
       "\nCopy the EXISTING Step 05 manifest into this location; do not reinvent its genes.")
}
manifest <- read.csv(discovery_manifest, stringsAsFactors = FALSE)
stopifnot(all(c("signature_version", "program", "gene", "method", "max_rank") %in% colnames(manifest)))
write.csv(manifest, file.path(out_dir, "frozen_discovery_signature_manifest_v1.csv"), row.names = FALSE)
writeLines(paste("MD5 of discovery manifest:", unname(tools::md5sum(discovery_manifest))),
           file.path(out_dir, "discovery_manifest_md5.txt"))
cat("\n================ FROZEN SIGNATURE ================\n")
print(as.data.frame(table(manifest$program)))
cat("Methods:", paste(unique(manifest$method), collapse = ", "),
    " | max.rank:", paste(unique(manifest$max_rank), collapse = ", "), "\n")

cat("\n========== STEP 07A COMPLETE ==========\n")
cat("CSV reports:", normalizePath(out_dir, winslash = "/"), "\n")
cat("Raw files:", normalizePath(raw_dir, winslash = "/"), "\n")
cat("NEXT: Inspect sample mapping, cell-type labels and count matrix layout before validation subclustering.\n")
cat("IMPORTANT: 6 patients, not 18 independent patients. Tumor Core and Border are repeated regions.\n")
