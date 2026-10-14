# CRC-ImmuneMap | Step 07B | GSE144735 RAW UMI sparse import and feasibility QC
# Run AFTER Step 07A. R >= 4.2; Seurat v5, Matrix, dplyr (no Python needed).
# Converts the supplied wide, gzipped gene x cell TSV into a sparse MatrixMarket
# file IN STREAMING CHUNKS; never allocates the complete dense count matrix.
# Saves a raw Seurat checkpoint and audit tables. Does NOT yet annotate CD8,
# macrophage or fibroblast subtypes, score UCell or perform hypothesis testing.

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
if (!dir.exists(project_dir)) stop("Missing project folder: edit project_dir at script top")
setwd(project_dir)

for (package in c("Seurat", "Matrix", "dplyr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing R package: ", package, ". Install it first.")
  }
}
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

raw_dir <- file.path("data", "raw", "GSE144735")
processed_dir <- file.path("data", "processed", "GSE144735")
out_dir <- file.path("results", "tables", "validation")
object_dir <- "objects"
for (folder in c(processed_dir, out_dir, object_dir)) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
}
annotation_path <- file.path(raw_dir, "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz")
counts_path <- file.path(raw_dir, "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz")
manifest_path <- file.path(out_dir, "frozen_discovery_signature_manifest_v1.csv")
checkpoint_path <- file.path(object_dir, "crc_validation_gse144735_raw_imported.rds")
for (path in c(annotation_path, counts_path, manifest_path)) {
  if (!file.exists(path)) stop("Required input missing: ", path, "\nRun Step 07A first.")
}

annotation <- read.delim(annotation_path, check.names = FALSE,
                         stringsAsFactors = FALSE, quote = "", comment.char = "")
required_cols <- c("Index", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype")
if (!all(required_cols %in% colnames(annotation))) {
  stop("Unexpected annotation columns. Missing: ",
       paste(setdiff(required_cols, colnames(annotation)), collapse = ", "))
}
if (anyDuplicated(annotation$Index) || anyNA(annotation$Index)) {
  stop("Annotation cell IDs are duplicated or missing")
}
if (nrow(annotation) != 27414L) stop("Expected 27,414 cells; found ", nrow(annotation))

annotation$Index <- as.character(annotation$Index)
annotation$Sample <- as.character(annotation$Sample)
annotation$Patient <- as.character(annotation$Patient)
annotation$Cell_type <- as.character(annotation$Cell_type)
annotation$Cell_subtype <- as.character(annotation$Cell_subtype)

suffix <- toupper(sub("^.*-([TBN])$", "\\1", annotation$Sample))
if (any(!suffix %in% c("T", "B", "N"))) stop("Unmapped sample-region suffix")
tissue_region <- c(T = "Tumor_Core", B = "Tumor_Border", N = "Normal")[suffix]
class_map <- c(Tumor = "Tumor_Core", Border = "Tumor_Border", Normal = "Normal")
if (any(is.na(class_map[annotation$Class])) ||
    any(unname(class_map[annotation$Class]) != unname(tissue_region))) {
  stop("Class / Sample suffix mismatch. Inspect annotation before proceeding.")
}
if (!all(annotation$Patient == sub("-[TBN]$", "", annotation$Sample))) {
  stop("Patient and Sample IDs disagree")
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest <- c("signature_version", "program", "gene", "method", "max_rank")
if (!all(required_manifest %in% names(manifest))) stop("Invalid frozen manifest columns")
if (length(unique(manifest$method)) != 1L || unique(manifest$method) != "UCell" ||
    length(unique(manifest$max_rank)) != 1L || unique(manifest$max_rank) != 1500L) {
  stop("Frozen signature method/max_rank differ from discovery")
}
if (anyDuplicated(manifest[, c("program", "gene")])) stop("Duplicated signature gene")
write.csv(manifest, file.path(out_dir, "validation_frozen_manifest_checked.csv"), row.names = FALSE)

# Build coordinate-format MatrixMarket file by reading a few whole gene-rows at a time.
# Use raw integer UMI counts, not TPM/lognormalized expression.
# Input dimensions: gene rows x 27,414 cell columns; first column contains gene names.
stream_count_matrix <- function(path, expected_cell_ids, mtx_path, genes_path, chunk_rows = 20L) {
  body_path <- paste0(mtx_path, ".body.partial")
  tmp_matrix <- paste0(mtx_path, ".partial")
  if (file.exists(body_path)) unlink(body_path)
  if (file.exists(tmp_matrix)) unlink(tmp_matrix)
  con <- gzfile(path, open = "rt")
  body <- file(body_path, open = "wt")
  input_open <- TRUE
  body_open <- TRUE
  on.exit({
    if (input_open) try(close(con), silent = TRUE)
    if (body_open) try(close(body), silent = TRUE)
  }, add = TRUE)

  header <- readLines(con, n = 1L, warn = FALSE)
  if (length(header) != 1L) stop("Empty count matrix")
  count_columns <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (count_columns[1] != "Index") stop("First matrix field should be Index")
  count_cells <- count_columns[-1]
  if (!identical(count_cells, expected_cell_ids)) {
    stop("Count header cells do not match annotation order. Stop before importing.")
  }
  n_cols <- length(count_cells)
  gene_names <- character()
  n_rows <- 0L
  n_nonzero <- 0
  cat("Streaming wide TSV to sparse MatrixMarket (may take several minutes)...\n")
  repeat {
    lines <- readLines(con, n = chunk_rows, warn = FALSE)
    if (!length(lines)) break
    parsed <- strsplit(lines, "\t", fixed = TRUE)
    triplets <- vector("list", length(parsed))
    for (k in seq_along(parsed)) {
      fields <- parsed[[k]]
      if (length(fields) != n_cols + 1L) {
        stop("Malformed count row ", n_rows + k, ": ", length(fields),
             " fields; expected ", n_cols + 1L)
      }
      gene <- fields[1]
      if (is.na(gene) || !nzchar(gene)) stop("Missing gene name at row ", n_rows + k)
      vals <- suppressWarnings(as.integer(fields[-1]))
      if (anyNA(vals) || any(vals < 0L)) {
        stop("Noninteger, missing or negative count at gene row ", n_rows + k)
      }
      pos <- which(vals > 0L)
      triplets[[k]] <- if (length(pos)) {
        paste(n_rows + k, pos, vals[pos])
      } else character()
      gene_names <- c(gene_names, gene)
      n_nonzero <- n_nonzero + length(pos)
    }
    out <- unlist(triplets, use.names = FALSE)
    if (length(out)) writeLines(out, body, sep = "\n")
    n_rows <- n_rows + length(parsed)
    if (n_rows %% 1000L == 0L) {
      cat(sprintf("  %d genes read; %.0f nonzero entries\n", n_rows, n_nonzero))
      flush(body)
    }
  }
  close(con); input_open <- FALSE
  close(body); body_open <- FALSE
  if (n_rows < 1000L || n_nonzero == 0) stop("Unexpectedly empty/small count matrix")
  if (anyDuplicated(gene_names)) {
    cat("WARNING: duplicated gene symbols: ", sum(duplicated(gene_names)),
        "; applying make.unique() only to duplicated names.\n", sep = "")
  }
  genes_unique <- make.unique(gene_names)
  writeLines(genes_unique, genes_path)

  # Second pass copies the temporary coordinates behind the MatrixMarket header.
  # This keeps peak RAM low even when the plain-text count table is huge.
  output <- file(tmp_matrix, open = "wb")
  writeChar(paste0("%%MatrixMarket matrix coordinate integer general\n",
                   "% CRC-ImmuneMap GSE144735 raw UMI counts\n",
                   n_rows, " ", n_cols, " ", format(n_nonzero, scientific = FALSE), "\n"),
            output, eos = NULL, useBytes = TRUE)
  input_body <- file(body_path, open = "rb")
  repeat {
    buffer <- readBin(input_body, what = "raw", n = 1024L * 1024L)
    if (!length(buffer)) break
    writeBin(buffer, output)
  }
  close(input_body)
  close(output)
  if (!file.rename(tmp_matrix, mtx_path)) stop("Could not finalize sparse matrix file")
  unlink(body_path)
  cat("Sparse conversion complete:", n_rows, "genes x", n_cols, "cells;",
      format(n_nonzero, scientific = FALSE), "nonzeros\n")
  invisible(list(n_genes = n_rows, n_cells = n_cols, nnz = n_nonzero))
}

if (file.exists(checkpoint_path)) {
  cat("Existing Seurat checkpoint found, loading:", checkpoint_path, "\n")
  validation <- readRDS(checkpoint_path)
} else {
  con_header <- gzfile(counts_path, open = "rt")
  header <- readLines(con_header, n = 1L, warn = FALSE)
  close(con_header)
  count_cells <- strsplit(header, "\t", fixed = TRUE)[[1]][-1]
  if (length(count_cells) != nrow(annotation) || anyDuplicated(count_cells)) {
    stop("Count matrix cell headers failed uniqueness/count checks")
  }
  # Preserve EXACT original count order; annotation is matched by Index.
  match_idx <- match(count_cells, annotation$Index)
  if (anyNA(match_idx) || !setequal(count_cells, annotation$Index)) {
    stop("Counts and annotation cell IDs do not match")
  }
  ann_sorted <- annotation[match_idx, , drop = FALSE]
  meta <- data.frame(
    row.names = ann_sorted$Index,
    patient_id = ann_sorted$Patient,
    sample_id = ann_sorted$Sample,
    class_original = ann_sorted$Class,
    tissue_region = unname(c(T = "Tumor_Core", B = "Tumor_Border", N = "Normal")[
      toupper(sub("^.*-([TBN])$", "\\1", ann_sorted$Sample))]),
    authors_cell_type = ann_sorted$Cell_type,
    authors_cell_subtype = ann_sorted$Cell_subtype,
    stringsAsFactors = FALSE
  )
  mtx_path <- file.path(processed_dir, "validation_raw_counts.mtx")
  genes_path <- file.path(processed_dir, "validation_genes.txt")
  if (!file.exists(mtx_path) || !file.exists(genes_path)) {
    stream_count_matrix(counts_path, count_cells, mtx_path, genes_path)
  }
  genes <- readLines(genes_path, warn = FALSE)
  cat("Loading sparse MatrixMarket counts into Matrix package...\n")
  counts <- methods::as(Matrix::readMM(mtx_path), "dgCMatrix")
  if (nrow(counts) != length(genes) || ncol(counts) != length(count_cells)) {
    stop("Sparse matrix dimension mismatch. Delete intermediate .mtx/.txt and rerun.")
  }
  rownames(counts) <- genes
  colnames(counts) <- count_cells
  stopifnot(identical(colnames(counts), rownames(meta)))
  validation <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "CRC_GSE144735",
    assay = "RNA",
    meta.data = meta,
    min.cells = 0L,
    min.features = 0L
  )
  rm(counts)
  invisible(gc())
  mt_genes <- grep("^MT-", rownames(validation), value = TRUE, ignore.case = TRUE)
  if (length(mt_genes)) {
    validation[["percent.mt"]] <- Seurat::PercentageFeatureSet(
      validation, assay = "RNA", features = mt_genes)
  } else {
    warning("No MT- genes detected; percent.mt set to NA, QC cannot use it.")
    validation$percent.mt <- rep(NA_real_, ncol(validation))
  }
  cat("Saving Seurat raw-import checkpoint (this can take a moment)...\n")
  saveRDS(validation, checkpoint_path, compress = TRUE)
}

# Audit the checkpoint and its annotation. No filtering/renormalization here.
DefaultAssay(validation) <- "RNA"
meta <- validation[[]]
if (ncol(validation) != 27414L ||
    !identical(colnames(validation), rownames(meta)) ||
    anyNA(meta$patient_id) || anyNA(meta$sample_id) ||
    !all(c("Tumor_Core", "Tumor_Border", "Normal") %in% unique(meta$tissue_region))) {
  stop("Checkpoint verification failed; inspect source and checkpoint")
}
expected_types <- c("Stromal cells", "Epithelial cells", "T cells",
                    "B cells", "Myeloids", "Mast cells")
if (!setequal(unique(meta$authors_cell_type), expected_types)) {
  stop("Unexpected author-defined broad cell types; inspect annotation")
}

# Original discovery QC thresholds displayed as reference, NOT applied to validation.
meta$qc_pass_discovery_reference <- with(meta,
  nFeature_RNA > 200 & nFeature_RNA < 6000 &
    nCount_RNA > 1000 & percent.mt < 20
)
validation$qc_pass_discovery_reference <- meta$qc_pass_discovery_reference

qc_summary <- meta %>%
  group_by(patient_id, sample_id, tissue_region) %>%
  summarise(
    n_cells = dplyr::n(),
    median_features = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_percent_mt = if (all(is.na(percent.mt))) NA_real_ else median(percent.mt, na.rm = TRUE),
    n_pass_reference = sum(qc_pass_discovery_reference, na.rm = TRUE),
    pct_pass_reference = round(100 * mean(qc_pass_discovery_reference, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>% arrange(patient_id, tissue_region)

lineage_by_sample <- meta %>%
  count(patient_id, sample_id, tissue_region, authors_cell_type, name = "n_cells") %>%
  arrange(patient_id, sample_id, authors_cell_type)

subtype_by_lineage <- meta %>%
  count(authors_cell_type, authors_cell_subtype, tissue_region, name = "n_cells") %>%
  arrange(authors_cell_type, desc(n_cells))

coverage <- manifest %>%
  mutate(present_in_validation = gene %in% rownames(validation))
if (any(!coverage$present_in_validation)) {
  cat("WARNING: Missing frozen gene(s): ",
      paste(coverage$gene[!coverage$present_in_validation], collapse = ", "), "\n")
  cat("DO NOT silently replace/drop genes; document coverage before scores.\n")
}

# Candidate base lineage coverage only. No claims about subtype identity yet.
base_lineages <- c("Myeloids", "Stromal cells", "T cells")
base_coverage <- lineage_by_sample %>%
  filter(authors_cell_type %in% base_lineages) %>%
  mutate(base_lineage = recode(authors_cell_type,
                              "Myeloids" = "Myeloids_unresolved",
                              "Stromal cells" = "Stromal_unresolved",
                              "T cells" = "T_cells_unresolved")) %>%
  select(patient_id, sample_id, tissue_region, base_lineage, n_cells)

write.csv(qc_summary, file.path(out_dir, "validation_07b_qc_by_sample.csv"), row.names = FALSE)
write.csv(lineage_by_sample, file.path(out_dir, "validation_07b_broad_lineage_by_sample.csv"), row.names = FALSE)
write.csv(subtype_by_lineage, file.path(out_dir, "validation_07b_author_subtypes.csv"), row.names = FALSE)
write.csv(coverage, file.path(out_dir, "validation_07b_frozen_signature_gene_coverage.csv"), row.names = FALSE)
write.csv(base_coverage, file.path(out_dir, "validation_07b_base_lineage_coverage.csv"), row.names = FALSE)

cat("\n========== STEP 07B SUMMARY ==========\n")
cat("Raw-count Seurat object:", nrow(validation), "genes x", ncol(validation), "cells\n")
cat("Patients:", n_distinct(meta$patient_id),
    "| samples:", n_distinct(meta$sample_id), "\n")
cat("Distinct features from the frozen signatures available:",
    sum(coverage$present_in_validation), "/", nrow(coverage), "\n")
cat("\nSAMPLE-LEVEL QC (DISCOVERY THRESHOLDS FOR REFERENCE ONLY):\n")
print(qc_summary, n = Inf, width = Inf)
cat("\nBROAD CELL-TYPE COVERAGE BY REGION:\n")
print(meta %>% count(tissue_region, authors_cell_type, name = "n_cells") %>%
        arrange(tissue_region, authors_cell_type), n = Inf, width = Inf)
cat("\nAUTHOR SUBTYPE DISTRIBUTION:\n")
print(subtype_by_lineage, n = Inf, width = Inf)
cat("\nGENES MISSING FROM FROZEN SIGNATURE:\n")
print(coverage %>% filter(!present_in_validation), n = Inf, width = Inf)
cat("\nCHECKPOINT:", normalizePath(checkpoint_path, winslash = "/"), "\n")
cat("REPORTS:", normalizePath(out_dir, winslash = "/"), "\n")
cat("NEXT: Review QC and subcluster Myeloids, Stromal and T cells; annotate specific macrophages, fibroblasts and CD8 prior to frozen scoring.\n")
cat("IMPORTANT: No inferential validation has been run. Six unique patients; Core and Border are repeated regions.\n")
