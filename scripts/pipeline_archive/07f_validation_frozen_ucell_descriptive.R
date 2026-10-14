# CRC-ImmuneMap | Step 07F | Validation frozen-signature UCell scoring
# Uses Step 07E pre-score candidate cell lists, and a PRE-SPECIFIED sensitivity review.
# Descriptive only; DO NOT call this a successful independent validation.
# It does not overwrite any source RDS, earlier results, or 07E files.
# GSE144735 has six patients, and very limited 3-compartment tumor coverage.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  if (!dir.exists(project_dir)) stop("Cannot find project directory: ", project_dir)
  setwd(project_dir)
  pkgs <- c("Seurat", "Matrix", "UCell", "dplyr", "tidyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing package(s): ", paste(missing_pkgs, collapse = ", "))
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
  })

  in_dir <- "results/tables/validation/07e_purity_coverage"
  out_dir <- "results/tables/validation/07f_frozen_scoring"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mf_path <- "results/tables/validation/frozen_discovery_signature_manifest_v1.csv"
  mf_origin <- "results/tables/programs/candidate_signature_manifest_v1.csv"
  stopifnot(file.exists(mf_path), file.exists(mf_origin))
  manifest <- read.csv(mf_path, stringsAsFactors = FALSE)
  original <- read.csv(mf_origin, stringsAsFactors = FALSE)
  req_mf <- c("signature_version", "program", "gene", "method", "max_rank")
  stopifnot(all(req_mf %in% names(manifest)), identical(manifest[req_mf], original[req_mf]))
  if (anyDuplicated(manifest[c("program", "gene")])) stop("Duplicate signature gene")
  if (!all(manifest$method == "UCell") || !all(manifest$max_rank == 1500)) {
    stop("Frozen method/max_rank differ from discovery. Stop and investigate.")
  }
  expected_programs <- c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")
  if (!setequal(unique(manifest$program), expected_programs)) stop("Unexpected frozen program names")
  signatures <- split(manifest$gene, factor(manifest$program, levels = expected_programs))
  signatures <- signatures[expected_programs]
  write.csv(manifest, file.path(out_dir, "signature_manifest_unchanged.csv"), row.names = FALSE)

  cfg <- list(
    myeloids = list(path = "objects/crc_validation_07c_myeloids_clustered_unreviewed.rds",
                    lineage = "Macrophages", programs = "SPP1_TAM",
                    epithelial = c("EPCAM", "KRT8", "KRT18"),
                    drop_clusters = c("0", "1")),
    stromal = list(path = "objects/crc_validation_07c_stromal_clustered_unreviewed.rds",
                   lineage = "Fibroblasts", programs = "FAP_CAF",
                   epithelial = c("EPCAM", "KRT8", "KRT18"),
                   drop_clusters = c("0", "12")),
    tcells = list(path = "objects/crc_validation_07c_tcells_clustered_unreviewed.rds",
                  lineage = "CD8", programs = c("CD8_Dysfunction", "CD8_Cytotoxicity"),
                  gammadelta = c("TRDC", "TRGC1", "TRGC2"),
                  myeloid = c("LYZ", "LST1", "TYROBP"),
                  drop_clusters = "10")
  )
  variant_desc <- data.frame(
    variant = c("primary_07e", "crosslineage_sensitivity", "omit_high_risk_clusters", "all_published_cd8_sensitivity"),
    description = c(
      "07E candidate primary clusters; CD8 also requires author CD8 subtype; NOT final confirmed identities",
      "Same primary, except exclude cells with >=2 epithelial panel genes (myeloid/fibro) or >=2 gamma-delta/myeloid panel genes (CD8); possible ambient RNA means this is NOT ground truth",
      "Same primary except omit macrophage clusters 0 and 1, fibroblast 0 and 12, and CD8 cluster 10; intentionally severe stress test",
      "Same selection as primary for macrophages and fibroblasts; all author-labeled CD8+ T cells in T-cell checkpoint, not just candidate clusters"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(variant_desc, file.path(out_dir, "selection_variant_manifest.csv"), row.names = FALSE)

  score_parts <- list()
  diag_parts <- list()
  counter <- 0L
  for (key in names(cfg)) {
    info <- cfg[[key]]
    input_csv <- file.path(in_dir, paste0(key, "_candidate_cells_NOT_FINAL.csv"))
    if (!file.exists(input_csv)) stop("Missing 07E candidate file: ", input_csv)
    if (!file.exists(info$path)) stop("Missing 07C RDS: ", info$path)
    cand <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("cell_id", "candidate_tier", "validation_cluster", "patient_id") %in% names(cand)))
    if (anyDuplicated(cand$cell_id)) stop(key, ": repeated 07E candidate cell ID")
    cand <- cand[cand$candidate_tier == "primary_candidate", , drop = FALSE]
    if (!nrow(cand)) stop(key, ": no primary candidate cells")

    cat("\n=== Step 07F: ", key, " (loading one RDS only) ===\n", sep = "")
    obj <- readRDS(info$path)
    stopifnot(inherits(obj, "Seurat"))
    DefaultAssay(obj) <- "RNA"
    md <- obj[[]]
    req_md <- c("validation_cluster", "patient_id", "sample_id", "tissue_region", "authors_cell_subtype")
    if (!all(req_md %in% names(md))) stop(key, ": missing metadata")
    if (!identical(rownames(md), colnames(obj))) stop(key, ": metadata alignment problem")
    if (!all(cand$cell_id %in% rownames(md))) stop(key, ": unknown 07E candidate cell(s)")
    md$cell_id <- rownames(md)
    md$validation_cluster <- as.character(md$validation_cluster)
    md$authors_cell_subtype <- as.character(md$authors_cell_subtype)
    check <- md[cand$cell_id, , drop = FALSE]
    if (!all(as.character(check$patient_id) == cand$patient_id) ||
        !all(check$validation_cluster == cand$validation_cluster)) {
      stop(key, ": 07E candidates disagree with RDS metadata")
    }
    primary_ids <- cand$cell_id
    if (key == "tcells") {
      if (any(check$authors_cell_subtype != "CD8+ T cells" | is.na(check$authors_cell_subtype))) {
        stop("Primary CD8 candidates not consistent with 07E subtype gate")
      }
      published_ids <- md$cell_id[!is.na(md$authors_cell_subtype) &
                                  md$authors_cell_subtype == "CD8+ T cells"]
      selected_ids <- union(primary_ids, published_ids)
    } else {
      selected_ids <- primary_ids
    }
    selected <- md[selected_ids, , drop = FALSE]
    selected$primary_mask <- selected$cell_id %in% primary_ids
    if (key == "tcells") {
      selected$published_cd8_mask <- selected$cell_id %in% published_ids
    } else {
      selected$published_cd8_mask <- FALSE
    }

    full_counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
    if (!identical(colnames(full_counts), md$cell_id)) stop("Count matrix cell order mismatch")
    these_signatures <- signatures[info$programs]
    signature_genes <- unique(unlist(these_signatures, use.names = FALSE))
    if (!all(signature_genes %in% rownames(full_counts))) {
      stop("Frozen signature genes missing in ", key, ": ",
           paste(setdiff(signature_genes, rownames(full_counts)), collapse = ", "))
    }
    panel_detect <- function(genes) {
      if (!all(genes %in% rownames(full_counts))) stop("Missing diagnostic gene(s): ", paste(setdiff(genes, rownames(full_counts)), collapse = ", "))
      as.integer(Matrix::colSums(full_counts[genes, selected_ids, drop = FALSE] > 0))
    }
    if (key %in% c("myeloids", "stromal")) {
      selected$crosslineage_flag <- panel_detect(info$epithelial) >= 2L
    } else {
      selected$crosslineage_flag <- panel_detect(info$gammadelta) >= 2L |
                                  panel_detect(info$myeloid) >= 2L
    }
    selected$strict_mask <- selected$primary_mask & !selected$crosslineage_flag
    selected$omit_risk_mask <- selected$primary_mask &
      !selected$validation_cluster %in% info$drop_clusters
    selected$all_published_cd8_mask <- if (key == "tcells")
      selected$published_cd8_mask else selected$primary_mask

    # IMPORTANT: UCell ranks genes WITHIN each cell. Keep ALL features and only
    # select columns, never score using just the signature genes.
    counts <- full_counts[, selected_ids, drop = FALSE]
    score_mat <- UCell::ScoreSignatures_UCell(
      counts, features = these_signatures, maxRank = 1500L, ncores = 1L
    )
    score_mat <- as.data.frame(score_mat)
    if (!setequal(rownames(score_mat), selected_ids)) stop("UCell returned unexpected cell IDs")
    score_mat <- score_mat[selected_ids, , drop = FALSE]
    required_scores <- paste0(info$programs, "_UCell")
    if (!all(required_scores %in% colnames(score_mat))) stop("Unexpected score names: ", paste(colnames(score_mat), collapse = ", "))
    if (any(!is.finite(as.matrix(score_mat[, required_scores, drop = FALSE])))) stop("Invalid scores")

    flags <- data.frame(
      cell_id = selected_ids, lineage = info$lineage,
      patient_id = as.character(selected$patient_id),
      sample_id = as.character(selected$sample_id),
      tissue_region = as.character(selected$tissue_region),
      validation_cluster = selected$validation_cluster,
      authors_cell_subtype = selected$authors_cell_subtype,
      primary_mask = selected$primary_mask,
      crosslineage_flag = selected$crosslineage_flag,
      strict_mask = selected$strict_mask,
      omit_risk_mask = selected$omit_risk_mask,
      all_published_cd8_mask = selected$all_published_cd8_mask,
      stringsAsFactors = FALSE
    )
    if (anyNA(flags[, c("cell_id", "patient_id", "sample_id", "tissue_region")])) stop(key, ": missing metadata")
    diag_parts[[key]] <- flags
    score_long <- cbind(flags, score_mat[, required_scores, drop = FALSE]) |>
      tidyr::pivot_longer(cols = dplyr::all_of(required_scores),
                          names_to = "program", values_to = "score")
    score_long$program <- sub("_UCell$", "", score_long$program)
    score_parts[[key]] <- score_long
    cat("Primary candidate cells: ", sum(flags$primary_mask),
        " | crosslineage flags in primary: ", sum(flags$primary_mask & flags$crosslineage_flag),
        " | programs: ", paste(info$programs, collapse = ", "), "\n", sep = "")
    rm(obj, full_counts, counts, score_mat, score_long, selected, md, cand)
    invisible(gc())
  }

  all_scores <- dplyr::bind_rows(score_parts)
  cell_flags <- dplyr::bind_rows(diag_parts)
  write.csv(cell_flags, file.path(out_dir, "scored_cell_selection_flags.csv"), row.names = FALSE)
  write.csv(all_scores, file.path(out_dir, "all_scored_cells_with_flags.csv"), row.names = FALSE)

  variants <- list(
    primary_07e = "primary_mask",
    crosslineage_sensitivity = "strict_mask",
    omit_high_risk_clusters = "omit_risk_mask",
    all_published_cd8_sensitivity = "all_published_cd8_mask"
  )
  scores_by_variant <- dplyr::bind_rows(lapply(names(variants), function(v) {
    ans <- all_scores[all_scores[[variants[[v]]]], , drop = FALSE]
    ans$variant <- v
    ans
  }))
  sample_scores <- scores_by_variant |>
    dplyr::group_by(variant, patient_id, sample_id, tissue_region, lineage, program) |>
    dplyr::summarise(n_cells = dplyr::n(), mean_score = mean(score),
                     median_score = median(score), .groups = "drop")
  patient_region <- scores_by_variant |>
    dplyr::group_by(variant, patient_id, tissue_region, lineage, program) |>
    dplyr::summarise(n_cells = dplyr::n(), n_samples = dplyr::n_distinct(sample_id),
                     mean_score = mean(score), median_score = median(score), .groups = "drop")
  # Tumor Core and Border constitute TWO regions from ONE patient. Combining
  # them is secondary and does not create independent replicates.
  combined <- scores_by_variant |>
    dplyr::filter(tissue_region %in% c("Tumor_Core", "Tumor_Border")) |>
    dplyr::group_by(variant, patient_id, lineage, program) |>
    dplyr::summarise(n_cells = dplyr::n(), n_samples = dplyr::n_distinct(sample_id),
                     mean_score = mean(score), median_score = median(score), .groups = "drop") |>
    dplyr::mutate(tissue_region = "Tumor_Core_plus_Border")
  write.csv(sample_scores, file.path(out_dir, "candidate_scores_by_sample.csv"), row.names = FALSE)
  write.csv(patient_region, file.path(out_dir, "candidate_scores_by_patient_region.csv"), row.names = FALSE)
  write.csv(combined, file.path(out_dir, "candidate_scores_core_plus_border_per_patient.csv"), row.names = FALSE)

  # For incomplete patient/region programs, keep n_cells=0 and score NA; NEVER zero-fill scores.
  all_pats <- sort(unique(cell_flags$patient_id))
  all_regions <- c("Tumor_Core", "Tumor_Border", "Normal", "Tumor_Core_plus_Border")
  expected_programs <- c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")
  all_n <- dplyr::bind_rows(patient_region, combined)
  base_grid <- tidyr::expand_grid(variant = names(variants), patient_id = all_pats,
                                  tissue_region = all_regions, program = expected_programs)
  base_grid <- base_grid |>
    dplyr::left_join(dplyr::select(all_n, variant, patient_id, tissue_region,
                                    program, n_cells, mean_score, median_score),
                     by = c("variant", "patient_id", "tissue_region", "program")) |>
    dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0L),
                  score_ge10 = dplyr::if_else(n_cells >= 10L, mean_score, NA_real_),
                  score_ge20 = dplyr::if_else(n_cells >= 20L, mean_score, NA_real_))
  wide_n <- base_grid |>
    dplyr::select(variant, patient_id, tissue_region, program, n_cells) |>
    tidyr::pivot_wider(names_from = program, values_from = n_cells) |>
    dplyr::mutate(all_three_ge10 = SPP1_TAM >= 10 & FAP_CAF >= 10 & CD8_Dysfunction >= 10,
                  all_three_ge20 = SPP1_TAM >= 20 & FAP_CAF >= 20 & CD8_Dysfunction >= 20) |>
    dplyr::arrange(variant, patient_id, tissue_region)
  totals <- wide_n |>
    dplyr::group_by(variant, tissue_region) |>
    dplyr::summarise(n_patients_ge10 = sum(all_three_ge10),
                     n_patients_ge20 = sum(all_three_ge20), .groups = "drop")
  write.csv(base_grid, file.path(out_dir, "patient_program_coverage_with_missing_scores.csv"), row.names = FALSE)
  write.csv(wide_n, file.path(out_dir, "patient_three_program_coverage.csv"), row.names = FALSE)
  write.csv(totals, file.path(out_dir, "coverage_by_selection_variant.csv"), row.names = FALSE)
  cat("\n========== STEP 07F SUMMARY (DESCRIPTIVE; NO ASSOCIATION TEST) ==========\n")
  print(totals, n = Inf, width = Inf)
  cat("\nOutputs:", normalizePath(out_dir, winslash = "/"), "\n")
  cat("Primary candidate groups are PROVISIONAL. Region-combined sample data are NOT extra patients.\n")
  cat("Do NOT report confirmatory p-values for 2-3 complete-case tumor patients.\n")
})
