# CRC-ImmuneMap | Step 07E | GSE144735 cell-level lineage-purity & candidate-coverage audit
# READ-ONLY relative to source RDS/07C/07D. No clustering, no gene program scoring,
# no final cell labels, no inferential tests. Uses only pre-existing Seurat checkpoints.
# Candidate masks are PROVISIONAL: inspect all outputs BEFORE freezing primary masks.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  if (!dir.exists(project_dir)) stop("Project not found: ", project_dir)
  setwd(project_dir)
  for (pkg in c("Seurat", "Matrix", "dplyr", "tidyr", "tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
  }
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
  })
  input_dir <- "objects"
  output_dir <- "results/tables/validation/07e_purity_coverage"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cfg <- list(
    myeloids = list(
      path = "crc_validation_07c_myeloids_clustered_unreviewed.rds",
      compartment = "Macrophages",
      primary = c("0", "1", "2", "5"),
      review = c("7", "9"),
      panels = list(
        macrophage = c("C1QA", "C1QC", "CD68"),
        monocyte = c("FCN1", "S100A8", "S100A9"),
        dc = c("FCER1A", "CD1C", "CLEC10A"),
        t = c("CD3D", "CD3E", "TRAC"),
        epithelial = c("EPCAM", "KRT8", "KRT18")
      )
    ),
    stromal = list(
      path = "crc_validation_07c_stromal_clustered_unreviewed.rds",
      compartment = "Fibroblasts",
      primary = c("0", "2", "3", "4", "9", "11", "12"),
      review = c("10"),
      panels = list(
        fibroblast = c("COL1A1", "DCN", "LUM"),
        endothelial = c("PECAM1", "VWF", "CLDN5"),
        mural = c("RGS5", "PDGFRB", "MYH11"),
        t = c("CD3D", "CD3E", "TRAC"),
        epithelial = c("EPCAM", "KRT8", "KRT18"),
        myeloid = c("LYZ", "LST1", "TYROBP")
      )
    ),
    tcells = list(
      path = "crc_validation_07c_tcells_clustered_unreviewed.rds",
      compartment = "CD8",
      primary = c("2", "4", "10"),
      review = character(),
      panels = list(
        t = c("CD3D", "CD3E", "TRAC"),
        cd8 = c("CD8A", "CD8B"),
        gammadelta = c("TRDC", "TRGC1", "TRGC2"),
        nk = c("KLRD1", "FCGR3A", "NCR1"),
        epithelial = c("EPCAM", "KRT8", "KRT18"),
        myeloid = c("LYZ", "LST1", "TYROBP")
      )
    )
  )

  manifest <- do.call(rbind, lapply(names(cfg), function(key) {
    x <- cfg[[key]]
    data.frame(
      lineage = key,
      compartment = x$compartment,
      cluster = c(x$primary, x$review),
      candidate_tier = c(rep("primary_candidate", length(x$primary)),
                         rep("sensitivity_review", length(x$review))),
      additional_subtype_gate = if (key == "tcells") "authors_cell_subtype == CD8+ T cells" else "none",
      stringsAsFactors = FALSE
    )
  }))
  write.csv(manifest, file.path(output_dir, "provisional_candidate_rules_not_final.csv"), row.names = FALSE)

  assess <- function(key) {
    conf <- cfg[[key]]
    path <- file.path(input_dir, conf$path)
    if (!file.exists(path)) stop("Missing checkpoint: ", path)
    cat("\n=== Loading ", key, " checkpoint (one at a time) ===\n", sep = "")
    obj <- readRDS(path)
    stopifnot(inherits(obj, "Seurat"))
    DefaultAssay(obj) <- "RNA"
    md <- obj[[]]
    required <- c("validation_cluster", "patient_id", "sample_id", "tissue_region", "authors_cell_subtype")
    if (!all(required %in% names(md))) stop("Missing metadata: ", paste(setdiff(required, names(md)), collapse = ", "))
    if (!identical(rownames(md), colnames(obj))) stop("Cell order disagreement in ", key)
    md$cell_id <- rownames(md)
    md$validation_cluster <- as.character(md$validation_cluster)
    md$authors_cell_subtype <- as.character(md$authors_cell_subtype)
    md$patient_id <- as.character(md$patient_id)
    md$tissue_region <- as.character(md$tissue_region)
    if (anyNA(md[, c("cell_id", "validation_cluster", "patient_id", "tissue_region")])) stop("Missing essential metadata")
    if (anyDuplicated(md$cell_id)) stop("Duplicated cell IDs")
    raw <- GetAssayData(obj, assay = "RNA", layer = "counts")
    needed <- unique(unlist(conf$panels, use.names = FALSE))
    present <- intersect(needed, rownames(raw))
    absent <- setdiff(needed, rownames(raw))
    if (length(absent)) cat("Absent diagnostic genes (not silently replaced): ", paste(absent, collapse = ", "), "\n", sep = "")
    raw <- raw[present, , drop = FALSE]
    if (!identical(colnames(raw), md$cell_id)) stop("Counts / metadata alignment failed")

    diag_data <- data.frame(row.names = md$cell_id)
    for (nm in names(conf$panels)) {
      genes <- intersect(conf$panels[[nm]], rownames(raw))
      if (length(genes) < 2L) stop("Insufficient genes in diagnostic panel: ", key, ":", nm)
      panel_counts <- raw[genes, , drop = FALSE]
      diag_data[[paste0("n_detect_", nm)]] <- as.integer(Matrix::colSums(panel_counts > 0))
      diag_data[[paste0("umis_", nm)]] <- as.numeric(Matrix::colSums(panel_counts))
    }
    diag_data <- diag_data[md$cell_id, , drop = FALSE]
    md <- cbind(md, diag_data)

    md$candidate_tier <- ifelse(md$validation_cluster %in% conf$primary, "primary_candidate",
                              ifelse(md$validation_cluster %in% conf$review, "sensitivity_review", "not_selected"))
    if (key == "tcells") {
      # Prevent gamma-delta / NK / CD4 minority cells in clusters 2/4/10
      # from being swept into the CD8 pool. These are provisional author labels;
      # cellular lineage diagnostics still need review.
      md$candidate_tier[md$candidate_tier == "primary_candidate" &
                          (is.na(md$authors_cell_subtype) |
                           md$authors_cell_subtype != "CD8+ T cells")] <- "not_selected"
    }

    # Diagnostic co-detection and UMI levels. Co-expression is NOT proof of doublets;
    # can also arise from ambient RNA or genuine cellular states.
    diag_names <- names(conf$panels)
    cluster_summary <- md |>
      dplyr::group_by(validation_cluster) |>
      dplyr::summarise(n_cells = dplyr::n(),
                       n_patients = dplyr::n_distinct(patient_id),
                       n_primary_candidates = sum(candidate_tier == "primary_candidate"),
                       n_review_candidates = sum(candidate_tier == "sensitivity_review"),
                       .groups = "drop")
    for (nm in diag_names) {
      ncol <- paste0("n_detect_", nm)
      ucol <- paste0("umis_", nm)
      x <- md |>
        dplyr::group_by(validation_cluster) |>
        dplyr::summarise(pct_panel_2genes = 100 * mean(.data[[ncol]] >= 2),
                         median_panel_umis = median(.data[[ucol]]),
                         .groups = "drop")
      names(x)[names(x) == "pct_panel_2genes"] <- paste0("pct_", nm, "_2genes")
      names(x)[names(x) == "median_panel_umis"] <- paste0("median_", nm, "_umis")
      cluster_summary <- dplyr::left_join(cluster_summary, x, by = "validation_cluster")
    }
    if (key == "myeloids") {
      mix <- md |> dplyr::group_by(validation_cluster) |>
        dplyr::summarise(pct_macrophage_and_epithelial_2genes =
                           100 * mean(n_detect_macrophage >= 2 & n_detect_epithelial >= 2), .groups = "drop")
    } else if (key == "stromal") {
      mix <- md |> dplyr::group_by(validation_cluster) |>
        dplyr::summarise(pct_fibroblast_and_epithelial_2genes =
                           100 * mean(n_detect_fibroblast >= 2 & n_detect_epithelial >= 2), .groups = "drop")
    } else {
      mix <- md |> dplyr::group_by(validation_cluster) |>
        dplyr::summarise(pct_T_and_myeloid_2genes =
                           100 * mean(n_detect_t >= 2 & n_detect_myeloid >= 2),
                         pct_CD8_panel_any_gene = 100 * mean(n_detect_cd8 >= 1),
                         pct_gammadelta_panel_any_gene = 100 * mean(n_detect_gammadelta >= 1),
                         .groups = "drop")
    }
    cluster_summary <- dplyr::left_join(cluster_summary, mix, by = "validation_cluster") |>
      dplyr::arrange(as.integer(validation_cluster))
    write.csv(cluster_summary,
              file.path(output_dir, paste0(key, "_cluster_cell_level_purity_summary.csv")), row.names = FALSE)

    # Persist an *auditable candidate list*, NOT a frozen validation inclusion mask.
    cand <- md[md$candidate_tier != "not_selected", , drop = FALSE]
    if (nrow(cand) == 0) stop("No candidate cells for ", key)
    keep <- c("cell_id", "patient_id", "sample_id", "tissue_region", "validation_cluster",
              "authors_cell_subtype", "candidate_tier")
    write.csv(cand[, keep, drop = FALSE],
              file.path(output_dir, paste0(key, "_candidate_cells_NOT_FINAL.csv")), row.names = FALSE)
    coverage <- cand |>
      dplyr::count(patient_id, tissue_region, candidate_tier, name = "n_cells") |>
      dplyr::mutate(lineage = conf$compartment)
    cat("Audit done for ", key, ": ", nrow(cand), " primary + review cells\n", sep = "")
    rm(obj, raw, diag_data, md)
    invisible(gc(verbose = FALSE))
    list(coverage = coverage)
  }

  coverages <- lapply(names(cfg), assess)
  counts <- dplyr::bind_rows(lapply(coverages, `[[`, "coverage"))
  write.csv(counts, file.path(output_dir, "candidate_counts_by_patient_region_and_tier.csv"), row.names = FALSE)

  pat <- sort(unique(counts$patient_id))
  regions <- c("Tumor_Core", "Tumor_Border", "Normal")
  lineages <- c("Macrophages", "Fibroblasts", "CD8")
  primary_counts <- counts |>
    dplyr::filter(candidate_tier == "primary_candidate") |>
    dplyr::group_by(patient_id, tissue_region, lineage) |>
    dplyr::summarise(n_cells = sum(n_cells), .groups = "drop")
  grid <- tidyr::expand_grid(patient_id = pat, tissue_region = regions, lineage = lineages) |>
    dplyr::left_join(primary_counts, by = c("patient_id", "tissue_region", "lineage")) |>
    dplyr::mutate(n_cells = tidyr::replace_na(n_cells, 0L))
  wide <- grid |>
    tidyr::pivot_wider(names_from = lineage, values_from = n_cells) |>
    dplyr::mutate(all_three_ge10 = Macrophages >= 10 & Fibroblasts >= 10 & CD8 >= 10,
                  all_three_ge20 = Macrophages >= 20 & Fibroblasts >= 20 & CD8 >= 20) |>
    dplyr::arrange(patient_id, factor(tissue_region, levels = regions))
  write.csv(wide, file.path(output_dir, "provisional_primary_coverage_by_patient_region.csv"), row.names = FALSE)

  # Secondary geography: Core+Border summed at PATIENT level, never counted as two patients.
  tumor <- grid |>
    dplyr::filter(tissue_region %in% c("Tumor_Core", "Tumor_Border")) |>
    dplyr::group_by(patient_id, lineage) |>
    dplyr::summarise(n_cells = sum(n_cells), .groups = "drop") |>
    tidyr::pivot_wider(names_from = lineage, values_from = n_cells) |>
    dplyr::mutate(all_three_ge10 = Macrophages >= 10 & Fibroblasts >= 10 & CD8 >= 10,
                  all_three_ge20 = Macrophages >= 20 & Fibroblasts >= 20 & CD8 >= 20)
  write.csv(tumor, file.path(output_dir, "provisional_core_plus_border_patient_coverage.csv"), row.names = FALSE)

  cat("\n========== STEP 07E — PROVISIONAL COVERAGE SUMMARY ==========\n")
  print(wide[, c("patient_id", "tissue_region", "Macrophages", "Fibroblasts", "CD8", "all_three_ge20")], row.names = FALSE)
  cat("\nPatients with all three >=20 in Tumor Core: ",
      sum(wide$all_three_ge20[wide$tissue_region == "Tumor_Core"]), " / ", length(pat), "\n", sep = "")
  cat("Patients with all three >=20 in Tumor Border: ",
      sum(wide$all_three_ge20[wide$tissue_region == "Tumor_Border"]), " / ", length(pat), "\n", sep = "")
  cat("Patients with all three >=20 in Core+Border: ",
      sum(tumor$all_three_ge20), " / ", nrow(tumor), "\n", sep = "")
  cat("\nOutput: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
  cat("IMPORTANT: These are candidate masks only. NO final annotation, NO UCell, NO association.\n")
})
