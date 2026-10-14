# CRC-ImmuneMap | Step 07D: validation annotation evidence audit
# Input: corrected donor reports and 07C marker/metadata tables.
# NO Seurat reload, NO annotation assigned to cells, NO UCell scoring.
# This script does not alter any source CSVs or RDS checkpoints.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  stopifnot(dir.exists(project_dir))
  setwd(project_dir)
  input_dir <- "results/tables/validation/07c"
  donor_dir <- file.path(input_dir, "patient_contribution_fixed")
  out_dir <- file.path(input_dir, "annotation_evidence_review")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  focus <- list(
    myeloids = c("LYZ", "LST1", "FCN1", "S100A8", "C1QA", "C1QC",
                 "CD68", "SPP1", "TREM2", "GPNMB", "FCER1A", "CD1C",
                 "CD3D", "EPCAM"),
    stromal  = c("COL1A1", "DCN", "PDGFRA", "FAP", "CTHRC1",
                 "COL11A1", "POSTN", "ACTA2", "PECAM1", "VWF",
                 "RGS5", "PDGFRB", "SOX10", "KRT8", "CD3D"),
    tcells   = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "GZMK",
                 "NKG7", "PDCD1", "LAG3", "TIGIT", "TOX", "CXCL13",
                 "FOXP3", "IL2RA", "TRDC", "KLRD1", "MKI67", "LYZ")
  )

  # Hypotheses for human review only: not final cell-type annotations.
  review_hints <- list(
    myeloids = c(
      "0" = "SPP1-associated macrophage candidate; donor-dominant",
      "1" = "SPP1-associated macrophage candidate; donor-dominant",
      "2" = "C1Q/FOLR2 macrophage candidate; review",
      "3" = "cDC candidate; review",
      "4" = "FCN1 monocyte candidate; review",
      "5" = "C1Q macrophage candidate; review",
      "6" = "FCN1 monocyte candidate; donor-dominant",
      "7" = "Myeloid state uncertain; review",
      "8" = "T-lineage cross-expression concern; review",
      "9" = "Cycling / rare myeloid cluster; review"
    ),
    stromal = c(
      "0" = "FAP-related fibroblast/CAF candidate; review",
      "10" = "FAP-related candidate; donor + keratin concerns",
      "12" = "FAP-related fibroblast/CAF candidate; review"
    ),
    tcells = c(
      "1" = "Regulatory T-cell candidate; review",
      "2" = "CD8-associated candidate; review",
      "4" = "CD8-associated candidate; review",
      "8" = "Gamma-delta candidate; review",
      "10" = "Checkpoint-associated CD8 candidate; donor-dominant"
    )
  )

  all_overviews <- list()
  for (lineage in names(focus)) {
    cat("\n=====", toupper(lineage), "=====\n")
    data_path <- file.path(input_dir, paste0(lineage, "_cluster_by_patient_sample_subtype.csv"))
    marker_path <- file.path(input_dir, paste0(lineage, "_marker_detection_by_cluster.csv"))
    donor_path <- file.path(donor_dir, paste0(lineage, "_cluster_patient_summary_corrected.csv"))
    input_paths <- c(data_path, marker_path, donor_path)
    if (!all(file.exists(input_paths))) {
      stop("Missing required file(s):\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"),
           "\nRun Step 07C and the donor-fix script before Step 07D.")
    }

    dat <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
    mark <- read.csv(marker_path, stringsAsFactors = FALSE, check.names = FALSE)
    donor <- read.csv(donor_path, stringsAsFactors = FALSE, check.names = FALSE)

    required_dat <- c("validation_cluster", "patient_id", "sample_id", "tissue_region",
                      "authors_cell_subtype", "n_cells")
    required_mark <- c("validation_cluster", "gene", "pct_detected")
    required_donor <- c("validation_cluster", "n_cells", "n_patients",
                        "dominant_patient", "dominant_patient_pct")
    stopifnot(all(required_dat %in% names(dat)),
              all(required_mark %in% names(mark)),
              all(required_donor %in% names(donor)),
              !anyNA(dat[, c("validation_cluster", "patient_id", "sample_id", "n_cells")]),
              all(dat$n_cells >= 0))

    dat$validation_cluster <- as.character(dat$validation_cluster)
    donor$validation_cluster <- as.character(donor$validation_cluster)
    mark$validation_cluster <- as.character(mark$validation_cluster)
    mark$gene <- as.character(mark$gene)

    if (anyDuplicated(mark[c("validation_cluster", "gene")])) {
      stop("Duplicate marker rows for ", lineage)
    }
    clids <- as.character(sort(unique(as.integer(dat$validation_cluster))))
    stopifnot(setequal(clids, donor$validation_cluster),
              setequal(clids, unique(mark$validation_cluster)))

    # Independently verify correction against cluster-by-patient counts.
    dat_total <- stats::aggregate(n_cells ~ validation_cluster, dat, sum)
    matched <- match(donor$validation_cluster, dat_total$validation_cluster)
    stopifnot(!anyNA(matched), all(donor$n_cells == dat_total$n_cells[matched]))

    # Tumor_Core / Tumor_Border / Normal counts and within-cluster percentages.
    regions <- c("Normal", "Tumor_Border", "Tumor_Core")
    region_tbl <- stats::aggregate(n_cells ~ validation_cluster + tissue_region, dat, sum)
    overview <- donor
    for (region in regions) {
      v <- region_tbl$n_cells[region_tbl$tissue_region == region]
      nm <- region_tbl$validation_cluster[region_tbl$tissue_region == region]
      overview[[paste0("n_", region)]] <- v[match(overview$validation_cluster, nm)]
      overview[[paste0("n_", region)]][is.na(overview[[paste0("n_", region)]])] <- 0
      overview[[paste0("pct_", region)]] <- round(100 * overview[[paste0("n_", region)]] / overview$n_cells, 1)
    }
    stopifnot(all(rowSums(overview[, paste0("n_", regions), drop = FALSE]) == overview$n_cells))

    # Author subtype mix as descriptive cross-check only (not gold-standard labels).
    dat$authors_cell_subtype[is.na(dat$authors_cell_subtype) |
                              dat$authors_cell_subtype == ""] <- "Unknown"
    subcounts <- stats::aggregate(n_cells ~ validation_cluster + authors_cell_subtype,
                                  dat, sum)
    subcounts <- subcounts[order(as.integer(subcounts$validation_cluster),
                                 -subcounts$n_cells), , drop = FALSE]
    subcounts$pct_of_cluster <- round(100 * subcounts$n_cells /
                                     overview$n_cells[match(subcounts$validation_cluster,
                                                             overview$validation_cluster)], 1)
    subtype_top <- lapply(overview$validation_cluster, function(k) {
      x <- subcounts[subcounts$validation_cluster == k, , drop = FALSE]
      x <- x[order(-x$n_cells, x$authors_cell_subtype), , drop = FALSE]
      top <- utils::head(x, 3)
      paste(paste0(top$authors_cell_subtype, " (", top$n_cells, "; ",
                   top$pct_of_cluster, "%)"), collapse = " | ")
    })
    overview$top_author_subtypes <- unlist(subtype_top)

    # Report detection frequencies (percent of cells with >=1 raw UMI).
    focus_genes <- focus[[lineage]]
    coverage <- data.frame(lineage = lineage, gene = focus_genes,
                           present_in_marker_report = focus_genes %in% unique(mark$gene))
    write.csv(coverage,
              file.path(out_dir, paste0(lineage, "_focus_marker_coverage.csv")), row.names = FALSE)
    for (gene in focus_genes) {
      x <- mark[mark$gene == gene, c("validation_cluster", "pct_detected"), drop = FALSE]
      overview[[paste0("pct_detect_", gene)]] <- x$pct_detected[
        match(overview$validation_cluster, x$validation_cluster)]
    }

    hints <- review_hints[[lineage]]
    hint <- unname(hints[overview$validation_cluster])
    hint[is.na(hint)] <- paste("Other", lineage, "cluster; inspect markers")
    overview$review_hint_not_final <- hint
    overview$donor_gt70 <- overview$dominant_patient_pct > 70
    overview$review_required <- TRUE
    overview <- overview[order(as.integer(overview$validation_cluster)), , drop = FALSE]

    write.csv(overview, file.path(out_dir, paste0(lineage, "_annotation_evidence_master.csv")),
              row.names = FALSE)
    write.csv(subcounts, file.path(out_dir, paste0(lineage, "_cluster_author_subtype_composition.csv")),
              row.names = FALSE)
    write.csv(region_tbl, file.path(out_dir, paste0(lineage, "_cluster_region_counts.csv")),
              row.names = FALSE)

    print(overview[, c("validation_cluster", "n_cells", "n_patients",
                       "dominant_patient", "dominant_patient_pct",
                       "n_Tumor_Core", "n_Tumor_Border", "n_Normal",
                       "review_hint_not_final")], row.names = FALSE)
    all_overviews[[lineage]] <- data.frame(
      lineage = lineage,
      n_clusters = nrow(overview),
      cells = sum(overview$n_cells),
      donor_dominant_clusters = sum(overview$donor_gt70),
      focus_marker_genes_missing = sum(!coverage$present_in_marker_report)
    )
  }

  full_summary <- do.call(rbind, all_overviews)
  write.csv(full_summary,
            file.path(out_dir, "validation_07d_evidence_audit_summary.csv"), row.names = FALSE)
  cat("\n========== STEP 07D EVIDENCE SUMMARY ==========\n")
  print(full_summary, row.names = FALSE)
  cat("\nReports: ", normalizePath(out_dir, winslash = "/"), "\n", sep = "")
  cat("No final labels or filtered RDS were created.\n")
  cat("NEXT: Check cluster author subtype composition, high-donor and cross-lineage evidence;\n")
  cat("      then freeze macrophage/fibroblast/CD8 inclusion masks BEFORE validation UCell scoring.\n")
})
