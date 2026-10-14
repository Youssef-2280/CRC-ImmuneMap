# CRC-ImmuneMap | Step 07C | GSE144735 broad-lineage subclustering and annotation REVIEW
# Starting point: Step 07B raw Seurat checkpoint. Does NOT repeat UMI import.
# IMPORTANT: This script does NOT assign definitive macrophage/CAF/CD8 identities,
# score signatures, or test correlations. Review cluster markers, donor contributions,
# cross-lineage signals, and author labels before selecting cells for validation.
# Author cell subtype labels are a REFERENCE, not ground-truth evidence on their own.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  stopifnot(dir.exists(project_dir))
  setwd(project_dir)

  pkgs <- c("Seurat", "Matrix", "dplyr", "ggplot2", "tibble", "tidyr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
  })
  set.seed(100)

  source_path <- "objects/crc_validation_gse144735_raw_imported.rds"
  if (!file.exists(source_path)) {
    stop("Step 07B Seurat checkpoint not found: ", source_path,
         "\nDo not rerun TSV conversion until checking the actual objects folder.")
  }
  table_dir <- "results/tables/validation/07c"
  figure_dir <- "results/figures/validation/07c"
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Marker panels include positive lineage markers and contaminating lineages.
  # Cross-lineage co-expression is a FLAG requiring review, not proof of doublets.
  cfg <- list(
    Myeloids = list(
      broad = "Myeloids", object_key = "myeloids", resolution = 0.4,
      panel = c("LYZ", "LST1", "TYROBP", "FCN1", "S100A8", "S100A9", "CD14",
                "C1QA", "C1QB", "C1QC", "CD68", "APOE", "SPP1", "TREM2",
                "GPNMB", "APOC1", "CCL18", "FOLR2", "LYVE1", "MARCO",
                "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CLEC4C",
                "LILRA4", "MKI67", "TOP2A", "CD3D", "EPCAM", "COL1A1")
    ),
    Stromal = list(
      broad = "Stromal cells", object_key = "stromal", resolution = 0.4,
      panel = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA",
                "FAP", "CTHRC1", "COL11A1", "POSTN", "INHBA", "TNC",
                "SERPINE1", "ACTA2", "TAGLN", "MYL9", "PI16", "CXCL12",
                "PECAM1", "VWF", "CLDN5", "RGS5", "PDGFRB", "MYH11",
                "SOX10", "PLP1", "S100B", "EPCAM", "KRT8", "CD3D", "LYZ")
    ),
    Tcells = list(
      broad = "T cells", object_key = "tcells", resolution = 0.4,
      panel = c("CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B",
                "GZMK", "NKG7", "GNLY", "PRF1", "GZMB", "CCL5", "XCL1",
                "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13",
                "LAYN", "ENTPD1", "FOXP3", "IL2RA", "CCR7", "SELL", "IL7R",
                "TRDC", "TRGC1", "KLRD1", "FCGR3A", "MKI67", "TOP2A",
                "LYZ", "EPCAM", "COL1A1")
    )
  )
  cat("=== Step 07C: loading saved validation object ===\n")
  validation <- readRDS(source_path)
  stopifnot(inherits(validation, "Seurat"))
  DefaultAssay(validation) <- "RNA"
  m <- validation[[]]
  required <- c("patient_id", "sample_id", "tissue_region",
                "authors_cell_type", "authors_cell_subtype",
                "nFeature_RNA", "nCount_RNA", "percent.mt")
  if (!all(required %in% names(m))) {
    stop("Missing metadata: ", paste(setdiff(required, names(m)), collapse = ", "))
  }
  stopifnot(identical(colnames(validation), rownames(m)))
  cat("Validation:", nrow(validation), "genes;", ncol(validation), "cells;",
      dplyr::n_distinct(m$patient_id), "patients.\n")

  audit_summary <- list()
  for (nm in names(cfg)) {
    info <- cfg[[nm]]
    key <- info$object_key
    cat("\n=== 07C:", nm, "| broad label:", info$broad, "===\n")
    selected_ids <- rownames(m)[!is.na(m$authors_cell_type) &
                                  m$authors_cell_type == info$broad]
    if (length(selected_ids) < 100L) stop("Too few cells in ", nm, "; verify metadata")
    obj <- subset(validation, cells = selected_ids)
    DefaultAssay(obj) <- "RNA"
    cat("Cells:", ncol(obj), "; patients:", dplyr::n_distinct(obj$patient_id), "\n")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                                verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj),
                     block.size = 500, verbose = FALSE)
    n_pcs <- min(30L, length(VariableFeatures(obj)), ncol(obj) - 1L)
    if (n_pcs < 20L) stop("Insufficient PCA components for ", nm)
    obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = n_pcs,
                  seed.use = 100, verbose = FALSE)
    dims <- 1:20
    graph_names <- paste0("val_", key, c("_nn", "_snn"))
    obj <- FindNeighbors(obj, reduction = "pca", dims = dims,
                         graph.name = graph_names, verbose = FALSE)
    obj <- FindClusters(obj, graph.name = graph_names[2],
                        resolution = info$resolution, random.seed = 100,
                        verbose = FALSE)
    obj$validation_cluster <- as.character(obj$seurat_clusters)
    obj <- RunUMAP(obj, reduction = "pca", dims = dims,
                   reduction.name = "umap_validation", reduction.key = "ValUMAP_",
                   seed.use = 100, verbose = FALSE)
    obj@misc$validation_07c <- list(source = source_path, source_label = info$broad,
                                    dims = dims, resolution = info$resolution,
                                    seed = 100, annotation = "UNREVIEWED")

    md <- tibble::as_tibble(obj[[]])
    cell_report <- md |>
      dplyr::count(validation_cluster, patient_id, sample_id, tissue_region,
                   authors_cell_subtype, name = "n_cells") |>
      dplyr::arrange(as.integer(validation_cluster), patient_id, sample_id)
    write.csv(cell_report, file.path(table_dir, paste0(key, "_cluster_by_patient_sample_subtype.csv")),
              row.names = FALSE)
    patient_cluster <- md |>
      dplyr::count(validation_cluster, patient_id, name = "n_cells")
    cluster_summary <- patient_cluster |>
      dplyr::group_by(validation_cluster) |>
      dplyr::summarise(n_cells = sum(n_cells), n_patients = dplyr::n_distinct(patient_id),
                       dominant_patient = patient_id[which.max(n_cells)],
                       dominant_patient_pct = round(100 * max(n_cells) / sum(n_cells), 1),
                       .groups = "drop")
    condition_sum <- md |>
      dplyr::count(validation_cluster, tissue_region, name = "n_cells") |>
      tidyr::pivot_wider(names_from = tissue_region, values_from = n_cells,
                         values_fill = 0)
    cluster_summary <- dplyr::left_join(cluster_summary, condition_sum,
                                        by = "validation_cluster") |>
      dplyr::arrange(as.integer(validation_cluster))
    write.csv(cluster_summary, file.path(table_dir, paste0(key, "_cluster_patient_summary.csv")),
              row.names = FALSE)
    write.csv(md |> dplyr::count(validation_cluster, authors_cell_subtype,
                                name = "n_cells") |>
                dplyr::group_by(validation_cluster) |>
                dplyr::mutate(pct_in_cluster = round(100 * n_cells / sum(n_cells), 1)) |>
                dplyr::ungroup() |>
                dplyr::arrange(as.integer(validation_cluster), dplyr::desc(n_cells)),
              file.path(table_dir, paste0(key, "_cluster_author_subtypes.csv")),
              row.names = FALSE)

    present <- intersect(info$panel, rownames(obj))
    missing_genes <- setdiff(info$panel, rownames(obj))
    if (length(present) < 12L) stop("Unexpected marker coverage for ", nm)
    write.csv(data.frame(lineage = nm, gene = info$panel,
                         present = info$panel %in% rownames(obj)),
              file.path(table_dir, paste0(key, "_marker_coverage.csv")), row.names = FALSE)
    if (length(missing_genes)) cat("Missing markers:", paste(missing_genes, collapse = ", "), "\n")
    raw_markers <- GetAssayData(obj, assay = "RNA", layer = "counts")[present, , drop = FALSE]
    cluster_levels <- sort(unique(obj$validation_cluster))
    detection <- do.call(rbind, lapply(cluster_levels, function(cl) {
      idx <- which(obj$validation_cluster == cl)
      data.frame(lineage = nm, validation_cluster = cl, gene = present,
                 pct_detected = round(100 * Matrix::rowMeans(raw_markers[, idx, drop = FALSE] > 0), 1),
                 n_cells = length(idx))
    }))
    write.csv(detection, file.path(table_dir, paste0(key, "_marker_detection_by_cluster.csv")),
              row.names = FALSE)

    # Exploratory markers: cell-level tests for annotation, NOT donor-level DE.
    Idents(obj) <- "validation_cluster"
    markers <- FindAllMarkers(obj, assay = "RNA", slot = "data", only.pos = TRUE,
                              test.use = "wilcox", min.pct = 0.25,
                              logfc.threshold = 0.25, densify = FALSE,
                              random.seed = 100, verbose = FALSE)
    write.csv(markers,
              file.path(table_dir, paste0(key, "_markers_all_exploratory.csv")),
              row.names = FALSE)
    if (nrow(markers)) {
      top <- markers |>
        dplyr::filter(p_val_adj < 0.05) |>
        dplyr::group_by(cluster) |>
        dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) |>
        dplyr::ungroup()
      write.csv(top, file.path(table_dir, paste0(key, "_markers_top10.csv")),
                row.names = FALSE)
    }

    p_cl <- DimPlot(obj, reduction = "umap_validation", group.by = "validation_cluster",
                    label = TRUE, repel = TRUE, pt.size = 0.15) +
      ggtitle(paste("GSE144735", nm, "clusters — NOT annotated"))
    p_region <- DimPlot(obj, reduction = "umap_validation", group.by = "tissue_region",
                        pt.size = 0.15) + ggtitle(paste(nm, "— tissue region"))
    p_patient <- DimPlot(obj, reduction = "umap_validation", group.by = "patient_id",
                         pt.size = 0.15) + ggtitle(paste(nm, "— patients"))
    p_dot <- DotPlot(obj, features = present, assay = "RNA",
                     group.by = "validation_cluster", dot.scale = 5) +
      RotatedAxis() + ggtitle(paste(nm, "— lineage / state / contamination markers"))
    ggsave(file.path(figure_dir, paste0(key, "_clusters.png")), p_cl,
           width = 11, height = 8, dpi = 220)
    ggsave(file.path(figure_dir, paste0(key, "_tissue_region.png")), p_region,
           width = 11, height = 8, dpi = 220)
    ggsave(file.path(figure_dir, paste0(key, "_patients.png")), p_patient,
           width = 11, height = 8, dpi = 220)
    ggsave(file.path(figure_dir, paste0(key, "_markers_dotplot.png")), p_dot,
           width = 17, height = 8, dpi = 240)

    out_object <- file.path("objects", paste0("crc_validation_07c_", key, "_clustered_unreviewed.rds"))
    saveRDS(obj, out_object, compress = TRUE)
    cat("Saved:", out_object, "\n")
    cat("Clusters:", nrow(cluster_summary), "; dominant-donor >70%:",
        sum(cluster_summary$dominant_patient_pct > 70), "\n")
    audit_summary[[nm]] <- data.frame(lineage = nm, n_cells = ncol(obj),
                                      n_clusters = nrow(cluster_summary),
                                      n_patients = dplyr::n_distinct(md$patient_id),
                                      high_donor_clusters = sum(cluster_summary$dominant_patient_pct > 70))
    rm(obj, raw_markers, markers, md, cell_report, patient_cluster, cluster_summary,
       condition_sum, detection, p_cl, p_region, p_patient, p_dot)
    invisible(gc())
  }
  summary_tbl <- dplyr::bind_rows(audit_summary)
  write.csv(summary_tbl, file.path(table_dir, "validation_07c_summary.csv"), row.names = FALSE)
  cat("\n========== STEP 07C SUMMARY ==========\n")
  print(tibble::as_tibble(summary_tbl), n = Inf, width = Inf)
  cat("\nSaved tables:", normalizePath(table_dir, winslash = "/"), "\n")
  cat("Saved plots:", normalizePath(figure_dir, winslash = "/"), "\n")
  cat("IMPORTANT: NO final macrophage/CAF/CD8 labels, NO UCell scores, NO associations.\n")
  cat("Next: review DotPlots, marker tables, and donor representation; THEN select lineage-validated cells.\n")
})
