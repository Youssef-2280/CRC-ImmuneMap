# CRC-ImmuneMap — R / Seurat single-cell QC and annotation figures
# ALL figures are plotted from existing Seurat RDS checkpoints (NO AI imagery).
# QC figures are explicitly POST-FILTER/PROCESSED, unless copied from the original
# before-filtering QC plots (which are marked as originals, not reconstructed).
# NO re-clustering, re-normalization, differential expression, or new UCell analysis.
# User-approved primary_07e is an ANALYSIS selection, not experimentally verified purity.

local({
  project_dir <- Sys.getenv(
    "CRC_PROJECT_DIR",
    unset = "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  )
  if (!dir.exists(project_dir)) stop("Project folder not found: ", project_dir)
  setwd(project_dir)
  for (p in c("Seurat", "ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Please install R package: ", p)
  }
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ggplot2))
  set.seed(100)

  out_dir <- file.path("results", "figures", "R_singlecell_QC_annotation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  k <- 0L
  add_entry <- function(fig, source, note, status) {
    k <<- k + 1L
    manifest[[k]] <<- data.frame(
      figure = fig, input_checkpoint = source, note = note, status = status,
      stringsAsFactors = FALSE
    )
  }

  save_plot <- function(p, basename, source, note, width = 11, height = 7) {
    img <- file.path(out_dir, paste0(basename, ".png"))
    vec <- file.path(out_dir, paste0(basename, ".pdf"))
    print(p)
    ggplot2::ggsave(img, plot = p, width = width, height = height,
                    units = "in", dpi = 300, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(vec, plot = p, width = width, height = height,
                    units = "in", bg = "white", limitsize = FALSE)
    add_entry(basename, source, note, "GENERATED_FROM_RDS")
    cat("SAVED:", img, "\n")
  }
  try_figure <- function(expr, fig_name, source, note) {
    tryCatch(force(expr), error = function(e) {
      msg <- paste0("FAILED: ", conditionMessage(e))
      warning(fig_name, " -> ", msg, call. = FALSE)
      add_entry(fig_name, source, note, msg)
    })
  }

  has_metadata <- function(obj, col) col %in% colnames(obj[[]])
  pick_reduction <- function(obj, requested) {
    options <- c(requested, "umap_before", "umap_validation", "umap")
    found <- options[options %in% names(obj@reductions)]
    if (!length(found)) return(NULL)
    found[1]
  }
  actual_umap <- function(obj, preferred, group, title,
                          source, basename, label = FALSE, width = 11, height = 7) {
    if (!has_metadata(obj, group)) {
      warning("Skipping ", basename, ": missing metadata column ", group)
      add_entry(basename, source, "Metadata column missing", "SKIPPED")
      return(invisible(NULL))
    }
    reduction <- pick_reduction(obj, preferred)
    if (is.null(reduction)) {
      warning("Skipping ", basename, ": NO stored UMAP embeddings")
      add_entry(basename, source, "No UMAP embeddings", "SKIPPED")
      return(invisible(NULL))
    }
    p <- Seurat::DimPlot(obj, reduction = reduction, group.by = group,
                         pt.size = if (ncol(obj) > 30000) 0.07 else 0.17,
                         shuffle = TRUE, seed = 100, raster = ncol(obj) > 30000,
                         label = label, repel = label) +
      ggplot2::labs(title = title,
                    subtitle = paste0("Genuine stored ", reduction,
                                      " embeddings; n=", format(ncol(obj), big.mark=","), " cells")) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(plot.title = element_text(face = "bold"))
    save_plot(p, basename, source, paste("Precomputed", reduction, "embedding"), width, height)
  }

  qc_plot <- function(obj, group_col, title, source, basename, max_groups = 25) {
    md <- obj[[]]
    required <- c("nCount_RNA", "nFeature_RNA", "percent.mt", group_col)
    missing <- setdiff(required, colnames(md))
    if (length(missing)) {
      warning("Skipping QC ", basename, ": missing ", paste(missing, collapse=", "))
      add_entry(basename, source, "QC fields missing", "SKIPPED")
      return(invisible(NULL))
    }
    groups <- as.character(md[[group_col]])
    if (length(unique(groups)) > max_groups) {
      warning("Skipping QC ", basename, ": too many groups")
      add_entry(basename, source, "Too many groups", "SKIPPED")
      return(invisible(NULL))
    }
    long <- data.frame(
      group = rep(groups, 3),
      metric = factor(rep(c("Detected genes", "UMI counts", "Mitochondrial %"),
                          each = nrow(md)),
                      levels = c("Detected genes", "UMI counts", "Mitochondrial %")),
      value = c(as.numeric(md$nFeature_RNA), as.numeric(md$nCount_RNA),
                as.numeric(md$percent.mt))
    )
    long <- long[is.finite(long$value) & !is.na(long$group), , drop = FALSE]
    p <- ggplot(long, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width", adjust = 1.0, trim = TRUE, linewidth = 0.2) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth=0.2) +
      facet_wrap(~metric, scales = "free_y", nrow = 1) +
      labs(title = title, x = group_col, y = NULL,
           subtitle = "QC metrics of cells stored in the checkpoint; not pre-filtering") +
      theme_classic(base_size=11) +
      theme(legend.position="none", axis.text.x=element_text(angle=35,hjust=1),
            strip.background=element_rect(fill="grey95",color=NA),
            plot.title=element_text(face="bold"))
    save_plot(p, basename, source, "Stored-object QC; not prefilter", 14, 5.5)
  }

  dotplot <- function(obj, features, group_col, title, source, basename,
                      width = 12, height = 6) {
    if (!has_metadata(obj, group_col)) {
      add_entry(basename, source, "Group metadata missing", "SKIPPED")
      return(invisible(NULL))
    }
    present <- unique(features[features %in% rownames(obj)])
    if (length(present) < 4) {
      warning("Too few marker genes for ", basename)
      add_entry(basename, source, "Too few markers available", "SKIPPED")
      return(invisible(NULL))
    }
    DefaultAssay(obj) <- "RNA"
    p <- Seurat::DotPlot(obj, features = present, assay = "RNA",
                         group.by = group_col, dot.scale = 5) +
      Seurat::RotatedAxis() +
      ggplot2::labs(title = title,
                    subtitle = "Expression and detection in saved cell groups; descriptive marker evidence",
                    x = "Marker genes", y = group_col) +
      ggplot2::theme_bw(base_size=10) +
      ggplot2::theme(plot.title = element_text(face="bold"))
    save_plot(p, basename, source, "Seurat RNA normalized-data DotPlot", width, height)
  }

  open_if_exists <- function(path, label, action) {
    if (!file.exists(path)) {
      cat("SKIPPED missing checkpoint:", path, "\n")
      add_entry(label, path, "Checkpoint missing", "SKIPPED")
      return(invisible(NULL))
    }
    cat("\n=== Loading", label, "===\n", path, "\n")
    obj <- readRDS(path)
    stopifnot(inherits(obj, "Seurat"))
    action(obj)
    rm(obj)
    invisible(gc())
  }

  # Original BEFORE-filtering images are only copied when already produced by
  # the user's original discovery script. These are NOT generated or imputed.
  original_pre_qc <- c(
    "results/figures/qc_before_filtering_condition.png",
    "results/figures/qc_before_filtering_samples.png"
  )
  for (old in original_pre_qc) {
    name <- paste0("S00_original_", basename(old))
    if (file.exists(old)) {
      ok <- file.copy(old, file.path(out_dir, name), overwrite = TRUE)
      add_entry(name, old, "Copied original pre-filtering plot; no recomputation", if (ok) "COPIED_ORIGINAL" else "COPY_FAILED")
      cat("COPIED original prefilter QC:", old, "\n")
    } else {
      add_entry(name, old, "Original pre-filter QC not available; not fabricated", "NOT_AVAILABLE")
    }
  }

  # Discovery post-QC overview and reviewed broad annotation.
  discovery <- "objects/crc_discovery_annotated_v2.rds"
  open_if_exists(discovery, "Discovery broad annotated", function(obj) {
    if (has_metadata(obj, "condition")) {
      try_figure(qc_plot(obj, "condition", "GSE132465 | QC of retained cells",
                         discovery, "S01_Discovery_retained_cell_QC"), "S01", discovery, "")
    }
    try_figure(actual_umap(obj, "umap_before", "cell_type_manual",
                           "GSE132465 | Broad manual cell annotation",
                           discovery, "S02_Discovery_broad_annotation_UMAP", label = TRUE, width=13, height=8),
               "S02", discovery, "")
    if (has_metadata(obj, "condition")) {
      try_figure(actual_umap(obj, "umap_before", "condition",
                             "GSE132465 | Tumor and normal cells",
                             discovery, "S03_Discovery_tissue_UMAP", width=11, height=7),
                 "S03", discovery, "")
    }
    try_figure(dotplot(obj, c("EPCAM","KRT8","KRT18", "LYZ","LST1","C1QA",
                              "COL1A1","DCN","LUM","CD3D","CD8A","MS4A1",
                              "PECAM1","VWF","TPSAB1"),
                       "cell_type_manual", "GSE132465 | Broad lineage marker evidence",
                       discovery, "S04_Discovery_broad_marker_DotPlot", width=14, height=6),
               "S04", discovery, "")
  })

  # Discovery lineage annotations: use existing embeddings and saved labels.
  discovery_lineages <- list(
    list(path="objects/crc_discovery_myeloid_annotated_v1.rds", key="myeloid_cell_type", red="umap_myeloid", name="Myeloid"),
    list(path="objects/crc_discovery_fibroblast_annotated_v1.rds", key="fibroblast_state_preliminary", red="umap_fibroblast", name="Fibroblast"),
    list(path="objects/crc_discovery_tcells_annotated_v1.rds", key="tcell_cell_type", red="umap_tcells", name="Tcell")
  )
  for (j in seq_along(discovery_lineages)) {
    spec <- discovery_lineages[[j]]
    fname <- paste0("S", sprintf("%02d", 4L+j), "_Discovery_", spec$name,"_annotation_UMAP")
    open_if_exists(spec$path, paste("Discovery",spec$name), function(obj) {
      try_figure(actual_umap(obj, spec$red, spec$key,
                             paste("GSE132465 |",spec$name,"preliminary lineage annotation"),
                             spec$path, fname, label=TRUE, width=12, height=8),
                 fname, spec$path, "")
    })
  }

  # GSE144735 original AUTHORS-PROCESSED data. We did NOT apply discovery QC
  # thresholds in validation and don't depict it as an unfiltered experiment.
  validation <- "objects/crc_validation_gse144735_raw_imported.rds"
  open_if_exists(validation, "Validation authors-processed", function(obj) {
    try_figure(qc_plot(obj, "tissue_region", "GSE144735 | QC of authors-processed cells",
                       validation, "S08_Validation_authors_processed_QC"), "S08", validation, "")
    md <- obj[[]]
    if (all(c("authors_cell_type","tissue_region") %in% colnames(md))) {
      try_figure({
        df <- as.data.frame(table(md$authors_cell_type, md$tissue_region), stringsAsFactors=FALSE)
        names(df) <- c("cell_type","region","n_cells")
        df <- df[df$n_cells>0, , drop=FALSE]
        p <- ggplot(df, aes(x=region, y=n_cells, fill=cell_type)) +
          geom_col(position="stack", width=0.75) +
          labs(title="GSE144735 | Authors' broad cell types by tissue region",
               x="Tissue region", y="Cell count", fill="Published cell type") +
          theme_classic(base_size=11) + theme(plot.title=element_text(face="bold"))
        save_plot(p,"S09_Validation_broad_celltype_composition",validation,
                  "Published annotation and source metadata", 11,7)
      }, "S09", validation, "")
    }
  })

  validation_cfg <- list(
    myeloids = list(panel = c("LYZ","LST1","CD68","C1QA","C1QC","FCN1","CD1C","SPP1","EPCAM","KRT8"), primary=c("0","1","2","5"), cmp="Macrophage"),
    stromal = list(panel = c("COL1A1","DCN","LUM","PDGFRA","FAP","CTHRC1","ACTA2","PECAM1","RGS5","KRT8"), primary=c("0","2","3","4","9","11","12"), cmp="Fibroblast"),
    tcells = list(panel = c("CD3D","CD3E","TRAC","CD8A","CD8B","PDCD1","LAG3","TRDC","TRGC1","FOXP3"), primary=c("2","4","10"), cmp="CD8")
  )
  for (key in names(validation_cfg)) {
    sp <- validation_cfg[[key]]
    path <- file.path("objects", paste0("crc_validation_07c_",key,"_clustered_unreviewed.rds"))
    # These are authentic unsupervised clusters; the raw Seurat object itself
    # stores UNREVIEWED labels, separate from user-approved analytical masks.
    open_if_exists(path, paste("Validation",key), function(obj) {
      try_figure(actual_umap(obj, "umap_validation", "validation_cluster",
                             paste("GSE144735 |",key,"unsupervised subclusters"),
                             path, paste0("S10_Validation_",key,"_clusters_UMAP"),
                             label=TRUE, width=11,height=7),
                 paste0("S10_",key), path, "")
      try_figure(dotplot(obj, sp$panel, "validation_cluster",
                         paste("GSE144735 |",key,"lineage and state markers"),
                         path, paste0("S11_Validation_",key,"_marker_DotPlot"),
                         width=14,height=6), paste0("S11_",key), path, "")
      # User-approved analytical selection: per-cell CD8 author subtype gate.
      md <- obj[[]]
      if ("validation_cluster" %in% colnames(md)) {
        selected <- as.character(md$validation_cluster) %in% sp$primary
        if (key=="tcells") {
          if (!"authors_cell_subtype" %in% colnames(md)) {
            stop("Missing authors_cell_subtype; cannot visualize CD8 mask")
          }
          selected <- selected & as.character(md$authors_cell_subtype)=="CD8+ T cells"
        }
        expected <- c(myeloids=1745L, stromal=4230L, tcells=1264L)[key]
        if (sum(selected, na.rm=TRUE) != expected) {
          warning("Primary selection mismatch for ",key,": ",sum(selected),
                  " cells, expected ", expected, ". Skipping mask visualization.")
          add_entry(paste0("S12_Validation_",key,"_primary_07e_mask_UMAP"),path,
                    "Selection count differs from audited primary_07e", "SKIPPED_MISMATCH")
        } else {
          obj$primary07e_viz <- ifelse(selected,"primary_07e","Other cells")
          try_figure(actual_umap(obj,"umap_validation","primary07e_viz",
                                 paste("GSE144735 |",key,"primary_07e analytical selection"),
                                 path, paste0("S12_Validation_",key,"_primary_07e_mask_UMAP"),
                                 width=10,height=7), paste0("S12_",key), path, "")
          cat("VERIFIED primary_07e analytical selection for ",key,": ",sum(selected)," cells\n",sep="")
        }
      }
    })
  }

  # Emit manifest, including any skips and failures, and make provenance explicit.
  rows <- do.call(rbind, manifest)
  write.csv(rows, file.path(out_dir,"FIGURE_SOURCE_MANIFEST.csv"), row.names=FALSE)
  writeLines(c(
    "CRC-ImmuneMap: Single-cell R-generated QC and annotation figures",
    "",
    "All newly generated figures in this directory come from Seurat RDS checkpoints",
    "using ggplot2/Seurat. No image generation or invented embeddings are used.",
    "The existing 5 main results figures in the GitHub package are separate.",
    "",
    "S00: when available, copies original BEFORE-filtering discovery QC images.",
    "S01: discovery QC metrics on cells ALREADY RETAINED after filtering.",
    "S02-S07: original discovery UMAP reductions and cell annotations.",
    "S08: author-processed validation cells (NOT claimed as pre-filtering).",
    "S09: published cell type distribution by tissue region.",
    "S10-S11: validation lineage-specific clusters and markers; clusters are unreviewed.",
    "S12: user-approved primary_07e ANALYTICAL mask, with count verification.",
    "",
    "CAVEAT: these plots do not prove clean doublets, causal interactions,",
    "statistical significance or experimentally verified cell-type purity."
  ), file.path(out_dir,"README_FIGURES.txt"))
  cat("\n========== R SINGLE-CELL FIGURES COMPLETE ==========\n")
  cat("Input root:", normalizePath(project_dir, winslash="/"),"\n")
  cat("Output:", normalizePath(out_dir, winslash="/",mustWork=FALSE), "\n")
  cat("Saved figures:",sum(rows$status %in% c("GENERATED_FROM_RDS","COPIED_ORIGINAL")),"\n")
  cat("Missing / skipped / failed:",sum(!rows$status %in% c("GENERATED_FROM_RDS","COPIED_ORIGINAL")),"\n")
  print(rows[,c("figure","status")], row.names=FALSE)
})
