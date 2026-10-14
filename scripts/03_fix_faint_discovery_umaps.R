# CRC-ImmuneMap | Fix faint Discovery UMAP figures using REAL saved Seurat coordinates
# Works directly from crc_discovery_annotated_v2.rds; does not run QC, UMAP,
# clustering, annotation, or modify any analysis object.
# Fixes S02 (broad annotation) and S03 (Tumor / Normal) by plotting
# the genuine 63,689 stored cell coordinates with larger, saturated dots.

local({
  root <- Sys.getenv("CRC_PROJECT_DIR", unset = "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap")
  stopifnot(dir.exists(root))
  setwd(root)
  for (pkg in c("Seurat", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("R package missing: ", pkg)
  }
  library(Seurat)
  library(ggplot2)

  input <- "objects/crc_discovery_annotated_v2.rds"
  if (!file.exists(input)) stop("Cannot find RDS: ", normalizePath(getwd()), "/", input)
  dest <- "results/figures/R_singlecell_QC_annotation/HIGH_CONTRAST_REPLACEMENTS"
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  cat("Loading discovery annotated RDS (large object, please wait): ", input, "\n", sep="")
  obj <- readRDS(input)
  stopifnot(inherits(obj, "Seurat"))
  available <- names(obj@reductions)
  reduction <- if ("umap_before" %in% available) "umap_before" else if ("umap" %in% available) "umap" else stop("No saved UMAP: ", paste(available, collapse = ", "))
  coord <- as.data.frame(Seurat::Embeddings(obj[[reduction]])[, 1:2, drop=FALSE])
  names(coord) <- c("UMAP_1", "UMAP_2")
  md <- obj[[]]
  ix <- match(rownames(coord), rownames(md))
  if (anyNA(ix)) stop("Stored embedding cell IDs do not match Seurat metadata")
  for (col in c("cell_type_manual", "condition")) {
    if (!col %in% colnames(md)) stop("Missing column: ", col, "\nAvailable: ", paste(colnames(md), collapse=", "))
    coord[[col]] <- as.character(md[[col]][ix])
    coord[[col]][is.na(coord[[col]]) | !nzchar(coord[[col]])] <- "Unannotated"
  }
  n <- nrow(coord)
  cat("Stored embedding:", reduction, "| cells:", format(n, big.mark=","), "\n")
  cat("Condition values:", paste(unique(coord$condition), collapse=", "), "\n")

  clean_theme <- theme_classic(base_size=13) +
    theme(plot.title=element_text(face="bold", size=17),
          plot.subtitle=element_text(color="#516171", size=10),
          legend.title=element_text(face="bold"),
          axis.title=element_text(size=12),
          panel.grid=element_blank(),
          plot.margin=margin(10, 12, 10, 10))

  save_both <- function(p, base, width=11.8, height=7.8) {
    f_png <- file.path(dest, paste0(base, ".png"))
    f_pdf <- file.path(dest, paste0(base, ".pdf"))
    ggplot2::ggsave(f_png, p, width=width, height=height, dpi=360, bg="white", limitsize=FALSE)
    ggplot2::ggsave(f_pdf, p, width=width, height=height, bg="white", limitsize=FALSE)
    cat("SAVED:", f_png, "\n")
    invisible(NULL)
  }

  # S02: real broad annotation, high-contrast qualitative palette.
  ann <- sort(unique(coord$cell_type_manual))
  colors <- grDevices::hcl.colors(length(ann), palette="Dark 3")
  names(colors) <- ann
  p_ann <- ggplot(coord, aes(x=UMAP_1, y=UMAP_2, color=cell_type_manual)) +
    geom_point(size=0.27, alpha=0.82, shape=16, stroke=0) +
    scale_color_manual(values=colors, breaks=ann, drop=FALSE,
                       guide=guide_legend(override.aes=list(size=3.5,alpha=1))) +
    coord_equal() +
    labs(title="GSE132465 | Broad cell annotation",
         subtitle=paste0("Stored ", reduction, " embedding | ", format(n,big.mark=","), " cells | R-generated high-contrast rendering"),
         x="UMAP 1", y="UMAP 2", color="Cell type") +
    clean_theme + theme(legend.position="right")
  if (requireNamespace("ggrepel", quietly=TRUE) && length(ann) <= 20) {
    centers <- aggregate(coord[,c("UMAP_1","UMAP_2")],
                         by=list(cell_type_manual=coord$cell_type_manual),
                         FUN=median)
    p_ann <- p_ann + ggrepel::geom_label_repel(data=centers,
           aes(x=UMAP_1,y=UMAP_2,label=cell_type_manual), inherit.aes=FALSE,
           size=3, fill="white", alpha=0.85, color="#1f2937",
           min.segment.length=0, seed=100, max.overlaps=Inf, box.padding=.35)
  }
  save_both(p_ann, "S02_Discovery_broad_annotation_UMAP_HIGH_CONTRAST", width=13.5, height=8.3)

  # S03: actual Normal/Tumor. Overlay and separate panels for clarity.
  cond <- sort(unique(coord$condition))
  known <- c("Normal"="#E15241", "Tumor"="#087E95",
             "normal"="#E15241", "tumor"="#087E95")
  missing <- setdiff(cond, names(known))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette="Dark 3")
    names(extra) <- missing
    known <- c(known, extra)
  }
  tissue_cols <- known[cond]
  p_tissue <- ggplot(coord, aes(UMAP_1, UMAP_2, color=condition)) +
    geom_point(size=.34, alpha=.78, shape=16, stroke=0) +
    scale_color_manual(values=tissue_cols,
                       guide=guide_legend(override.aes=list(size=4,alpha=1))) +
    coord_equal() +
    labs(title="GSE132465 | Tumor and normal cells",
         subtitle=paste0("Stored ", reduction, " embedding | ", format(n,big.mark=","), " cells | high-contrast overlay"),
         x="UMAP 1", y="UMAP 2", color="Tissue") +
    clean_theme
  save_both(p_tissue, "S03_Discovery_tissue_UMAP_HIGH_CONTRAST", width=12, height=8)

  # Optional alternate view: a panel per tissue, maintaining UMAP axes.
  p_split <- ggplot(coord, aes(UMAP_1, UMAP_2, color=condition)) +
    geom_point(size=.38, alpha=.85, shape=16, stroke=0) +
    scale_color_manual(values=tissue_cols, guide="none") +
    facet_wrap(~condition, nrow=1) + coord_equal() +
    labs(title="GSE132465 | Tumor and normal UMAP (separate panels)",
         subtitle=paste0("Same stored ", reduction, " embedding; panels share coordinate limits"),
         x="UMAP 1", y="UMAP 2") + clean_theme +
    theme(strip.text=element_text(face="bold",size=12),
          strip.background=element_rect(fill="#f0f4f7",color=NA))
  save_both(p_split, "S03_Discovery_tissue_UMAP_SPLIT_HIGH_CONTRAST", width=14, height=6.7)

  rm(obj, md, coord)
  invisible(gc())
  cat("\nDONE. Three figures were exported as PNG + PDF.\n")
  cat("Folder: ", normalizePath(dest), "\n", sep="")
  cat("IMPORTANT: Changes are visualization only. No analysis results changed.\n")
})
