# CRC-ImmuneMap | Step 11B | Targeted CELL-LEVEL annotation evidence
# Read-only review. NO approval, NO masks, NO UCell re-scoring.
# Process 3 Seurat checkpoints sequentially to limit memory on Windows.
# Does not use SPP1/FAP/checkpoint program genes to DEFINE lineages.
local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  if (!dir.exists(project_dir)) stop("Project not found: ", project_dir)
  setwd(project_dir)
  pkgs <- c("Seurat", "Matrix", "dplyr")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
  suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
  in_review <- "results/tables/final_annotation_review/step11_adjudication/all_36_clusters_REVIEW_AND_SIGNOFF.csv"
  if (!file.exists(in_review)) stop("Missing Step 11A worksheet: ", in_review)
  worksheet <- read.csv(in_review, stringsAsFactors = FALSE, check.names = FALSE)
  expected <- c("lineage", "validation_cluster", "n_cells", "proposal_NOT_APPROVED", "reviewer_decision")
  if (!all(expected %in% names(worksheet))) stop("Worksheet columns missing: ", paste(setdiff(expected, names(worksheet)), collapse=", "))
  if (nrow(worksheet) != 36L) stop("Expected exactly 36 cluster rows in Step 11A worksheet")
  if (anyDuplicated(paste(worksheet$lineage, worksheet$validation_cluster))) stop("Repeated cluster ID in worksheet")
  out_dir <- "results/tables/final_annotation_review/step11b_targeted_evidence"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- list(
    myeloids = list(
      path="objects/crc_validation_07c_myeloids_clustered_unreviewed.rds",
      focus=c("0","3","7"), controls=c("1","2","5"),
      panels=list(
        primary=c("C1QA","C1QC","CD68"), myeloid=c("LYZ","LST1","TYROBP"),
        dc=c("FCER1A","CD1C","CLEC10A"), mono=c("FCN1","S100A8","S100A9"),
        epithelial=c("EPCAM","KRT8","KRT18"),
        tcell=c("CD3D","CD3E","TRAC")) ,
      cross="epithelial"
    ),
    stromal = list(
      path="objects/crc_validation_07c_stromal_clustered_unreviewed.rds",
      focus=c("0","12"), controls=c("2","3","4","9","11","10"),
      panels=list(
        primary=c("COL1A1","DCN","LUM"), epithelial=c("EPCAM","KRT8","KRT18"),
        endothelial=c("PECAM1","VWF","CLDN5"), mural=c("RGS5","PDGFRB","MYH11"),
        myeloid=c("LYZ","LST1","TYROBP"), tcell=c("CD3D","CD3E","TRAC")),
      cross="epithelial"
    ),
    tcells = list(
      path="objects/crc_validation_07c_tcells_clustered_unreviewed.rds",
      focus=c("2","4","5","10"), controls=c("1","8"),
      panels=list(
        primary=c("CD3D","CD3E","TRAC"), cd8=c("CD8A","CD8B"),
        gamma_delta=c("TRDC","TRGC1","TRGC2"),
        nk=c("KLRD1","FCGR3A","NCR1"), cd4=c("CD4","IL7R","CCR7"),
        treg=c("FOXP3","IL2RA"), epithelial=c("EPCAM","KRT8","KRT18"),
        myeloid=c("LYZ","LST1","TYROBP")),
      cross="gamma_delta"
    )
  )
  all_summaries <- list()
  all_audits <- list()
  safe_med <- function(x) if (all(is.na(x)) || !length(x)) NA_real_ else median(x, na.rm=TRUE)
  safe_p90 <- function(x) if (all(is.na(x)) || !length(x)) NA_real_ else as.numeric(quantile(x, 0.9, na.rm=TRUE, names=FALSE))

  for (key in names(cfg)) {
    conf <- cfg[[key]]
    cat("\n=== 11B: loading ", key, " RDS ===\n", sep="")
    if (!file.exists(conf$path)) stop("RDS missing: ", conf$path)
    obj <- readRDS(conf$path)
    if (!inherits(obj,"Seurat")) stop("Checkpoint is not a Seurat object: ", key)
    DefaultAssay(obj) <- "RNA"
    md <- obj[[]]
    require_md <- c("validation_cluster","patient_id","sample_id","tissue_region","authors_cell_subtype","nCount_RNA","nFeature_RNA")
    if (!all(require_md %in% names(md))) stop("Missing metadata in ", key, ": ", paste(setdiff(require_md,names(md)),collapse=", "))
    if (!identical(rownames(md),colnames(obj))) stop("Metadata cell order mismatch: ", key)
    md$cell_id <- rownames(md)
    md$validation_cluster <- as.character(md$validation_cluster)
    md$authors_cell_subtype <- as.character(md$authors_cell_subtype)
    md$patient_id <- as.character(md$patient_id)
    md$sample_id <- as.character(md$sample_id)
    md$tissue_region <- as.character(md$tissue_region)
    if (anyNA(md[,c("cell_id","validation_cluster","patient_id","sample_id","tissue_region")])) stop("Missing essential metadata: ",key)
    if (anyDuplicated(md$cell_id)) stop("Duplicate cell IDs in ",key)

    ws <- worksheet[worksheet$lineage==key,,drop=FALSE]
    actual_tab <- table(md$validation_cluster)
    exp_counts <- setNames(as.integer(ws$n_cells), as.character(ws$validation_cluster))
    if (!setequal(names(actual_tab), names(exp_counts)) || !identical(as.integer(actual_tab[sort(names(exp_counts))]),as.integer(exp_counts[sort(names(exp_counts))]))) {
      # Explicit name matching independent of lexical/numeric ordering.
      if (!setequal(names(actual_tab),names(exp_counts)) || any(as.integer(actual_tab[names(exp_counts)]) != unname(exp_counts)))
        stop("Cluster membership/counts changed since 11A in ",key,"; STOP, do not transfer decisions")
    }

    selected <- md$validation_cluster %in% unique(c(conf$focus, conf$controls))
    # Keep all focus & comparison cells; avoid whole-object diagnostics where not needed.
    ids <- md$cell_id[selected]
    mini <- md[selected, c("cell_id","validation_cluster","patient_id","sample_id","tissue_region","authors_cell_subtype","nCount_RNA","nFeature_RNA"),drop=FALSE]
    mini$review_role <- ifelse(mini$validation_cluster %in% conf$focus,"PRIORITY_REVIEW","COMPARATOR")
    raw <- Seurat::GetAssayData(obj, assay="RNA", layer="counts")
    if (!identical(colnames(raw), md$cell_id)) stop("Counts order mismatch in ",key)
    use <- unique(unlist(conf$panels, use.names=FALSE))
    absent <- setdiff(use, rownames(raw))
    if (length(absent)) cat("Markers absent and NOT imputed: ",paste(absent,collapse=", "),"\n")
    for (panel in names(conf$panels)) {
      genes <- intersect(conf$panels[[panel]],rownames(raw))
      if (length(genes)<1L) stop("Panel completely absent: ",key,"/",panel)
      part <- raw[genes, ids, drop=FALSE]
      mini[[paste0("detect_",panel)]] <- as.integer(Matrix::colSums(part>0))
      mini[[paste0("umi_",panel)]] <- as.numeric(Matrix::colSums(part))
    }
    mini$primary_and_cross_ge2 <- mini$detect_primary >= 2 & mini[[paste0("detect_",conf$cross)]] >= 2
    mini$primary_lineage_ge2 <- mini$detect_primary >= 2
    mini$cross_panel_umi_fraction <- ifelse(mini$nCount_RNA>0,
          mini[[paste0("umi_",conf$cross)]] / mini$nCount_RNA, NA_real_)
    mini$primary_panel_umi_fraction <- ifelse(mini$nCount_RNA>0,mini$umi_primary/mini$nCount_RNA,NA_real_)
    mini$published_cd8 <- if (key=="tcells") mini$authors_cell_subtype=="CD8+ T cells" else NA

    group_cols <- c("validation_cluster","review_role")
    summaries <- mini |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        n_cells=dplyr::n(),
        n_patients=dplyr::n_distinct(patient_id),
        n_samples=dplyr::n_distinct(sample_id),
        n_tumor_core=sum(tissue_region=="Tumor_Core"),
        n_tumor_border=sum(tissue_region=="Tumor_Border"),
        n_normal=sum(tissue_region=="Normal"),
        pct_primary_ge2=round(100*mean(primary_lineage_ge2),1),
        pct_primary_and_cross_ge2=round(100*mean(primary_and_cross_ge2),1),
        pct_cross_any=round(100*mean(.data[[paste0("detect_",conf$cross)]]>0),1),
        median_total_umi=safe_med(nCount_RNA),
        median_nfeatures=safe_med(nFeature_RNA),
        median_cross_panel_umi=safe_med(.data[[paste0("umi_",conf$cross)]]),
        p90_cross_panel_umi=safe_p90(.data[[paste0("umi_",conf$cross)]]),
        median_cross_fraction=safe_med(cross_panel_umi_fraction),
        p90_cross_fraction=safe_p90(cross_panel_umi_fraction),
        median_total_umi_cross_positive=safe_med(nCount_RNA[primary_and_cross_ge2]),
        median_total_umi_not_cross_positive=safe_med(nCount_RNA[!primary_and_cross_ge2]),
        pct_authors_cd8=if(key=="tcells") round(100*mean(published_cd8,na.rm=TRUE),1) else NA_real_,
        pct_cd8_detect_any=if(key=="tcells") round(100*mean(detect_cd8>=1),1) else NA_real_,
        pct_gamma_detect_any=if(key=="tcells") round(100*mean(detect_gamma_delta>=1),1) else NA_real_,
        .groups="drop"
      ) |> dplyr::arrange(as.integer(validation_cluster))
    # Avoid using group-wise mutable n_cells alias when calculating donor proportions.
    patient_ct <- mini |>
      dplyr::count(validation_cluster,patient_id,name="patient_n") |>
      dplyr::group_by(validation_cluster) |>
      dplyr::mutate(patient_fraction=patient_n/sum(patient_n)) |>
      dplyr::ungroup()
    high <- patient_ct |>
      dplyr::group_by(validation_cluster) |>
      dplyr::slice_max(patient_n,n=1,with_ties=FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(validation_cluster,dominant_patient=patient_id,
                       dominant_patient_n=patient_n,
                       dominant_patient_pct=round(100*patient_fraction,1))
    summaries <- dplyr::left_join(summaries,high,by="validation_cluster")
    summaries$lineage <- key
    summaries$cross_panel_interpretation <- if (key=="tcells")
       "Gamma-delta transcript detection is not definitive lineage exclusion" else
       "Epithelial co-detection is not proof of doublets"
    if (any(summaries$n_cells != unname(exp_counts[as.character(summaries$validation_cluster)]))) stop("Counts unexpectedly mismatch",key)
    write.csv(summaries,file.path(out_dir,paste0(key,"_focus_and_comparators_cluster_summary.csv")),row.names=FALSE,na="",fileEncoding="UTF-8")

    # Focus cell diagnostics: supports deeper investigation, not auto cell exclusion.
    focus <- mini[mini$review_role=="PRIORITY_REVIEW",,drop=FALSE]
    write.csv(focus,file.path(out_dir,paste0(key,"_priority_per_cell_diagnostics_NOT_MASK.csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
    focused_subtypes <- focus |>
      dplyr::count(validation_cluster,patient_id,sample_id,tissue_region,authors_cell_subtype,name="n_cells") |>
      dplyr::arrange(as.integer(validation_cluster),patient_id,sample_id,authors_cell_subtype)
    write.csv(focused_subtypes,file.path(out_dir,paste0(key,"_priority_patient_sample_subtype.csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
    review_by_patient <- focus |>
      dplyr::group_by(validation_cluster,patient_id,tissue_region) |>
      dplyr::summarise(n_cells=dplyr::n(),
        pct_primary_ge2=round(100*mean(primary_lineage_ge2),1),
        pct_primary_and_cross_ge2=round(100*mean(primary_and_cross_ge2),1),
        median_cross_umi=median(.data[[paste0("umi_",conf$cross)]],na.rm=TRUE),
        median_cross_frac=safe_med(cross_panel_umi_fraction),
        pct_authors_cd8=if(key=="tcells") round(100*mean(published_cd8,na.rm=TRUE),1) else NA_real_,
        .groups="drop") |>
      dplyr::arrange(as.integer(validation_cluster),patient_id,tissue_region)
    write.csv(review_by_patient,file.path(out_dir,paste0(key,"_priority_patient_region_diagnostics.csv")),row.names=FALSE,na="",fileEncoding="UTF-8")
    all_summaries[[key]] <- summaries
    all_audits[[key]] <- data.frame(lineage=key, n_priority_clusters=length(conf$focus),
       n_priority_cells=nrow(focus),n_control_clusters=length(conf$controls),
       n_assessed_cells=nrow(mini),missing_genes=paste(absent,collapse="; "))
    cat("Audit: ",key," | priority cells ",nrow(focus)," | comparator cells ",nrow(mini)-nrow(focus),"\n",sep="")
    print(as.data.frame(summaries[summaries$review_role=="PRIORITY_REVIEW",
      c("validation_cluster","n_cells","n_patients","dominant_patient_pct",
        "pct_primary_ge2","pct_primary_and_cross_ge2","median_cross_panel_umi",
        "median_total_umi","pct_authors_cd8","pct_gamma_detect_any")]),row.names=FALSE)
    rm(obj,md,raw,mini,focus,part)
    invisible(gc(verbose=FALSE))
  }
  combined <- dplyr::bind_rows(all_summaries)
  write.csv(combined,file.path(out_dir,"all_focus_and_comparator_cluster_summary.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
  write.csv(dplyr::bind_rows(all_audits),file.path(out_dir,"step11b_audit_inventory.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
  readme <- c(
    "CRC-ImmuneMap | Step 11B targeted cell-level annotation evidence",
    "", "STATUS: DIAGNOSTIC ONLY — not reviewer sign-off or final mask.",
    "Exactly nine prioritized clusters reviewed alongside lineage-appropriate comparators.",
    "Cross-lineage 2-marker co-detection is an exploratory diagnostic, not proof of doublets.",
    "Interpret marker UMI fractions and total nCount/nFeature distributions together with original author subtype labels.",
    "T-cell gamma/delta panel is not a perfect exclusivity marker; CD8 and gamma/delta labels can overlap biologically.",
    "No lineage was selected on the basis of SPP1, FAP, PDCD1 or other program scores.",
    "No masking, annotation, score recalculation, correlation test or RDS modification was performed.",
    "After manual scientific adjudication: edit Step11A REVIEW_AND_SIGNOFF worksheet with rationales/initials/date.",
    "ONLY THEN implement final per-cell masks and recompute any changed summaries; do not treat cluster INCLUDE as all cells INCLUDE.",
    "No triple-program correlation is supportable with 2–3 complete-case validation tumor patients."
  )
  writeLines(readme,file.path(out_dir,"README_11B_INTERPRETATION.txt"),useBytes=TRUE)
  cat("\n========== STEP 11B TARGETED CELL LEVEL REVIEW COMPLETE ==========\n")
  cat("Priority clusters: ",sum(vapply(cfg,function(x) length(x$focus),integer(1)))," / 36\n",sep="")
  cat("Evidence tables saved to: ",normalizePath(out_dir,winslash="/"),"\n",sep="")
  cat("No decisions approved, no final masks, no UCell/re-analysis.\n")
})
