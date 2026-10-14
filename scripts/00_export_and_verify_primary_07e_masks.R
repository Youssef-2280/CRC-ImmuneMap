# CRC-ImmuneMap — user-approved primary_07e masks: verify and freeze from existing R tables.
# Does not load large RDS, rerun clustering, compute UCell, or overwrite original 07E/07F.
local({
  project_dir <- Sys.getenv('CRC_PROJECT_DIR', unset='C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap')
  if (!dir.exists(project_dir)) stop('Project folder not found: ', project_dir)
  setwd(project_dir)
  flag_file <- 'results/tables/validation/07f_frozen_scoring/scored_cell_selection_flags.csv'
  score_file <- 'results/tables/validation/07f_frozen_scoring/all_scored_cells_with_flags.csv'
  if (!file.exists(flag_file)) stop('Missing Step07F per-cell selection flags: ', flag_file)
  flags <- read.csv(flag_file, stringsAsFactors=FALSE, check.names=FALSE)
  required <- c('cell_id','lineage','patient_id','sample_id','tissue_region',
                'validation_cluster','authors_cell_subtype','primary_mask','crosslineage_flag',
                'strict_mask','omit_risk_mask','all_published_cd8_mask')
  if (!all(required %in% names(flags))) stop('Unexpected 07F flags schema; do not continue')
  if (anyDuplicated(flags$cell_id)) stop('Repeated cell IDs across lineage tables')
  flags$primary_mask <- as.logical(flags$primary_mask)
  if (anyNA(flags$primary_mask)) stop('Missing/invalid primary mask')
  wanted <- c(Macrophages=1745L,Fibroblasts=4230L,CD8=1264L)
  counts <- table(factor(flags$lineage[flags$primary_mask],levels=names(wanted)))
  if (!identical(as.integer(counts),unname(wanted))) stop('Primary cell counts changed: ',paste(names(counts),counts,collapse='; '))
  expected_clusters <- list(Macrophages=c('0','1','2','5'),
                            Fibroblasts=c('0','2','3','4','9','11','12'),CD8=c('2','4','10'))
  lineages <- c(myeloids='Macrophages',stromal='Fibroblasts',tcells='CD8')
  for (key in names(lineages)) {
    label <- unname(lineages[[key]])
    f <- sprintf('results/tables/validation/07e_purity_coverage/%s_candidate_cells_NOT_FINAL.csv',key)
    if (!file.exists(f)) stop('Missing original 07E candidate table: ',f)
    cand <- read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
    if (!all(c('cell_id','candidate_tier','patient_id','validation_cluster') %in% names(cand))) stop('07E candidate schema incomplete: ',f)
    ids07e <- sort(cand$cell_id[cand$candidate_tier=='primary_candidate'])
    sub <- flags[flags$lineage==label & flags$primary_mask,,drop=FALSE]
    if (!identical(ids07e,sort(sub$cell_id))) stop('07E and 07F primary cell IDs disagree for ',key)
    if (!setequal(as.character(unique(sub$validation_cluster)),expected_clusters[[label]])) stop('Unexpected selected clusters in ',label)
    matchix <- match(sub$cell_id,cand$cell_id)
    if (anyNA(matchix) || any(sub$patient_id!=cand$patient_id[matchix]) ||
        any(as.character(sub$validation_cluster)!=as.character(cand$validation_cluster[matchix]))) {
      stop('Metadata disagreement between 07E and 07F for ',label)
    }
    if (label=='CD8') {
      if (anyNA(sub$authors_cell_subtype) || any(sub$authors_cell_subtype!='CD8+ T cells')) stop('CD8 author subtype per-cell gate failed')
      k <- table(factor(as.character(sub$validation_cluster),levels=c('2','4','10')))
      if (!identical(as.integer(k),c(723L,381L,160L))) stop('CD8 723/381/160 cluster counts changed')
    }
  }
  origmf <- read.csv('results/tables/programs/candidate_signature_manifest_v1.csv',stringsAsFactors=FALSE)
  valmf  <- read.csv('results/tables/validation/frozen_discovery_signature_manifest_v1.csv',stringsAsFactors=FALSE)
  reqmf <- c('signature_version','program','gene','method','max_rank')
  if (!all(reqmf %in% names(origmf)) || !all(reqmf %in% names(valmf)) ||
      !identical(origmf[reqmf],valmf[reqmf]) || nrow(origmf)!=29 ||
      anyDuplicated(origmf[c('program','gene')]) ||
      !all(origmf$method=='UCell' & origmf$max_rank==1500)) stop('Frozen signatures mismatch')
  if (!file.exists(score_file)) stop('Missing original 07F score table: ',score_file)
  sc <- read.csv(score_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!all(c('cell_id','patient_id','tissue_region','lineage','program','score','primary_mask') %in% names(sc))) stop('07F scores schema changed')
  if (!all(is.finite(sc$score))) stop('Nonfinite UCell scores')
  sc <- sc[as.logical(sc$primary_mask),,drop=FALSE]
  for (label in names(wanted)) {
    obs <- sort(unique(sc$cell_id[sc$lineage==label]))
    ids <- sort(flags$cell_id[flags$lineage==label & flags$primary_mask])
    if (!identical(obs,ids)) stop('Saved scores and primary mask disagree: ',label)
  }
  out <- 'results/tables/final_annotation_review/primary_07e_LOCKED'
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  primary <- flags[flags$primary_mask,required,drop=FALSE]
  primary <- primary[order(primary$lineage,primary$patient_id,primary$cell_id),,drop=FALSE]
  write.csv(primary,file.path(out,'primary_07e_FINAL_per_cell_mask.csv'),row.names=FALSE)
  write.csv(data.frame(lineage=names(wanted),n_final_primary=unname(wanted)),
            file.path(out,'primary_07e_FINAL_cell_counts.csv'),row.names=FALSE)
  m <- tools::md5sum(c(flag_file, score_file,
                     file.path(out,'primary_07e_FINAL_per_cell_mask.csv')))
  write.csv(data.frame(path=names(m),md5=unname(m)),
            file.path(out,'primary_07e_source_checksums.csv'),row.names=FALSE)
  writeLines(c('# PRIMARY_07E ANALYSIS SET LOCKED',
               'User approved primary_07e as final analysis selection on 2026-10-14.',
               'This does not certify biological identity or eliminate donor/cross-lineage uncertainty.',
               sprintf('Macrophages %d | Fibroblasts %d | CD8 %d',wanted[[1]],wanted[[2]],wanted[[3]]),
               'No UCell scores recomputed; original Step 07F scores were cross-checked.'),
             file.path(out,'PRIMARY_07E_LOCK_REPORT.txt'))
  cat('\n===== PRIMARY_07E LOCKED / VERIFIED =====\n')
  print(data.frame(lineage=names(wanted),cells=unname(wanted)),row.names=FALSE)
  cat('Written to: ',normalizePath(out,winslash='/'),'\n',sep='')
  cat('NOTE: Patient/sample source metadata and independent biological QC are separate checks.\n')
})
