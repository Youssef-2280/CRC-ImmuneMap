# CRC-ImmuneMap | Step 11A: Annotation adjudication workbench
# PURPOSE: Create ONE concise review table for 36 validation clusters.
# NO automatic annotation approval, cell-mask changes, UCell or statistics.
# Uses existing Step 10 templates; does not modify them or original RDS files.
# REVIEWER_DECISION is deliberately PENDING_REVIEW until explicitly completed.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  if (!dir.exists(project_dir)) stop("Project directory does not exist: ", project_dir)
  setwd(project_dir)
  in_dir <- "results/tables/final_annotation_review"
  out_dir <- file.path(in_dir, "step11_adjudication")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  keys <- c("myeloids", "stromal", "tcells")
  # These are independent PROPOSALS, not a final assignment or human sign-off.
  # The source is the earlier Step 10 evidence tables and review proposals.
  propose <- list(
    myeloids = c(
      "0"="NEEDS_CELL_LEVEL_REVIEW", "1"="CONDITIONAL_INCLUDE",
      "2"="CONDITIONAL_INCLUDE", "3"="NEEDS_CELL_LEVEL_REVIEW",
      "4"="PROPOSE_EXCLUDE", "5"="CONDITIONAL_INCLUDE",
      "6"="PROPOSE_EXCLUDE", "7"="NEEDS_CELL_LEVEL_REVIEW",
      "8"="SENSITIVITY_ONLY", "9"="SENSITIVITY_ONLY"
    ),
    stromal = c(
      "0"="NEEDS_CELL_LEVEL_REVIEW", "1"="PROPOSE_EXCLUDE",
      "2"="CONDITIONAL_INCLUDE", "3"="CONDITIONAL_INCLUDE",
      "4"="CONDITIONAL_INCLUDE", "5"="PROPOSE_EXCLUDE",
      "6"="PROPOSE_EXCLUDE", "7"="PROPOSE_EXCLUDE",
      "8"="PROPOSE_EXCLUDE", "9"="CONDITIONAL_INCLUDE",
      "10"="SENSITIVITY_ONLY", "11"="CONDITIONAL_INCLUDE",
      "12"="NEEDS_CELL_LEVEL_REVIEW", "13"="PROPOSE_EXCLUDE",
      "14"="PROPOSE_EXCLUDE"
    ),
    tcells = c(
      "0"="PROPOSE_EXCLUDE", "1"="PROPOSE_EXCLUDE",
      "2"="NEEDS_CELL_LEVEL_REVIEW", "3"="PROPOSE_EXCLUDE",
      "4"="NEEDS_CELL_LEVEL_REVIEW", "5"="NEEDS_CELL_LEVEL_REVIEW",
      "6"="PROPOSE_EXCLUDE", "7"="PROPOSE_EXCLUDE",
      "8"="PROPOSE_EXCLUDE", "9"="PROPOSE_EXCLUDE",
      "10"="NEEDS_CELL_LEVEL_REVIEW"
    )
  )

  # Evidence prompts, deliberately LINEAGE-BASED rather than selected on
  # SPP1/FAP/checkpoint program score itself.
  reasons <- list(
    myeloids = c(
      "0"="Macrophage and epithelial co-detection; KUL01 contributes 94.7%",
      "1"="Macrophage lineage; KUL21 contributes 77.4%; donor sensitivity needed",
      "2"="C1Q macrophage-associated markers; mostly Normal; six donors",
      "3"="DC-associated markers and mixed published subtypes",
      "4"="FCN1/monocyte-associated; not macrophage primary pool",
      "5"="C1Q macrophage lineage; proliferating annotations should be documented",
      "6"="FCN1/monocyte-associated; donor-dominant",
      "7"="Macrophage markers but mostly published Unknown; cell-level review",
      "8"="Rare and mixed myeloid; cross-lineage evidence",
      "9"="Rare cycling myeloid, only 24 cells"
    ),
    stromal = c(
      "0"="Myofibroblast markers; epithelial co-detection and donor effects",
      "1"="Endothelial-associated, exclude from fibroblast pool",
      "2"="Stromal 3/fibroblast markers; limited cross-lineage signal",
      "3"="Stromal 3/fibroblast markers; limited cross-lineage signal",
      "4"="Stromal 1/fibroblast markers; limited cross-lineage signal",
      "5"="Pericyte/mural-associated rather than fibroblast primary",
      "6"="Enteric glial-associated rather than fibroblast primary",
      "7"="Endothelial-associated rather than fibroblast primary",
      "8"="Smooth muscle/mural-associated rather than fibroblast primary",
      "9"="Fibroblast markers, some epithelial co-detection; review ambient signal",
      "10"="High epithelial UMI/co-detection and donor dominance",
      "11"="Fibroblast/stromal markers; mixed myofibroblast subtype",
      "12"="Fibroblast-associated, but high epithelial co-detection; no Normal cells",
      "13"="Smooth muscle/mural-associated rather than fibroblast primary",
      "14"="Lymphatic endothelial rather than fibroblast primary"
    ),
    tcells = c(
      "0"="Mostly CD4 T, not CD8 primary",
      "1"="Treg-like, not CD8 primary",
      "2"="723 of 855 author CD8; mixed cluster; apply PER-CELL CD8 gate",
      "3"="Mostly CD4 T, not CD8 primary",
      "4"="381 author CD8 plus 155 gamma-delta; PER-CELL review needed",
      "5"="Unknown/activated/Treg mix; donor-dominant, not CD8 wholesale",
      "6"="Mostly Th17, not CD8 primary",
      "7"="Mostly CD4 T, not CD8 primary",
      "8"="Mostly NK rather than CD8 T",
      "9"="Mostly Tfh rather than CD8 T",
      "10"="160 of 166 author CD8 but KUL19 = 92.8%; PER-CELL review"
    )
  )

  must_have <- c("validation_cluster", "n_cells", "n_primary_candidates",
                 "n_review_candidates", "n_patients", "dominant_patient_pct",
                 "top_author_subtypes", "screening_flags")
  parts <- lapply(keys, function(key) {
    path <- file.path(in_dir, paste0(key, "_decision_register_TEMPLATE.csv"))
    if (!file.exists(path)) stop("Missing Step 10 template: ", path)
    dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    miss <- setdiff(must_have, names(dat))
    if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse=", "))
    dat$validation_cluster <- as.character(dat$validation_cluster)
    if (anyDuplicated(dat$validation_cluster)) stop("Duplicate cluster in ", path)
    if (!setequal(dat$validation_cluster, names(propose[[key]]))) {
      stop("Cluster IDs differ from reviewed proposals for ", key,
           ". STOP: do not silently transfer decisions to changed clusters.")
    }
    ord <- match(names(propose[[key]]), dat$validation_cluster)
    dat <- dat[ord, , drop=FALSE]
    data.frame(
      lineage = key,
      validation_cluster = dat$validation_cluster,
      n_cells = dat$n_cells,
      n_primary_candidates = dat$n_primary_candidates,
      n_review_candidates = dat$n_review_candidates,
      n_patients = dat$n_patients,
      dominant_patient_pct = dat$dominant_patient_pct,
      top_author_subtypes = dat$top_author_subtypes,
      screening_flags = dat$screening_flags,
      proposal_NOT_APPROVED = unname(propose[[key]][dat$validation_cluster]),
      rationale_for_review = unname(reasons[[key]][dat$validation_cluster]),
      reviewer_decision = "PENDING_REVIEW",
      reviewer_rationale = "",
      reviewer_initials = "",
      review_date_YYYY_MM_DD = "",
      stringsAsFactors = FALSE
    )
  })
  proposed <- do.call(rbind, parts)
  stopifnot(nrow(proposed) == 36L, sum(proposed$n_primary_candidates)==7239L)

  work_file <- file.path(out_dir, "all_36_clusters_REVIEW_AND_SIGNOFF.csv")
  if (!file.exists(work_file)) {
    write.csv(proposed, work_file, row.names=FALSE, na="", fileEncoding="UTF-8")
    cat("CREATED review table (NOT APPROVED): ", work_file, "\n", sep="")
  } else {
    cat("PRESERVED your existing edited review table: ", work_file, "\n", sep="")
  }

  # A separate always-refreshable, read-only list of priority questions.
  priority <- proposed[proposed$proposal_NOT_APPROVED == "NEEDS_CELL_LEVEL_REVIEW",
                       c("lineage", "validation_cluster", "n_cells",
                         "n_primary_candidates", "dominant_patient_pct",
                         "rationale_for_review"), drop=FALSE]
  write.csv(priority, file.path(out_dir, "priority_9_clusters_for_cell_level_review.csv"),
            row.names=FALSE, fileEncoding="UTF-8")

  current <- read.csv(work_file, stringsAsFactors=FALSE, check.names=FALSE)
  if (nrow(current) != 36L || anyDuplicated(paste(current$lineage, current$validation_cluster)))
    stop("Edited review table has missing or duplicated cluster rows; no sign-off accepted")
  reference_ids <- paste(proposed$lineage, proposed$validation_cluster)
  if (!setequal(reference_ids, paste(current$lineage, current$validation_cluster)))
    stop("Edited review table has mismatched cluster IDs; no sign-off accepted")
  allowed <- c("PENDING_REVIEW", "INCLUDE", "EXCLUDE", "SENSITIVITY_ONLY", "NEEDS_FURTHER_REVIEW")
  if (!all(current$reviewer_decision %in% allowed))
    stop("reviewer_decision must be one of: ", paste(allowed, collapse=", "))
  if (any(current$reviewer_decision == "INCLUDE" & current$n_primary_candidates == 0L))
    stop("An INCLUDE decision has zero provisional primary candidate cells. Requires an expanded cell-level review; not supported by 07E mask")
  if (any(current$reviewer_decision != "PENDING_REVIEW" &
          (is.na(current$reviewer_rationale) | !nzchar(trimws(current$reviewer_rationale)) |
           is.na(current$reviewer_initials) | !nzchar(trimws(current$reviewer_initials)) |
           is.na(current$review_date_YYYY_MM_DD) |
           !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", current$review_date_YYYY_MM_DD))))
    stop("Every reviewed decision needs rationale, reviewer initials and YYYY-MM-DD date")

  tallies <- as.data.frame(table(lineage=current$lineage, decision=current$reviewer_decision),
                           stringsAsFactors=FALSE)
  tallies <- tallies[tallies$Freq > 0, , drop=FALSE]
  write.csv(tallies, file.path(out_dir, "review_decision_status.csv"),
            row.names=FALSE, fileEncoding="UTF-8")
  writeLines(c(
    "CRC-ImmuneMap | Step 11A annotation review workbench",
    "",
    "EDIT: all_36_clusters_REVIEW_AND_SIGNOFF.csv (do NOT edit original Step 10 templates).",
    "proposal_NOT_APPROVED is just a suggestion, not evidence of reviewer sign-off.",
    "Use reviewer_decision: INCLUDE, EXCLUDE, SENSITIVITY_ONLY or NEEDS_FURTHER_REVIEW.",
    "Also enter reviewer_rationale, reviewer_initials and review_date_YYYY_MM_DD for each reviewed row.",
    "Investigate the nine priority cell-level review clusters before deciding.",
    "A cluster-level INCLUDE is not sufficient to define CD8: the 07E per-cell author CD8 gate",
    "and independent lineage diagnostics still apply to clusters 2, 4 and 10.",
    "Donor dominance and epithelial co-detection are review flags, NOT automatic doublet proof.",
    "SPP1, FAP and checkpoint program scores must not determine which cells belong to a lineage.",
    "This script does NOT freeze cell masks or update any prior scores or results.",
    "After scientific review and explicit sign-off, run a separate Step 11B mask-freezing script.",
    "Validation only has 2 Core and 3 Core+Border complete-case tumor patients at the 20-cell cutoff.",
    "Do not use these data for confirmatory triple-program correlation statistics."
  ), file.path(out_dir, "HOW_TO_REVIEW.txt"), useBytes=TRUE)

  cat("\n========== STEP 11A: REVIEW WORKBENCH ==========\n")
  cat("Cluster proposals: ", nrow(proposed), " | Priority detailed reviews: ", nrow(priority), "\n", sep="")
  cat("Proposals (NOT final decisions):\n")
  print(as.data.frame(table(proposed$proposal_NOT_APPROVED)), row.names=FALSE)
  cat("Current reviewer decision status:\n")
  print(tallies, row.names=FALSE)
  cat("Nine priority clusters:\n")
  print(priority[, c("lineage", "validation_cluster", "rationale_for_review")], row.names=FALSE)
  cat("Review sheet: ", normalizePath(work_file, winslash="/"), "\n", sep="")
  cat("NOTE: NO masks were frozen and NO results recalculated.\n")
})
