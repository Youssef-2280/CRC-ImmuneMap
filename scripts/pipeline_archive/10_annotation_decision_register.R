# CRC-ImmuneMap — Step 10: FINAL ANNOTATION DECISION REGISTER (REVIEW ONLY)
# This does NOT assign final labels, change inclusion masks, run UCell, or alter RDS.
# Based on the evidence/purity CSVs from Steps 07D/07E and corrected donor counts.
# Requires base R only. Intended for use in RStudio on Windows.

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
stopifnot(dir.exists(project_dir))
setwd(project_dir)

base_evidence <- "results/tables/validation/07c/annotation_evidence_review"
base_purity <- "results/tables/validation/07e_purity_coverage"
out_dir <- "results/tables/final_annotation_review"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lineages <- c("myeloids", "stromal", "tcells")
files <- lapply(lineages, function(z) list(
  evidence = file.path(base_evidence, paste0(z, "_annotation_evidence_master.csv")),
  purity = file.path(base_purity, paste0(z, "_cluster_cell_level_purity_summary.csv"))
))
names(files) <- lineages
missing_inputs <- unlist(files, use.names = FALSE)
missing_inputs <- missing_inputs[!file.exists(missing_inputs)]
if (length(missing_inputs)) stop("Missing upstream files:\n", paste(missing_inputs, collapse="\n"))

read_table <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
fmt <- function(x) ifelse(is.na(x), "unknown", as.character(x))
collect_flags <- function(row, lineage) {
  flags <- character(0)
  if (isTRUE(row$dominant_patient_pct > 70)) flags <- c(flags, "donor_dominance_gt70pct")
  if (isTRUE(row$n_patients < 4)) flags <- c(flags, "represented_in_fewer_than_4_patients")
  if (identical(lineage,"myeloids") && isTRUE(row$pct_macrophage_and_epithelial_2genes > 20))
    flags <- c(flags, "myeloid_epithelial_codetection_gt20pct")
  if (identical(lineage,"stromal") && isTRUE(row$pct_fibroblast_and_epithelial_2genes > 20))
    flags <- c(flags, "fibroblast_epithelial_codetection_gt20pct")
  if (identical(lineage,"tcells") && isTRUE(row$pct_gammadelta_panel_any_gene > 20))
    flags <- c(flags, "gamma_delta_signal_gt20pct")
  if (length(flags) == 0) "none_of_screening_flags" else paste(flags, collapse = "; ")
}

summaries <- list()
for (lineage in lineages) {
  evidence <- read_table(files[[lineage]]$evidence)
  purity <- read_table(files[[lineage]]$purity)
  for (tbl in list(evidence, purity)) {
    stopifnot("validation_cluster" %in% names(tbl))
    stopifnot(!anyDuplicated(as.character(tbl$validation_cluster)))
  }
  evidence$validation_cluster <- as.character(evidence$validation_cluster)
  purity$validation_cluster <- as.character(purity$validation_cluster)
  stopifnot(setequal(evidence$validation_cluster, purity$validation_cluster))
  combined <- merge(evidence, purity, by = "validation_cluster", all = FALSE,
                    suffixes = c("", "_purity"), sort = FALSE)
  combined <- combined[order(as.integer(combined$validation_cluster)), , drop=FALSE]
  stopifnot(nrow(combined) == nrow(evidence))
  if ("n_cells_purity" %in% names(combined))
    stopifnot(all(combined$n_cells == combined$n_cells_purity))
  stopifnot(all(combined$n_primary_candidates + combined$n_review_candidates <= combined$n_cells))
  expected <- switch(lineage, myeloids = 1745L, stromal = 4230L, tcells = 1264L)
  stopifnot(sum(combined$n_primary_candidates) == expected)

  combined$lineage <- lineage
  combined$provisional_group <- ifelse(combined$n_primary_candidates > 0,
                                "primary_candidate",
                                ifelse(combined$n_review_candidates > 0,
                                       "review_candidate", "not_selected_provisionally"))
  combined$screening_flags <- vapply(seq_len(nrow(combined)), function(i)
    collect_flags(as.list(combined[i, , drop=FALSE]), lineage), character(1))
  # Decisions MUST be filled and signed off by a human reviewer, never assumed.
  combined$reviewer_decision <- "PENDING_REVIEW" # e.g. INCLUDE / EXCLUDE / SENSITIVITY_ONLY
  combined$reviewer_rationale <- ""
  combined$reviewer_initials <- ""
  combined$review_date_YYYY_MM_DD <- ""

  lead <- c("lineage", "validation_cluster", "provisional_group", "n_cells",
            "n_primary_candidates", "n_review_candidates", "n_patients",
            "dominant_patient", "dominant_patient_pct", "n_Normal",
            "n_Tumor_Core", "n_Tumor_Border", "top_author_subtypes",
            "screening_flags", "review_hint_not_final",
            "reviewer_decision", "reviewer_rationale", "reviewer_initials",
            "review_date_YYYY_MM_DD")
  combined <- combined[, c(intersect(lead, names(combined)),
                           setdiff(names(combined), lead)), drop=FALSE]
  output <- file.path(out_dir, paste0(lineage,"_decision_register_TEMPLATE.csv"))
  if (!file.exists(output)) {
    write.csv(combined, output, row.names=FALSE, na="")
    cat("CREATED:", output, "\n")
  } else {
    cat("PRESERVED EXISTING TEMPLATE (not overwritten):", output, "\n")
  }
  summaries[[lineage]] <- data.frame(
    lineage=lineage, clusters=nrow(combined), total_cells=sum(combined$n_cells),
    primary_candidate_cells=sum(combined$n_primary_candidates),
    review_candidate_cells=sum(combined$n_review_candidates),
    donor_dominant_clusters=sum(combined$dominant_patient_pct > 70, na.rm=TRUE),
    pending_review=nrow(combined), stringsAsFactors=FALSE
  )
}
summary_df <- do.call(rbind, summaries)
write.csv(summary_df, file.path(out_dir,"annotation_decision_register_summary.csv"), row.names=FALSE)

notes <- c(
  "# CRC-ImmuneMap — Annotation decision register (human review required)",
  "",
  "This is a REVIEW TEMPLATE, not a final approval, re-annotation or selection mask.",
  "Inputs: Step 07D evidence tables, Step 07E raw-count cell-level purity summaries.",
  "All input files are read-only; existing RDS and scores are left unchanged.",
  "",
  "## Human decision protocol",
  "1. Evaluate lineage markers and published subtypes independently of SPP1, FAP, or checkpoint program scores.",
  "2. For each cluster, inspect donor imbalance, sample-region coverage, cross-lineage UMI counts and plausible ambient RNA or doublets.",
  "3. Review mixed clusters at CELL LEVEL: cluster-level decisions cannot recreate partial CD8 inclusion masks (e.g. cluster 2/4/10).",
  "4. Write reviewer_decision INCLUDE / EXCLUDE / SENSITIVITY_ONLY / NEEDS_FURTHER_REVIEW, and an evidence-based rationale.",
  "5. Sign and date completed decisions; save edited versions under new filenames, preserving the templates.",
  "6. After review, generate final per-cell selection masks and RECOMPUTE scores/summary tables for any changed selections.",
  "7. Re-run integrity audit and verify paired-region sample identifiers from source clinical/sample metadata.",
  "",
  "## Flag meanings",
  "- donor_dominance_gt70pct is an exploratory review flag, NOT an exclusion criterion.",
  "- co-detection percentages alone do not prove doublets; consider RNA ambient signals and UMI burden.",
  "- FAP/SPP1/checkpoint gene expression must not be used as the sole criterion for defining the lineage pool.",
  "- Validation Core and Border represent regions from the same patient, not independent observations.",
  "- Only 2 Core and 3 Core+Border patients had three-program >=20-cell coverage in provisional analysis.",
  "- Neither patient-level triple-program association nor causal signaling was validated.",
  "",
  "No automatic final-annotation step is run by this script."
)
writeLines(notes, file.path(out_dir,"REVIEW_INSTRUCTIONS.md"), useBytes=TRUE)
cat("\n========== STEP 10: ANNOTATION DECISION REGISTER ==========\n")
print(summary_df, row.names=FALSE)
cat("Outputs:", normalizePath(out_dir, winslash="/"), "\n")
cat("STATUS: REVIEW TEMPLATES ONLY. No final labels, no RDS modification.\n")
