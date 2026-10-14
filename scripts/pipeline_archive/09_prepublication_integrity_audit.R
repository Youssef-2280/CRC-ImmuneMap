# CRC-ImmuneMap — Step 09: Prepublication integrity & reproducibility audit
# Base R only. READ-ONLY with respect to prior analysis files.
# Outputs new files under results/tables/final_audit.

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
stopifnot(dir.exists(project_dir))
setwd(project_dir)
outdir <- "results/tables/final_audit"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  discovery_associations = "results/tables/integrated_evidence/discovery_primary_associations.csv",
  validation_paired = "results/tables/integrated_evidence/validation_primary_paired_ge20.csv",
  validation_sensitivity = "results/tables/integrated_evidence/validation_selection_sensitivity_combined_ge20.csv",
  evidence_report = "results/tables/integrated_evidence/integrated_evidence_summary.md",
  discovery_manifest = "results/tables/programs/candidate_signature_manifest_v1.csv",
  validation_manifest = "results/tables/validation/frozen_discovery_signature_manifest_v1.csv",
  validation_myeloids_rds = "objects/crc_validation_07c_myeloids_clustered_unreviewed.rds",
  validation_stromal_rds = "objects/crc_validation_07c_stromal_clustered_unreviewed.rds",
  validation_tcells_rds = "objects/crc_validation_07c_tcells_clustered_unreviewed.rds",
  discovery_macrophages_rds = "objects/crc_discovery_macrophages_selected.rds",
  discovery_fibroblasts_rds = "objects/crc_discovery_fibroblasts_selected.rds",
  discovery_cd8_rds = "objects/crc_discovery_cd8_selected.rds"
)

file_inventory <- data.frame(
  item = names(input_paths), path = unname(input_paths),
  exists = file.exists(input_paths), stringsAsFactors = FALSE
)
file_inventory$size_MB <- round(file.info(input_paths)$size / 1024^2, 2)
file_inventory$md5 <- NA_character_
small <- file_inventory$exists & !is.na(file_inventory$size_MB) & file_inventory$size_MB < 100
if (any(small)) file_inventory$md5[small] <- unname(tools::md5sum(file_inventory$path[small]))
write.csv(file_inventory, file.path(outdir, "input_file_inventory.csv"), row.names = FALSE)

checks <- data.frame(check = character(), status = character(), details = character())
add_check <- function(name, ok, details = "") {
  checks <<- rbind(checks, data.frame(check = name,
                    status = if (isTRUE(ok)) "PASS" else "FAIL",
                    details = as.character(details), stringsAsFactors = FALSE))
}
for (i in seq_len(nrow(file_inventory))) {
  add_check(paste0("File exists: ", file_inventory$item[i]), file_inventory$exists[i],
            file_inventory$path[i])
}

required <- c("discovery_associations", "validation_paired", "validation_sensitivity", "evidence_report")
if (!all(file.exists(input_paths[required]))) {
  write.csv(checks, file.path(outdir, "audit_checks.csv"), row.names = FALSE)
  stop("One or more required integration files are missing. See audit_checks.csv")
}

discovery <- read.csv(input_paths["discovery_associations"], check.names = FALSE)
paired <- read.csv(input_paths["validation_paired"], check.names = FALSE)
sensitivity <- read.csv(input_paths["validation_sensitivity"], check.names = FALSE)

cols_discovery <- c("comparison", "n_patients", "rho", "p_permutation", "p_BH_permutation")
cols_paired <- c("variant", "tumor_region", "program", "n_matched_patients",
                 "n_tumor_higher", "n_normal_higher", "n_equal", "median_paired_delta")
add_check("Discovery schema", all(cols_discovery %in% names(discovery)))
add_check("Validation paired schema", all(cols_paired %in% names(paired)))
add_check("Validation sensitivity schema", all(cols_paired %in% names(sensitivity)))
if (!all(vapply(tail(checks$status, 3), identical, logical(1), "PASS"))) {
  write.csv(checks, file.path(outdir, "audit_checks.csv"), row.names = FALSE)
  stop("A required table schema did not match. See audit_checks.csv")
}

add_check("Three discovery hypotheses", nrow(discovery) == 3)
add_check("Discovery patient n within cohort", all(discovery$n_patients > 0 & discovery$n_patients <= 23))
add_check("Discovery rho in [-1,1]", all(is.finite(discovery$rho) & abs(discovery$rho) <= 1))
add_check("Discovery adjusted p in [0,1]", all(!is.na(discovery$p_BH_permutation) &
          discovery$p_BH_permutation >= 0 & discovery$p_BH_permutation <= 1))
add_check("BH p not below raw p", all(discovery$p_BH_permutation + 1e-8 >= discovery$p_permutation))
add_check("Validation groups sum to matched n", all(paired$n_tumor_higher + paired$n_normal_higher + paired$n_equal == paired$n_matched_patients))
add_check("Validation sensitivity groups sum to matched n", all(sensitivity$n_tumor_higher + sensitivity$n_normal_higher + sensitivity$n_equal == sensitivity$n_matched_patients))
add_check("Validation matched n <= 6", all(paired$n_matched_patients >= 0 & paired$n_matched_patients <= 6))
add_check("No repeated program/region in primary summary", !anyDuplicated(paste(paired$tumor_region, paired$program, paired$variant)))
add_check("Core and Border treated as regions", all(paired$tumor_region %in% c("Tumor_Core", "Tumor_Border", "Tumor_Core_plus_Border")))
add_check("Four frozen program names", setequal(unique(paired$program),
           c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")))

# Count manifest entries without presuming source-specific column spellings.
manifest_counts <- data.frame(source = c("discovery", "validation"),
                              n_rows = c(NA_integer_, NA_integer_))
if (file.exists(input_paths["discovery_manifest"])) {
  manifest_counts$n_rows[1] <- nrow(read.csv(input_paths["discovery_manifest"], check.names = FALSE))
}
if (file.exists(input_paths["validation_manifest"])) {
  manifest_counts$n_rows[2] <- nrow(read.csv(input_paths["validation_manifest"], check.names = FALSE))
}
write.csv(manifest_counts, file.path(outdir, "manifest_row_counts.csv"), row.names = FALSE)
# Step 08 already validated the 29 unique program-gene entries and settings.
# Do NOT reinterpret manifest row counts as distinct program-gene counts.

write.csv(checks, file.path(outdir, "audit_checks.csv"), row.names = FALSE)
write.csv(discovery, file.path(outdir, "discovery_summary_for_figures.csv"), row.names = FALSE)
write.csv(paired, file.path(outdir, "validation_paired_summary_for_figures.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(outdir, "validation_sensitivity_summary_for_figures.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(outdir, "R_session_info.txt"))

primary <- subset(paired, variant == "primary_07e" & tumor_region == "Tumor_Core_plus_Border")
primary <- primary[match(c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity"), primary$program), ]
fmt <- function(x, digits = 3) format(round(x, digits), nsmall = digits, trim = TRUE)
report <- c(
  "# CRC-ImmuneMap — Step 09 integrity audit",
  "",
  sprintf("- Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- Technical checks passed: %d/%d", sum(checks$status == "PASS"), nrow(checks)),
  "- Discovery patient-level associations: see discovery_summary_for_figures.csv (n=11 for primary comparisons).",
  "- Validation: 6 patients in same published study, not an external-study cohort.",
  "- Validation Tumor Core+Border is one within-patient aggregate, not two independent observations.",
  "- High donor dominance/cross-lineage concerns remain; final lineage masks are NOT automatically approved by this audit.",
  "- Do not claim validated cross-program correlation or causal signaling.",
  "",
  "## Validation primary matched descriptive summary (>=20 cells each side)",
  "",
  "| Program | Covered patients | Tumor higher | Median tumor-normal UCell delta |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(primary))) {
  if (is.na(primary$program[i])) next
  report <- c(report, sprintf("| %s | %d | %d | %s |", primary$program[i],
    primary$n_matched_patients[i], primary$n_tumor_higher[i], fmt(primary$median_paired_delta[i])))
}
report <- c(report, "", "## Required human review before publishing", "",
  "- Confirm and document final cluster/cell inclusion masks for macrophages, fibroblasts and CD8.",
  "- Explain exclusion/review of epithelial cross-lineage, donor-dominated and rare clusters.",
  "- Assess whether all paired tissue labels genuinely represent matched samples for each patient.",
  "- State explicitly that validation patient coverage is insufficient for a meaningful triple-program correlation.",
  "- Add QC/doublet limitations and frozen signature provenance to the README.",
  "- Include software versions, code, small derived tables and figures; exclude raw data and large RDS.",
  "")
writeLines(report, file.path(outdir, "prepublication_integrity_report.md"))
cat("\n========== STEP 09 INTEGRITY AUDIT ==========\n")
print(checks, row.names = FALSE, right = FALSE)
cat("\nValidation primary Core+Border summary:\n")
print(primary[, c("program", "n_matched_patients", "n_tumor_higher", "median_paired_delta")], row.names = FALSE)
cat("\nAudit reports saved in: ", normalizePath(outdir, winslash = "/"), "\n", sep = "")
cat("Checks identify technical consistency only; annotation decisions require manual review.\n")
