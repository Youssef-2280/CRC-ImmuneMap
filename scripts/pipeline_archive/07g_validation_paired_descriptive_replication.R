# CRC-ImmuneMap | Step 07G
# Validation: patient-matched, descriptive tumor-vs-normal replication
# This script READS saved Step 07F CSV tables only. It DOES NOT rerun UCell,
# clustering, download, or import. It does not change existing results.
# IMPORTANT: Descriptive analysis ONLY, with six patients and limited CD8 coverage.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  if (!dir.exists(project_dir)) stop("Project directory not found: ", project_dir)
  setwd(project_dir)

  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the packages 'dplyr' and 'ggplot2' first.")
  }
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  input_dir <- "results/tables/validation/07f_frozen_scoring"
  output_dir <- "results/tables/validation/07g_paired_descriptive"
  figures_dir <- "results/figures/validation/07g_paired_descriptive"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  file_region <- file.path(input_dir, "candidate_scores_by_patient_region.csv")
  file_combined <- file.path(input_dir, "candidate_scores_core_plus_border_per_patient.csv")
  file_coverage <- file.path(input_dir, "patient_three_program_coverage.csv")
  missing_files <- c(file_region, file_combined, file_coverage)[
    !file.exists(c(file_region, file_combined, file_coverage))
  ]
  if (length(missing_files)) stop("Missing Step 07F inputs: ", paste(missing_files, collapse = "; "))

  region <- read.csv(file_region, stringsAsFactors = FALSE, check.names = FALSE)
  combined <- read.csv(file_combined, stringsAsFactors = FALSE, check.names = FALSE)
  cover <- read.csv(file_coverage, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("variant", "patient_id", "tissue_region", "lineage", "program",
                "n_cells", "mean_score", "median_score")
  if (!all(required %in% names(region)) || !all(required %in% names(combined))) {
    stop("Unexpected Step 07F score-table schema")
  }
  if (!all(c("variant", "patient_id", "tissue_region", "all_three_ge20") %in% names(cover))) {
    stop("Unexpected Step 07F coverage-table schema")
  }
  programs <- c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")
  variants <- c("primary_07e", "crosslineage_sensitivity",
                "omit_high_risk_clusters", "all_published_cd8_sensitivity")
  tissue_regions <- c("Tumor_Core", "Tumor_Border", "Tumor_Core_plus_Border")
  if (!setequal(unique(region$variant), variants) ||
      !setequal(unique(region$program), programs) ||
      !setequal(unique(combined$variant), variants)) {
    stop("Unexpected variants/program names; stop to check the frozen analysis")
  }
  if (!all(combined$tissue_region == "Tumor_Core_plus_Border")) {
    stop("Combined-region labels are not as expected")
  }
  if (anyNA(region[c("variant", "patient_id", "tissue_region", "program", "n_cells", "mean_score")]) ||
      anyNA(combined[c("variant", "patient_id", "program", "n_cells", "mean_score")])) {
    stop("Unexpected missing values in scored rows")
  }
  if (any(region$n_cells <= 0 | combined$n_cells <= 0)) stop("Nonpositive cell counts")

  keycols <- c("variant", "patient_id", "tissue_region", "program")
  if (anyDuplicated(region[keycols]) || anyDuplicated(combined[keycols])) {
    stop("Repeated patient-region-program entries")
  }
  if (length(unique(region$patient_id)) != 6L) {
    warning("Expected 6 validation patients; found ", length(unique(region$patient_id)))
  }

  # "Core + Border" is ONE within-patient aggregate, not a new sample or patient.
  all_scores <- bind_rows(region[, required], combined[, required])
  normal <- all_scores |>
    filter(tissue_region == "Normal") |>
    select(variant, patient_id, program,
           n_normal = n_cells, mean_normal = mean_score, median_normal = median_score)
  tumor <- all_scores |>
    filter(tissue_region %in% tissue_regions) |>
    transmute(variant, patient_id, program, tumor_region = tissue_region,
              n_tumor = n_cells, mean_tumor = mean_score, median_tumor = median_score)

  paired_all <- inner_join(tumor, normal,
                           by = c("variant", "patient_id", "program")) |>
    mutate(delta_mean = mean_tumor - mean_normal,
           delta_median = median_tumor - median_normal)
  if (anyDuplicated(paired_all[c("variant", "patient_id", "tumor_region", "program")])) {
    stop("Unexpected duplicate paired observations")
  }
  if (any(!is.finite(paired_all$delta_mean))) stop("Non-finite patient score differences")

  # Thresholds are sensitivity descriptions, NOT p-value optimization.
  paired <- bind_rows(lapply(c(10L, 20L), function(k) {
    paired_all |>
      filter(n_normal >= k, n_tumor >= k) |>
      mutate(min_cells_per_region = k)
  })) |>
    arrange(variant, min_cells_per_region, tumor_region, program, patient_id)
  write.csv(paired_all,
            file.path(output_dir, "all_observed_matched_deltas_UNFILTERED.csv"),
            row.names = FALSE)
  write.csv(paired,
            file.path(output_dir, "matched_patient_deltas_ge10_ge20.csv"),
            row.names = FALSE)

  summary <- paired |>
    group_by(variant, min_cells_per_region, tumor_region, program) |>
    summarise(n_matched_patients = n(),
              n_tumor_higher = sum(delta_mean > 0),
              n_normal_higher = sum(delta_mean < 0),
              n_equal = sum(delta_mean == 0),
              median_paired_delta = median(delta_mean),
              min_paired_delta = min(delta_mean),
              max_paired_delta = max(delta_mean),
              .groups = "drop") |>
    arrange(variant, min_cells_per_region, tumor_region, program)
  write.csv(summary,
            file.path(output_dir, "descriptive_paired_summary_by_variant.csv"),
            row.names = FALSE)

  # Report complete-three-program coverage separately; no correlation p-values.
  coverage_summary <- cover |>
    group_by(variant, tissue_region) |>
    summarise(n_patients_ge20_all_three = sum(as.logical(all_three_ge20)),
              .groups = "drop") |>
    arrange(variant, tissue_region)
  write.csv(coverage_summary,
            file.path(output_dir, "three_compartment_patient_coverage.csv"),
            row.names = FALSE)

  # Explicitly compare sensitivity variants on the SAME matched patients
  # (a score change can otherwise just reflect a different eligible cohort).
  comparison <- paired |>
    filter(min_cells_per_region == 20L) |>
    select(variant, patient_id, tumor_region, program, delta_mean, n_tumor, n_normal) |>
    rename(variant_delta = delta_mean)
  primary_deltas <- comparison |>
    filter(variant == "primary_07e") |>
    transmute(patient_id, tumor_region, program, primary_delta = variant_delta)
  sensitivity_common <- comparison |>
    filter(variant != "primary_07e") |>
    inner_join(primary_deltas,
               by = c("patient_id", "tumor_region", "program")) |>
    mutate(delta_difference_from_primary = variant_delta - primary_delta,
           same_sign_as_primary = sign(variant_delta) == sign(primary_delta))
  write.csv(sensitivity_common,
            file.path(output_dir, "common_patient_selection_sensitivity_ge20.csv"),
            row.names = FALSE)

  primary_summary <- summary |>
    filter(variant == "primary_07e", min_cells_per_region == 20L)
  cat("\n========== STEP 07G — PAIRED DESCRIPTIVE VALIDATION ==========\n")
  cat("Frozen program scores read from Step 07F; NO Seurat objects loaded.\n")
  cat("6 patients total; thresholds count cells within each compared region.\n\n")
  print(as.data.frame(primary_summary), row.names = FALSE)
  cat("\nTHREE-COMPARTMENT COVERAGE:\n")
  print(as.data.frame(coverage_summary |>
                        filter(variant == "primary_07e")), row.names = FALSE)
  cat("\nNo confirmatory correlation tests or causal claims from 2–3 complete-case patients.\n")

  # Plot 1: matched patient Normal-to-Tumor score differences, 20+ cells in BOTH.
  paired_primary <- paired |>
    filter(variant == "primary_07e", min_cells_per_region == 20L) |>
    mutate(program = factor(program, levels = programs),
           tumor_region = factor(tumor_region, levels = tissue_regions))
  if (nrow(paired_primary) > 0) {
    plot_data <- bind_rows(
      paired_primary |>
        transmute(patient_id, program, tumor_region,
                  category = "Normal", mean_score = mean_normal),
      paired_primary |>
        transmute(patient_id, program, tumor_region,
                  category = "Tumor", mean_score = mean_tumor)
    ) |>
      mutate(category = factor(category, levels = c("Normal", "Tumor")))
    p_paired <- ggplot(plot_data,
                       aes(x = category, y = mean_score, group = patient_id)) +
      geom_line(alpha = 0.5, linewidth = 0.45, color = "gray55") +
      geom_point(aes(color = patient_id), size = 2) +
      facet_grid(program ~ tumor_region, scales = "free_y") +
      labs(title = "GSE144735: matched patient gene-program scores",
           subtitle = "Only >=20 cells in BOTH compared regions; Core+Border is one aggregate per patient",
           x = "Matched tissue", y = "Mean UCell score", color = "Patient") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            legend.position = "bottom")
    ggsave(file.path(figures_dir, "paired_tumor_normal_ge20.png"),
           p_paired, width = 14, height = 10, dpi = 300)
  }

  # Plot 2: patient-specific primary tumor-normal deltas.
  if (nrow(paired_primary) > 0) {
    p_delta <- ggplot(paired_primary,
                      aes(x = patient_id, y = delta_mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_point(size = 2) +
      facet_grid(program ~ tumor_region, scales = "free_y") +
      labs(title = "GSE144735: matched patient score differences",
           subtitle = "Tumor minus Normal, mean UCell; >=20 cells in both groups; descriptive",
           x = "Patient", y = "Paired delta") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1))
    ggsave(file.path(figures_dir, "matched_patient_deltas_ge20.png"),
           p_delta, width = 14, height = 10, dpi = 300)
  }

  cat("\nSaved tables: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
  cat("Saved figures: ", normalizePath(figures_dir, winslash = "/"), "\n", sep = "")
  cat("========== STEP 07G COMPLETE ==========\n")
})
