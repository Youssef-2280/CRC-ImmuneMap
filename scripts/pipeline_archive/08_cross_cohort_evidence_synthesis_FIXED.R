# CRC-ImmuneMap — STEP 08
# Transparent discovery/validation evidence synthesis (no new hypothesis tests)
# Reads EXISTING CSV files only; never loads Seurat or changes the prior results.
# Evidence categories are distinct:
#  (1) Discovery patient-level cross-program associations
#  (2) Validation matched Tumor-vs-Normal program differences
# They are NOT interchangeable evidence for a three-program patient-level axis.

local({
  project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
  stopifnot(dir.exists(project_dir))
  setwd(project_dir)

  discovery_file <- "results/tables/associations/primary_three_correlations_ge20.csv"
  validation_file <- "results/tables/validation/07g_paired_descriptive/descriptive_paired_summary_by_variant.csv"
  frozen_manifest <- "results/tables/validation/frozen_discovery_signature_manifest_v1.csv"
  discovery_manifest <- "results/tables/programs/candidate_signature_manifest_v1.csv"
  out_tables <- "results/tables/integrated_evidence"
  out_figures <- "results/figures/integrated_evidence"
  dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)

  input_files <- c(discovery_file, validation_file, frozen_manifest, discovery_manifest)
  not_found <- input_files[!file.exists(input_files)]
  if (length(not_found)) stop("Required saved files missing:\n", paste(not_found, collapse = "\n"))

  disc <- read.csv(discovery_file, stringsAsFactors = FALSE, check.names = FALSE)
  val <- read.csv(validation_file, stringsAsFactors = FALSE, check.names = FALSE)
  sig_disc <- read.csv(discovery_manifest, stringsAsFactors = FALSE, check.names = FALSE)
  sig_val <- read.csv(frozen_manifest, stringsAsFactors = FALSE, check.names = FALSE)

  disc_cols <- c("comparison", "n_patients", "rho", "p_permutation", "p_BH_permutation")
  val_cols <- c("variant", "min_cells_per_region", "tumor_region", "program",
                "n_matched_patients", "n_tumor_higher", "n_normal_higher", "n_equal",
                "median_paired_delta")
  if (!all(disc_cols %in% names(disc))) stop("Unexpected discovery correlations CSV columns")
  if (!all(val_cols %in% names(val))) stop("Unexpected validation paired summary CSV columns")

  # Safety: ascertain that the signature in validation really is the same gene list.
  # Compare keys only, not filename/version labels, which can legitimately differ.
  sig_fields <- c("program", "gene", "method", "max_rank")
  if (!all(sig_fields %in% names(sig_disc)) || !all(sig_fields %in% names(sig_val))) {
    stop("Signature manifests are missing a gene/method/max_rank field; inspect before proceeding")
  }
  sig_key <- function(df) {
    x <- unique(df[, sig_fields, drop = FALSE])
    for (nm in sig_fields) x[[nm]] <- as.character(x[[nm]])
    x <- x[do.call(order, unname(x)), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  signature_identical <- identical(sig_key(sig_disc), sig_key(sig_val))
  if (!signature_identical) {
    write.csv(sig_key(sig_disc), file.path(out_tables, "discovery_signature_check.csv"), row.names = FALSE)
    write.csv(sig_key(sig_val), file.path(out_tables, "validation_signature_check.csv"), row.names = FALSE)
    stop("Discovery and validation signatures differ — stop, inspect signature_check files")
  }
  stopifnot(nrow(sig_key(sig_disc)) == 29L)

  # Preserve exact primary discovery results — no testing or refitting.
  discovery_primary <- disc[, disc_cols, drop = FALSE]
  stopifnot(nrow(discovery_primary) == 3L,
            all(discovery_primary$n_patients == 11L),
            all(is.finite(discovery_primary$rho)))
  write.csv(discovery_primary,
            file.path(out_tables, "discovery_primary_associations.csv"), row.names = FALSE)

  # Match only primary pre-defined validation selection and >=20 cells in both tissues.
  programs <- c("SPP1_TAM", "FAP_CAF", "CD8_Dysfunction", "CD8_Cytotoxicity")
  regions <- c("Tumor_Core", "Tumor_Border", "Tumor_Core_plus_Border")
  primary <- val[val$variant == "primary_07e" & val$min_cells_per_region == 20L, , drop = FALSE]
  expected_key <- expand.grid(tumor_region = regions, program = programs, stringsAsFactors = FALSE)
  key <- paste(primary$tumor_region, primary$program, sep = "___")
  if (nrow(primary) != nrow(expected_key) || anyDuplicated(key) ||
      !setequal(key, paste(expected_key$tumor_region, expected_key$program, sep = "___"))) {
    stop("Primary validation summary does not have expected program-region combinations")
  }
  if (any(primary$n_tumor_higher + primary$n_normal_higher + primary$n_equal !=
          primary$n_matched_patients)) stop("Direction counts do not sum to matched patients")
  primary <- primary[order(match(primary$tumor_region, regions),
                           match(primary$program, programs)), , drop = FALSE]
  write.csv(primary,
            file.path(out_tables, "validation_primary_paired_ge20.csv"), row.names = FALSE)

  # Sensitivity comparisons without optimizing the inclusion mask or threshold.
  sensitivity <- val[val$min_cells_per_region == 20L &
                       val$tumor_region == "Tumor_Core_plus_Border", , drop = FALSE]
  sensitivity <- sensitivity[order(sensitivity$program, sensitivity$variant), , drop = FALSE]
  write.csv(sensitivity,
            file.path(out_tables, "validation_selection_sensitivity_combined_ge20.csv"),
            row.names = FALSE)

  # A compact console view of what is scientifically supported.
  cat("\n========== STEP 08: CROSS-COHORT EVIDENCE ==========")
  cat("\nFrozen signature manifest check: PASS (29 unique program-gene pairs, settings matched)\n")
  cat("\nDISCOVERY: associations across patients (n = 11)\n")
  print(discovery_primary, row.names = FALSE)
  cat("\nVALIDATION: paired within-patient Tumor vs Normal (primary; >=20 cells)\n")
  print(primary[, c("tumor_region", "program", "n_matched_patients", "n_tumor_higher",
                    "n_normal_higher", "median_paired_delta")], row.names = FALSE)
  cat("\nVALIDATION: Core+Border sensitivity to inclusion masks\n")
  print(sensitivity[, c("variant", "program", "n_matched_patients", "n_tumor_higher",
                        "n_normal_higher", "median_paired_delta")], row.names = FALSE)
  cat("\nIMPORTANT: Validation Tumor-vs-Normal direction does not test SPP1-FAP-CD8 correlations.\n")
  cat("Core and Border are regions within the same patient; do not count as independent patients.\n")

  fmt <- function(x, digits = 3) format(round(x, digits), nsmall = digits, trim = TRUE)
  primary_combined <- primary[primary$tumor_region == "Tumor_Core_plus_Border", ]
  main_lines <- vapply(programs, function(p) {
    x <- primary_combined[primary_combined$program == p, ]
    paste0("- **", p, "**: tumor higher in ", x$n_tumor_higher, "/",
           x$n_matched_patients, " matched patients; median tumor-normal UCell delta = ",
           fmt(x$median_paired_delta), ".")
  }, character(1))
  main_disc <- vapply(seq_len(nrow(discovery_primary)), function(i) {
    x <- discovery_primary[i, ]
    paste0("- **", x$comparison, "**: Spearman rho = ", fmt(x$rho),
           "; permutation p = ", fmt(x$p_permutation),
           "; BH-adjusted p = ", fmt(x$p_BH_permutation),
           "; n = ", x$n_patients, ".")
  }, character(1))
  report <- c(
    "# CRC-ImmuneMap — Evidence synthesis (descriptive; provisional annotation)",
    "", "## Analysis units and cohorts", "",
    "- Discovery: GSE132465, 23 patients available; primary complete-case association analysis includes 11 patients with at least 20 eligible cells per compartment.",
    "- Validation: GSE144735, 6 different patients, 18 regions/samples (Normal, Tumor Core, Tumor Border). Core and Border are **not** two independent patients.",
    "- Validation patients belong to a separate cohort from the **same underlying published study**, not an entirely separate independent publication.",
    "- Candidate inclusion masks are **provisional**. This report does not convert them to final curated labels.",
    "- Frozen UCell signatures checked against both manifests: 29 program-gene pairs; max_rank=1500 in the saved manifests.",
    "", "## Discovery: cross-program associations", "", main_disc,
    "", "**Interpretation:** No prespecified association was significant after BH correction. SPP1-TAM and FAP-CAF showed a positive exploratory trend but no validated patient-level association.",
    "", "## Validation: within-patient tissue comparisons", "",
    "Primary selection, at least 20 cells in *each* matched compared group; Core+Border is one within-patient aggregate:",
    "", main_lines,
    "", "**Interpretation:** Increased SPP1-TAM and FAP-CAF scores in tumors were directionally consistent among covered matched patients, but this does NOT establish cross-program correlation, cell-cell interaction or causation. The CD8 result is particularly sparse and sensitive to cluster selection.",
    "", "## Sensitivity and limitations", "",
    "- Cross-lineage exclusion preserves the positive SPP1/FAP paired direction among eligible patients, but changes the sample composition and estimates.",
    "- Exclusion of high-risk clusters leaves just 2 matched patients for SPP1 in combined Core+Border at >=20 cells, reducing interpretability.",
    "- Complete three-program tumor coverage in GSE144735 is only 2 patients for Core and 3 for Core+Border at >=20 cells. **No confirmatory three-program correlation test is supportable.**",
    "- Tumor-minus-normal paired estimates for different programs use differing patients and cell numbers, so magnitude comparisons between programs should not be treated as standardized effect sizes.",
    "- Differential expression, ligand-receptor communication, and causal mechanisms have not been established by these results.",
    "", "## Next reproducibility checks", "",
    "- Freeze final lineage inclusion/exclusion rules with marker evidence and explicit handling of cross-lineage, donor-dominated and rare clusters.",
    "- Save sample/patient mapping and confirm no repeated sampling is treated as independent subjects.",
    "- Keep raw counts and large RDS outside public GitHub; publish small derived tables, scripts, dependency versions, and a limitations section.",
    "- Report no significant cross-program replication claim unless supported by a sufficiently powered independent cohort.",
    ""
  )
  writeLines(report, con = file.path(out_tables, "integrated_evidence_summary.md"), useBytes = TRUE)

  # Optional compact figure, genuinely descriptive and with explicit denominators.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    plot_data <- primary
    plot_data$program <- factor(plot_data$program, levels = rev(programs))
    plot_data$tumor_region <- factor(plot_data$tumor_region, levels = regions)
    plot_data$proportion_higher <- plot_data$n_tumor_higher / plot_data$n_matched_patients
    plot_data$direction_label <- paste0(plot_data$n_tumor_higher, "/", plot_data$n_matched_patients)
    p <- ggplot(plot_data, aes(x = tumor_region, y = program)) +
      geom_point(aes(size = n_matched_patients, fill = proportion_higher),
                 shape = 21, colour = "gray20") +
      geom_text(aes(label = direction_label), size = 3, vjust = -1.2) +
      scale_size_continuous(name = "Matched patients", range = c(4, 10)) +
      scale_fill_gradient(low = "white", high = "gray25", name = "Fraction tumor higher", limits = c(0, 1)) +
      coord_cartesian(clip = "off") +
      labs(title = "GSE144735: directions of paired tumor-normal scores",
           subtitle = "Primary candidate inclusion; >=20 cells in both matched groups; labels = tumor higher / eligible patients",
           x = "Tumor region compared with Normal", y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.margin = margin(12, 10, 10, 10))
    ggsave(file.path(out_figures, "validation_matched_direction_overview_ge20.png"),
           p, width = 12, height = 6.5, dpi = 300)
  } else {
    cat("ggplot2 not installed: summary tables and markdown are saved; plot skipped.\n")
  }
  cat("\nSaved tables/report to: ", normalizePath(out_tables, winslash = "/"), "\n", sep = "")
  cat("========== STEP 08 COMPLETE ==========\n")
})
