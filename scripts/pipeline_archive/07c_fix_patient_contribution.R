# CRC-ImmuneMap | 07C Repair: correct patient contribution summaries
# Reads existing Step 07C CSVs. Does NOT rerun clustering, markers or Seurat.
# Fixes a summarise() masking bug which falsely reported KUL01 = 100% in every cluster.

project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
stopifnot(dir.exists(project_dir))
setwd(project_dir)

source_dir <- "results/tables/validation/07c"
output_dir <- file.path(source_dir, "patient_contribution_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (lineage in c("myeloids", "stromal", "tcells")) {
  cat("\n=== Correcting", lineage, "===\n")
  p <- file.path(source_dir,
                 paste0(lineage, "_cluster_by_patient_sample_subtype.csv"))
  stopifnot(file.exists(p))
  dat <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  cols <- c("validation_cluster", "patient_id", "sample_id",
            "tissue_region", "n_cells")
  stopifnot(all(cols %in% names(dat)))
  stopifnot(!anyNA(dat[, cols]), all(dat$n_cells >= 0))

  # Each entry is already aggregated by cluster, patient, sample, subtype;
  # sum across subtype/sample to obtain cluster-patient counts.
  donor <- stats::aggregate(n_cells ~ validation_cluster + patient_id,
                            data = dat, FUN = sum)
  donor <- donor[order(as.integer(as.character(donor$validation_cluster)),
                       donor$patient_id), , drop = FALSE]
  totals <- stats::aggregate(n_cells ~ validation_cluster, data = donor, FUN = sum)
  names(totals)[2] <- "cluster_total"
  donor <- merge(donor, totals, by = "validation_cluster", sort = FALSE)
  donor$pct_of_cluster <- round(100 * donor$n_cells / donor$cluster_total, 2)
  donor <- donor[order(as.integer(as.character(donor$validation_cluster)),
                       -donor$n_cells, donor$patient_id), , drop = FALSE]

  splits <- split(donor, donor$validation_cluster)
  summary <- do.call(rbind, lapply(splits, function(x) {
    x <- x[order(-x$n_cells, x$patient_id), , drop = FALSE]
    data.frame(
      validation_cluster = as.character(x$validation_cluster[1]),
      n_cells = sum(x$n_cells),
      n_patients = sum(x$n_cells > 0),
      dominant_patient = x$patient_id[1],
      dominant_patient_n = x$n_cells[1],
      dominant_patient_pct = round(100 * x$n_cells[1] / sum(x$n_cells), 1),
      second_patient_pct = if (nrow(x) >= 2) round(100 * x$n_cells[2] / sum(x$n_cells), 1) else NA_real_,
      flagged_donor_gt70 = (100 * x$n_cells[1] / sum(x$n_cells)) > 70,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary <- summary[order(as.integer(summary$validation_cluster)), , drop = FALSE]

  # Verify against original cluster sizes and region totals.
  old_path <- file.path(source_dir, paste0(lineage, "_cluster_patient_summary.csv"))
  if (file.exists(old_path)) {
    old <- read.csv(old_path, stringsAsFactors = FALSE)
    a <- summary[order(as.integer(summary$validation_cluster)), ]
    b <- old[order(as.integer(old$validation_cluster)), ]
    stopifnot(identical(as.character(a$validation_cluster),
                        as.character(b$validation_cluster)),
              all(a$n_cells == b$n_cells),
              sum(a$n_cells) == sum(old$n_cells))
    region_cols <- intersect(c("Normal", "Tumor_Border", "Tumor_Core"), names(old))
    if (length(region_cols)) {
      stopifnot(all(rowSums(old[, region_cols, drop = FALSE]) == old$n_cells))
    }
  }

  write.csv(donor, file.path(output_dir,
           paste0(lineage, "_cluster_patient_counts_corrected.csv")), row.names = FALSE)
  write.csv(summary, file.path(output_dir,
           paste0(lineage, "_cluster_patient_summary_corrected.csv")), row.names = FALSE)
  cat("Total cells:", sum(summary$n_cells), "| clusters:", nrow(summary),
      "| >70% clusters:", sum(summary$flagged_donor_gt70), "\n")
  print(summary, row.names = FALSE)
}
cat("\n=== DONE: corrected donor summaries ===\n")
cat("Saved to: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
cat("Do not use original *_cluster_patient_summary.csv dominant_patient columns.\n")
