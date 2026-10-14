# CRC-ImmuneMap | Step 04: annotation and patient eligibility audit
# Requires Seurat, dplyr, ggplot2. Run from the CRC-ImmuneMap project root.
# This script does not change the saved Seurat objects.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

# EDIT ONLY if your working directory is not the CRC-ImmuneMap root:
project_dir <- "C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap"
if (!dir.exists(project_dir)) {
  stop("Project folder not found. Edit project_dir at the top of this script.")
}
setwd(project_dir)

outdir <- "results/tables/audit"
figdir <- "results/figures/audit"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

configs <- list(
  Macrophages = list(
    path = "objects/crc_discovery_myeloid_annotated_v1.rds",
    cluster = "myeloid_clusters", state = "myeloid_state_preliminary",
    status = "myeloid_annotation_status", eligible = "macrophage_primary_eligible",
    markers = c("LYZ", "LST1", "C1QA", "C1QC", "APOE", "CD68", "SPP1", "TREM2", "GPNMB", "APOC1", "FOLR2", "LYVE1", "FCN1", "CD3D", "EPCAM")
  ),
  Fibroblasts = list(
    path = "objects/crc_discovery_fibroblast_annotated_v1.rds",
    cluster = "fibroblast_clusters", state = "fibroblast_state_preliminary",
    status = "fibroblast_annotation_status", eligible = "fibroblast_primary_eligible",
    markers = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP", "CTHRC1", "COL11A1", "POSTN", "INHBA", "ACTA2", "EPCAM", "CD3D")
  ),
  CD8 = list(
    path = "objects/crc_discovery_tcells_annotated_v1.rds",
    cluster = "tcell_clusters", state = "tcell_state_preliminary",
    status = "tcell_annotation_status", eligible = "cd8_primary_eligible",
    markers = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "CD4", "GZMK", "GZMB", "PRF1", "NKG7", "CCL5", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "LAYN", "ENTPD1", "FOXP3", "TRDC")
  )
)

missing_files <- vapply(configs, function(x) !file.exists(x$path), logical(1))
if (any(missing_files)) {
  stop("Missing checkpoint(s): ", paste(vapply(configs[missing_files], `[[`, "", "path"), collapse = ", "))
}

all_meta <- list()
all_clusters <- list()
all_markers <- list()
coverage <- list()

for (lineage in names(configs)) {
  conf <- configs[[lineage]]
  cat("\n=== Loading", lineage, "===\n")
  obj <- readRDS(conf$path)
  DefaultAssay(obj) <- "RNA"
  meta <- obj[[]]
  required <- c("patient_id", "sample_id", "condition", conf$cluster, conf$state, conf$status, conf$eligible,
                "nCount_RNA", "nFeature_RNA", "percent.mt")
  absent <- setdiff(required, colnames(meta))
  if (length(absent)) stop(lineage, ": missing metadata columns: ", paste(absent, collapse = ", "))
  if (!all(as.character(meta$condition) %in% c("Normal", "Tumor"))) {
    stop(lineage, ": unexpected condition labels. Inspect unique(obj$condition).")
  }
  if (anyNA(meta[, c("patient_id", "sample_id", "condition")])) {
    stop(lineage, ": NA in patient_id, sample_id or condition; inspect before proceeding.")
  }
  if (any(is.na(meta[[conf$eligible]]))) stop(lineage, ": NA in primary eligibility column")
  eligible <- as.logical(meta[[conf$eligible]])
  if (anyNA(eligible)) stop(lineage, ": eligibility values cannot be interpreted as TRUE/FALSE")

  meta_short <- data.frame(
    cell_id = rownames(meta), lineage = lineage,
    patient_id = as.character(meta$patient_id), sample_id = as.character(meta$sample_id),
    condition = as.character(meta$condition),
    cluster = as.character(meta[[conf$cluster]]),
    state = as.character(meta[[conf$state]]),
    status = as.character(meta[[conf$status]]),
    primary_eligible = eligible,
    nCount_RNA = meta$nCount_RNA,
    nFeature_RNA = meta$nFeature_RNA,
    percent.mt = meta$percent.mt,
    stringsAsFactors = FALSE
  )
  all_meta[[lineage]] <- meta_short

  # Cluster-level distribution, including clusters excluded from primary analysis
  cluster_counts <- meta_short %>%
    count(cluster, patient_id, name = "patient_cells")
  cluster_summary <- meta_short %>%
    group_by(cluster, state, status, primary_eligible) %>%
    summarise(n_cells = n(), n_patients = n_distinct(patient_id),
              tumor_cells = sum(condition == "Tumor"), normal_cells = sum(condition == "Normal"),
              median_features = median(nFeature_RNA), median_umis = median(nCount_RNA),
              median_mt = median(percent.mt), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(cluster_total = sum(n_cells),
           dominant_patient_pct = round(100 * max(cluster_counts$patient_cells[cluster_counts$cluster == first(cluster)]) / cluster_total, 1)) %>%
    ungroup() %>%
    mutate(lineage = lineage, .before = 1)
  all_clusters[[lineage]] <- cluster_summary

  # Raw count detection rates for a marker panel. This is annotation QC, not gene-program scoring.
  present <- intersect(conf$markers, rownames(obj))
  coverage[[lineage]] <- data.frame(lineage = lineage, gene = conf$markers,
                                    present = conf$markers %in% rownames(obj))
  if (length(present)) {
    for (cl in sort(unique(meta_short$cluster))) {
      cell_ids <- meta_short$cell_id[meta_short$cluster == cl]
      # Extract a small marker panel to avoid expanding the full count matrix.
      counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[present, cell_ids, drop = FALSE]
      rates <- 100 * Matrix::rowMeans(counts > 0)
      all_markers[[paste(lineage, cl, sep = "_")]] <- data.frame(
        lineage = lineage, cluster = cl, gene = names(rates),
        pct_detected = round(as.numeric(rates), 1), n_cells = length(cell_ids)
      )
    }
  }
  cat("Cells:", nrow(meta_short), "| primary eligible:", sum(eligible),
      "| patients:", n_distinct(meta_short$patient_id), "\n")
  rm(obj, meta, meta_short)
  invisible(gc())
}

cells <- bind_rows(all_meta)
# Use the broad annotated discovery object to include samples with zero selected cells.
broad_file <- "objects/crc_discovery_annotated_v2.rds"
if (!file.exists(broad_file)) stop("Missing broad annotation checkpoint: ", broad_file)
cat("\n=== Loading broad discovery metadata for complete sample denominators ===\n")
broad <- readRDS(broad_file)
broad_meta <- broad[[]]
if (!all(c("patient_id", "sample_id", "condition") %in% names(broad_meta))) {
  stop("Broad metadata lacks patient_id, sample_id, or condition.")
}
broad_mapping <- broad_meta %>%
  distinct(patient_id, sample_id, condition) %>%
  mutate(across(everything(), as.character))
broad_sample_totals <- broad_meta %>%
  count(patient_id, sample_id, condition, name = "n_total_qc_cells") %>%
  mutate(across(c(patient_id, sample_id, condition), as.character))
rm(broad, broad_meta)
invisible(gc())
cluster_summary <- bind_rows(all_clusters)
marker_summary <- bind_rows(all_markers)
coverage_table <- bind_rows(coverage)

# Consistency: one sample must map to one patient and one condition
mapping <- bind_rows(
  broad_mapping %>% select(sample_id, patient_id, condition),
  cells %>% distinct(sample_id, patient_id, condition)
) %>%
  distinct() %>%
  group_by(sample_id) %>%
  mutate(n_patient_ids_for_sample = n_distinct(patient_id),
         n_conditions_for_sample = n_distinct(condition)) %>%
  ungroup() %>% arrange(patient_id, sample_id)
if (any(mapping$n_patient_ids_for_sample != 1 | mapping$n_conditions_for_sample != 1)) {
  write.csv(mapping, file.path(outdir, "sample_mapping_check.csv"), row.names = FALSE)
  stop("Sample-to-patient/condition mapping conflict. Review sample_mapping_check.csv.")
}

# Main descriptive table: selected cells per sample and compartment.
selected <- cells %>% filter(primary_eligible)
sample_counts <- selected %>%
  count(patient_id, sample_id, condition, lineage, name = "n_cells")

# Show all three compartments even when a sample has no selected cells.
sample_base <- mapping %>% select(patient_id, sample_id, condition)
sample_complete <- tidyr::crossing(sample_base, lineage = names(configs)) %>%
  left_join(sample_counts, by = c("patient_id", "sample_id", "condition", "lineage")) %>%
  left_join(broad_sample_totals, by = c("patient_id", "sample_id", "condition")) %>%
  mutate(n_cells = tidyr::replace_na(n_cells, 0L)) %>%
  arrange(patient_id, condition, sample_id, lineage)

# Screen for samples with enough cells to compute reasonably stable *descriptive* means.
# A 20-cell threshold is a pragmatic audit flag, NOT a validated hard cutoff.
min_cells_flag <- 20L
sample_complete <- sample_complete %>% mutate(audit_at_least_20 = n_cells >= min_cells_flag)
patient_counts <- sample_complete %>%
  group_by(patient_id, condition, lineage) %>%
  summarise(n_cells = sum(n_cells), n_samples = sum(n_cells > 0), .groups = "drop")
patient_wide <- patient_counts %>%
  tidyr::pivot_wider(names_from = lineage, values_from = n_cells, values_fill = 0L)
for (nm in names(configs)) if (!nm %in% names(patient_wide)) patient_wide[[nm]] <- 0L
patient_wide <- patient_wide %>%
  mutate(all_three_present = Macrophages > 0 & Fibroblasts > 0 & CD8 > 0,
         all_three_at_least_20 = Macrophages >= 20 & Fibroblasts >= 20 & CD8 >= 20) %>%
  arrange(patient_id, condition)

paired_candidates <- mapping %>%
  distinct(patient_id, condition) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = condition, values_from = present, values_fill = FALSE)
for (nm in c("Normal", "Tumor")) if (!nm %in% names(paired_candidates)) paired_candidates[[nm]] <- FALSE
paired_candidates <- paired_candidates %>%
  mutate(potential_paired_patient = Normal & Tumor) %>% arrange(patient_id)

write.csv(mapping, file.path(outdir, "sample_mapping_check.csv"), row.names = FALSE)
write.csv(cluster_summary, file.path(outdir, "cluster_patient_qc_summary.csv"), row.names = FALSE)
write.csv(marker_summary, file.path(outdir, "marker_detection_by_cluster.csv"), row.names = FALSE)
write.csv(coverage_table, file.path(outdir, "marker_gene_coverage.csv"), row.names = FALSE)
write.csv(sample_complete, file.path(outdir, "primary_cells_by_sample.csv"), row.names = FALSE)
write.csv(patient_wide, file.path(outdir, "primary_cells_by_patient_condition.csv"), row.names = FALSE)
write.csv(paired_candidates, file.path(outdir, "potential_paired_patients.csv"), row.names = FALSE)

# Plot: aggregate selected cell coverage by patient and condition (not cell fractions).
coverage_for_plot <- patient_counts
p <- ggplot(coverage_for_plot, aes(x = patient_id, y = lineage, fill = log1p(n_cells))) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~condition, scales = "free_x") +
  scale_fill_viridis_c(name = "log(1 + cells)") +
  labs(title = "Primary compartments: cell coverage by patient",
       subtitle = "Counts summed across a patient's samples within each condition",
       x = "Patient", y = "Primary compartment") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(figdir, "primary_cell_coverage.png"), p,
       width = 13, height = 6, dpi = 300)

cat("\n========== AUDIT SUMMARY ==========\n")
print(cells %>% group_by(lineage) %>%
        summarise(total_cells = n(), primary_cells = sum(primary_eligible),
                  n_patients = n_distinct(patient_id), .groups = "drop"))
cat("\nPRIMARY CELLS BY PATIENT AND CONDITION:\n")
print(patient_wide, n = Inf, width = Inf)
cat("\nPOTENTIAL PAIRED PATIENTS (identity match only; check clinical metadata):\n")
print(paired_candidates, n = Inf, width = Inf)
cat("\nFLAGGED CLUSTERS (dominant donor > 70%, or annotation not marker-supported):\n")
print(cluster_summary %>%
        filter(dominant_patient_pct > 70 | status != "Marker_supported") %>%
        select(lineage, cluster, state, status, primary_eligible, n_cells,
               n_patients, dominant_patient_pct), n = Inf, width = Inf)
cat("\nMISSING PANEL GENES:\n")
print(coverage_table %>% filter(!present), n = Inf, width = Inf)
cat("\nSAVED: ", normalizePath(outdir, winslash = "/"), "\n", sep = "")
cat("NOTE: This script evaluates annotation and sample coverage only; no hypothesis test yet.\n")
