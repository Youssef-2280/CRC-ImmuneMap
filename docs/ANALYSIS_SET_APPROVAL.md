# Final primary analysis set — approved selection policy

**Approval on 2026-10-14:** The project user explicitly approved retaining **`primary_07e`** as the final *primary analysis selection* while preserving donor-dominance concerns, cross-lineage concerns, sensitivity analyses, and limitations. This records a choice of *analytical cohort*, **not** proof of histological cell identity, absence of doublets, or clinical validation.

| Compartment | Inclusion policy (validation) | Selected cells in original 07E/07F scoring |
| --- | --- | ---: |
| Macrophage candidate pool | Myeloid clusters **0, 1, 2, 5**, defined before scoring; macrophage-marker/published-subtype evidence considered; do not select on SPP1 score | 1,745 |
| Fibroblast candidate pool | Stromal clusters **0, 2, 3, 4, 9, 11, 12**; do not select on FAP score | 4,230 |
| CD8 candidate pool | T clusters **2, 4, 10** with an additional per-cell author-subtype **`CD8+ T cells`** gate (723 + 381 + 160); do not select on checkpoint score | 1,264 |

The source 07E/07F per-cell masks were built and scored on the user's machine. In this release, the *policy* and aggregated numerical results are frozen; the *complete per-cell ID list* is **not present in these conversation uploads**. Execute `scripts/00_export_and_verify_primary_07e_masks.R` in the original project to materialize/checksum the actual cell-ID masks **without recalculating UCell**. Until that is run, do not claim the repository contains cell-by-cell masks.

## Interpretation of cluster-level policy

- `final_analysis_set_cluster_policy.csv` records the approved **selection action** for each of 36 clusters; it does **not** certify that every cell in an included cluster is a pure cell type.
- For T cells, selected cluster numbers alone are **insufficient**: require the published `CD8+ T cells` gate for individual cells.
- Source flags are retained. Donor-dominated myeloid 0 (KUL01, 94.7%), myeloid 1 (KUL21, 77.4%), and T 10 (KUL19, 92.8%) remain concerns.
- High epithelial co-detection in macrophage cluster 0 and stromal clusters 0/12 does not establish either doublets or ambient RNA as the explanation. 
- High-risk cluster omission and cross-lineage exclusion are *sensitivity* analyses; they do not replace the approved primary set.
- User approval is recorded from the conversation; no unrelated reviewer signature, initials, or independent biological adjudication is claimed.
