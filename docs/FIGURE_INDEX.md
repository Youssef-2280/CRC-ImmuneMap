# Figure index — PNG only

All 26 images in this release are **PNG only**: 5 main results figures and 21 supplementary single-cell figures. There are no SVG copies and no AI-generated graphical workflow. The old faint discovery UMAPs were replaced by R-rendered high-contrast plots; a split tumor/normal view was added.

## Main results figures

| PNG file | Simple title |
|---|---|
| [01_discovery_associations.png](../figures/main/01_discovery_associations.png) | Discovery cross-program correlations |
| [02_validation_scores.png](../figures/main/02_validation_scores.png) | Matched tumor versus normal scores |
| [03_paired_differences.png](../figures/main/03_paired_differences.png) | Paired patient differences by region |
| [04_sensitivity_analysis.png](../figures/main/04_sensitivity_analysis.png) | Cell-selection sensitivity |
| [05_patient_coverage.png](../figures/main/05_patient_coverage.png) | Validation patient coverage |

## Supplementary single-cell figures

| PNG file | Simple title |
|---|---|
| [s01_qc_before_condition.png](../figures/supplementary/s01_qc_before_condition.png) | Pre-filter QC by condition |
| [s02_qc_before_sample.png](../figures/supplementary/s02_qc_before_sample.png) | Pre-filter QC by sample |
| [s03_qc_after_filtering.png](../figures/supplementary/s03_qc_after_filtering.png) | Post-filter QC (discovery) |
| [s04_cell_types_umap.png](../figures/supplementary/s04_cell_types_umap.png) | Discovery broad cell-type UMAP |
| [s05_tumor_normal_umap.png](../figures/supplementary/s05_tumor_normal_umap.png) | Discovery tumor and normal UMAP |
| [s06_tumor_normal_split.png](../figures/supplementary/s06_tumor_normal_split.png) | Discovery tumor and normal split UMAP |
| [s07_cell_type_markers.png](../figures/supplementary/s07_cell_type_markers.png) | Discovery broad marker DotPlot |
| [s08_myeloid_umap.png](../figures/supplementary/s08_myeloid_umap.png) | Discovery myeloid UMAP |
| [s09_fibroblast_umap.png](../figures/supplementary/s09_fibroblast_umap.png) | Discovery fibroblast UMAP |
| [s10_tcell_umap.png](../figures/supplementary/s10_tcell_umap.png) | Discovery T-cell UMAP |
| [s11_validation_qc.png](../figures/supplementary/s11_validation_qc.png) | Validation authors-processed QC |
| [s12_validation_cell_types.png](../figures/supplementary/s12_validation_cell_types.png) | Validation cell-type composition |
| [s13_validation_myeloid_clusters.png](../figures/supplementary/s13_validation_myeloid_clusters.png) | Validation myeloid clusters |
| [s14_validation_stromal_clusters.png](../figures/supplementary/s14_validation_stromal_clusters.png) | Validation stromal clusters |
| [s15_validation_tcell_clusters.png](../figures/supplementary/s15_validation_tcell_clusters.png) | Validation T-cell clusters |
| [s16_validation_myeloid_markers.png](../figures/supplementary/s16_validation_myeloid_markers.png) | Validation myeloid marker DotPlot |
| [s17_validation_stromal_markers.png](../figures/supplementary/s17_validation_stromal_markers.png) | Validation stromal marker DotPlot |
| [s18_validation_tcell_markers.png](../figures/supplementary/s18_validation_tcell_markers.png) | Validation T-cell marker DotPlot |
| [s19_macrophage_selection.png](../figures/supplementary/s19_macrophage_selection.png) | Primary myeloid analysis selection UMAP |
| [s20_fibroblast_selection.png](../figures/supplementary/s20_fibroblast_selection.png) | Primary stromal analysis selection UMAP |
| [s21_cd8_selection.png](../figures/supplementary/s21_cd8_selection.png) | Primary CD8 analysis selection UMAP |

## Provenance and caveats

- Main 5 plots: rendered in Python/matplotlib from numerical CSV summaries exported from R. Not newly computed statistics.
- Supplementary: 2 researcher-supplied original pre-filter QC PNGs; 18 Seurat-generated panels from the supplied RDS images with the two low-contrast UMAPs replaced; plus one split high-contrast view. Thus 21 distinct supplementary PNG views in this release.
- A few cluster annotation panels show **preliminary/unreviewed** subcluster labels. Do not interpret them as experimentally validated cell identity.
- A post-QC plot is not evidence of the unfiltered starting-cell distribution; validation processed QC is based on the original published processed data.
- Exact per-cell masks are not packaged with this repo; see the R mask-export script.
- Original R output filenames and new paths are in [FIGURE_FILE_MAPPING.csv](FIGURE_FILE_MAPPING.csv).
- Original R figure generation source manifest is preserved as [R_FIGURE_SOURCE_MANIFEST_ORIGINAL.csv](R_FIGURE_SOURCE_MANIFEST_ORIGINAL.csv).
