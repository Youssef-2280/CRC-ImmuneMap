# CRC-ImmuneMap

**Single-Cell Analysis of Macrophage, Fibroblast, and CD8 T-Cell Programs in Colorectal Cancer**

## Project Overview

CRC-ImmuneMap is an exploratory single-cell RNA-seq project investigating SPP1-associated macrophages, FAP-associated fibroblasts, and CD8 T-cell dysfunction in colorectal cancer (CRC).

The project includes cell-type annotation, gene-program scoring, patient-level analysis, and cross-cohort descriptive validation.

## Datasets

| Dataset | Cohort | Patients |
|---|---|---|
| GSE132465 | Discovery (SMC) | 23 patients; 11 in primary association analysis |
| GSE144735 | Validation (KUL3) | 6 patients |

Both cohorts originate from the same published study and represent different patient groups.

## Analysis Workflow

1. Single-cell data preprocessing and quality control
2. Broad cell-type annotation
3. Myeloid, fibroblast, and T-cell subclustering
4. Lineage-specific cell-selection review
5. Frozen UCell gene-program scoring
6. Patient-level association and paired tissue analyses
7. Sensitivity analyses and cross-cohort evidence synthesis

The four evaluated programs were **SPP1_TAM**, **FAP_CAF**, **CD8_Dysfunction**, and **CD8_Cytotoxicity**.

## Key Findings

**Discovery cohort (GSE132465):**

- SPP1_TAM versus FAP_CAF: Spearman ρ = 0.400; BH-adjusted p = 0.669.
- SPP1_TAM versus CD8_Dysfunction: ρ = 0.082; BH-adjusted p = 0.817.
- FAP_CAF versus CD8_Dysfunction: ρ = −0.182; BH-adjusted p = 0.817.

No tested cross-program association reached statistical significance.

**Validation cohort (GSE144735):**

- SPP1_TAM scores were higher in tumor tissue in 4/4 eligible matched patients.
- FAP_CAF scores were higher in tumor tissue in 6/6 eligible matched patients.
- CD8_Dysfunction scores were higher in tumor tissue in 3/3 eligible matched patients.

These are descriptive observations. Limited patient coverage prevents a reliable confirmatory test of the proposed three-program relationship.

## Main Figures

All figures are provided in **PNG format only**.

| Figure | Description |
|---|---|
| [Figure 01](figures/main/01_discovery_associations.png) | Discovery patient-level associations |
| [Figure 02](figures/main/02_validation_scores.png) | Matched validation program scores |
| [Figure 03](figures/main/03_paired_differences.png) | Paired tumor-normal score differences |
| [Figure 04](figures/main/04_sensitivity_analysis.png) | Cell-selection sensitivity analysis |
| [Figure 05](figures/main/05_patient_coverage.png) | Patient coverage |

## Supplementary Single-Cell Figures

The repository contains **21 supplementary PNG figures**, including quality-control plots, cell-type UMAPs, marker DotPlots, validation subclustering, and primary cell-selection visualizations.

[Browse Supplementary Figures](figures/supplementary/)

[Complete Figure Index](docs/FIGURE_INDEX.md)

## Repository Structure

- **figures/main/** — Five main results figures.
- **figures/supplementary/** — Twenty-one single-cell figures.
- **scripts/** — Analysis scripts and figure-generation code.
- **data/derived/** — Processed summary tables.
- **data/manifests/** — Frozen gene-program signatures.
- **docs/** — English results, methods, figure legends, and limitations.
- **evidence/** — Annotation review and quality-control evidence.
- **references/** — Bibliographic references.

## Reproducibility

Statistical analyses and UCell scoring were performed in R. The five main figures were generated in Python using numerical results exported from R. Supplementary single-cell figures were produced using R/Seurat or retained from original quality-control outputs.

The repository contains derived tables and analysis scripts, but does not include large raw sequencing matrices or Seurat RDS objects.

## Study Limitations

This is an exploratory computational analysis.

The main limitations include small patient-level sample sizes, donor imbalance, potential cross-lineage transcript contamination, uncertain residual doublets, and incomplete coverage of all three cell compartments in the validation cohort.

Increased tumor-associated SPP1 and FAP program activity does not establish a direct biological interaction or a statistically validated SPP1–FAP–CD8 dysfunction axis.

## References

Lee et al. (2020), *Nature Genetics*.

https://doi.org/10.1038/s41588-020-0636-z

Andreatta & Carmona (2021), UCell.

https://doi.org/10.1016/j.csbj.2021.06.043

---

**Project status:** Exploratory research repository with documented analyses, figures, and methodological limitations.
