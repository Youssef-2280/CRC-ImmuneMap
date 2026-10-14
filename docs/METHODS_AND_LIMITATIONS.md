# Methods and limitations — evidence-grounded release notes

## Cohorts and unit of analysis

- Discovery: GEO **GSE132465**, 23 available patients; main cross-program analysis uses **11 complete-case tumor patients** with at least 20 eligible cells per compartment.
- Validation: GEO **GSE144735**, **6 different patients**, sampled in Normal, Tumor Core and Tumor Border (18 regional samples). It is a distinct patient cohort **within the same published study**, not a separate independent publication. Core and Border are never counted as independent patients.
- Patient-level rather than cell-level summaries define statistical units, guarding against pseudoreplication from thousands of cells.

## Program scores and comparisons

The same frozen **29 distinct program–gene entries** are used in both datasets: SPP1_TAM (7), FAP_CAF (7), CD8_Dysfunction (8), CD8_Cytotoxicity (7). Scores use **UCell** with **maxRank=1500**; the supplied manifest has the exact genes. The validation Step 07F scored candidate cells using raw counts, all genes retained for within-cell ranking, and then summarized scores per patient/region. The final primary mask is the previously scored `primary_07e`.

- Discovery: patient-level Spearman correlations across 11 complete-case tumor patients; permutation p-values, Benjamini–Hochberg adjustment for the three hypotheses.
- Validation: within-patient tumor-minus-normal *mean* UCell score, requiring **>=20 eligible cells on each compared side**. Paired comparisons of Tumor Core, Tumor Border and combined Core+Border are alternative **within-patient** views; do not pool as independent replications. Median paired delta summarizes patients, not cells.
- Sensitivities: remove cross-lineage-flagged cells; omit donor/crosslineage high-risk clusters; broader published-CD8 sensitivity.
- **No confirmatory triad correlation test is supportable in validation**: only 2/6 patients have >=20 cells from all three relevant pools in Tumor Core; 3/6 in combined Core+Border.

## Interpretation restrictions

- None of three Discovery cross-program hypotheses reached a significant BH-adjusted p value.
- Directionally repeated SPP1_TAM and FAP_CAF tumor-normal contrasts do **not** establish patient-level SPP1/FAP correlation, co-localization, cellular communication, causation or a new external-study validation.
- Cell annotations are an *approved analysis set* with uncertainty, not microscopy- or genotype-verified identities.
- Strong donor concentration may confound cluster phenotype; sensitivity exclusions can reduce effective patient n (SPP1 from 4 to 2 matched patients in combined tumor comparison).
- Expression of cross-lineage genes does not alone prove a doublet; possible ambient RNA, genuine expression and contamination cannot be separated unequivocally from these summary tables.
- No direct spatial interaction, ligand–receptor communication, functional perturbation or treatment outcome effect was tested here.
- Technical Step 09 integrity checks were **26/26 PASS** for supplied version, but do not independently validate biological annotation or source pairing claims. Final source patient/sample mapping and raw QC should be reviewed before manuscript submission.

## Software environment (captured in user's session)

R 4.6.1 (Windows 11); Seurat 5.5.1; SeuratObject 5.4.0; UCell 2.16.0; Matrix 1.7-5; dplyr 1.2.1; ggplot2 4.0.3. Complete user-side session info: `evidence/R_session_info.txt`. Figures in this package use Python packages in `requirements-figures.txt`.

## Source-verification addendum (2026-10-14)

The original source publication is **Lee et al., Nature Genetics (2020), DOI 10.1038/s41588-020-0636-z**. GEO 51 sample IDs were checked against project-derived IDs (see `docs/SOURCE_METADATA_VERIFICATION.md`). KUL3 source QC included UMIs >1000, genes 200–6000, mitochondrial expression <20%, followed by source-author lineage-consensus/outlier filtering. These source-author procedures cannot independently certify the absence of residual doublets in the downstream selected cell populations. See `docs/ANNOTATION_QC_DOUBLET_ASSESSMENT.md` and the publication text in `docs/MANUSCRIPT_READY_METHODS_LIMITATIONS.md`.
