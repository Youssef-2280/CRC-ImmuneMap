# CRC-ImmuneMap — Step 09 integrity audit

- Timestamp: 2026-10-14 20:36:06 EEST
- Technical checks passed: 26/26
- Discovery patient-level associations: see discovery_summary_for_figures.csv (n=11 for primary comparisons).
- Validation: 6 patients in same published study, not an external-study cohort.
- Validation Tumor Core+Border is one within-patient aggregate, not two independent observations.
- High donor dominance/cross-lineage concerns remain; final lineage masks are NOT automatically approved by this audit.
- Do not claim validated cross-program correlation or causal signaling.

## Validation primary matched descriptive summary (>=20 cells each side)

| Program | Covered patients | Tumor higher | Median tumor-normal UCell delta |
|---|---:|---:|---:|
| SPP1_TAM | 4 | 4 | 0.344 |
| FAP_CAF | 6 | 6 | 0.228 |
| CD8_Dysfunction | 3 | 3 | 0.027 |
| CD8_Cytotoxicity | 3 | 2 | 0.024 |

## Required human review before publishing

- Confirm and document final cluster/cell inclusion masks for macrophages, fibroblasts and CD8.
- Explain exclusion/review of epithelial cross-lineage, donor-dominated and rare clusters.
- Assess whether all paired tissue labels genuinely represent matched samples for each patient.
- State explicitly that validation patient coverage is insufficient for a meaningful triple-program correlation.
- Add QC/doublet limitations and frozen signature provenance to the README.
- Include software versions, code, small derived tables and figures; exclude raw data and large RDS.

