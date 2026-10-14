# CRC-ImmuneMap — Evidence synthesis (descriptive; provisional annotation)

## Analysis units and cohorts

- Discovery: GSE132465, 23 patients available; primary complete-case association analysis includes 11 patients with at least 20 eligible cells per compartment.
- Validation: GSE144735, 6 different patients, 18 regions/samples (Normal, Tumor Core, Tumor Border). Core and Border are **not** two independent patients.
- Validation patients belong to a separate cohort from the **same underlying published study**, not an entirely separate independent publication.
- Candidate inclusion masks are **provisional**. This report does not convert them to final curated labels.
- Frozen UCell signatures checked against both manifests: 29 program-gene pairs; max_rank=1500 in the saved manifests.

## Discovery: cross-program associations

- **SPP1_TAM__FAP_CAF**: Spearman rho = 0.400; permutation p = 0.223; BH-adjusted p = 0.669; n = 11.
- **SPP1_TAM__CD8_Dysfunction**: Spearman rho = 0.082; permutation p = 0.817; BH-adjusted p = 0.817; n = 11.
- **FAP_CAF__CD8_Dysfunction**: Spearman rho = -0.182; permutation p = 0.596; BH-adjusted p = 0.817; n = 11.

**Interpretation:** No prespecified association was significant after BH correction. SPP1-TAM and FAP-CAF showed a positive exploratory trend but no validated patient-level association.

## Validation: within-patient tissue comparisons

Primary selection, at least 20 cells in *each* matched compared group; Core+Border is one within-patient aggregate:

- **SPP1_TAM**: tumor higher in 4/4 matched patients; median tumor-normal UCell delta = 0.344.
- **FAP_CAF**: tumor higher in 6/6 matched patients; median tumor-normal UCell delta = 0.228.
- **CD8_Dysfunction**: tumor higher in 3/3 matched patients; median tumor-normal UCell delta = 0.027.
- **CD8_Cytotoxicity**: tumor higher in 2/3 matched patients; median tumor-normal UCell delta = 0.024.

**Interpretation:** Increased SPP1-TAM and FAP-CAF scores in tumors were directionally consistent among covered matched patients, but this does NOT establish cross-program correlation, cell-cell interaction or causation. The CD8 result is particularly sparse and sensitive to cluster selection.

## Sensitivity and limitations

- Cross-lineage exclusion preserves the positive SPP1/FAP paired direction among eligible patients, but changes the sample composition and estimates.
- Exclusion of high-risk clusters leaves just 2 matched patients for SPP1 in combined Core+Border at >=20 cells, reducing interpretability.
- Complete three-program tumor coverage in GSE144735 is only 2 patients for Core and 3 for Core+Border at >=20 cells. **No confirmatory three-program correlation test is supportable.**
- Tumor-minus-normal paired estimates for different programs use differing patients and cell numbers, so magnitude comparisons between programs should not be treated as standardized effect sizes.
- Differential expression, ligand-receptor communication, and causal mechanisms have not been established by these results.

## Next reproducibility checks

- Freeze final lineage inclusion/exclusion rules with marker evidence and explicit handling of cross-lineage, donor-dominated and rare clusters.
- Save sample/patient mapping and confirm no repeated sampling is treated as independent subjects.
- Keep raw counts and large RDS outside public GitHub; publish small derived tables, scripts, dependency versions, and a limitations section.
- Report no significant cross-program replication claim unless supported by a sufficiently powered independent cohort.

