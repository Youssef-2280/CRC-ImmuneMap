# CRC-ImmuneMap

**Patient-level single-cell evaluation of macrophage SPP1, fibroblast FAP, and CD8 dysfunction programs in colorectal cancer**

> **Evidence status:** exploratory. This repository contains a **user-approved final primary analysis selection** (`primary_07e`), derived summaries and presentation figures. It **does not** claim that the SPP1–FAP–CD8 interaction hypothesis is statistically validated. There is **no independently published external-study validation** here.

## At a glance

| Analysis | Evidence |
| --- | --- |
| Discovery GSE132465 | 23 patients available; 11 complete-case tumor patients for primary program associations. SPP1–FAP Spearman ρ = +0.400, permutation p = 0.223, BH p = 0.669 (**non-significant**). |
| Validation GSE144735 | 6 different patients from **same published study**; 18 regional samples, Normal/Core/Border. No independent-study claim. |
| Tumor Core+Border vs matched Normal (>=20 cells per side) | SPP1_TAM increased in 4/4 (median delta +0.344); FAP_CAF increased in 6/6 (+0.228); CD8_Dysfunction 3/3 (+0.027), all **descriptive**. |
| Testing the three-program association in validation | Not adequately supported: complete-case coverage only 2/6 Core, 3/6 combined Core+Border (>=20 cells). |
| Frozen signatures | 29 unique program–gene pairs; UCell maxRank=1500; see `data/manifests/`. |
| Integrity | Original Step 09 technical checks: 26/26 PASS; not independent biological validation. |

## Figures

1. [Figure 1 — Discovery correlations](figures/main/01_discovery_associations.png)
2. [Figure 2 — Matched normal / tumor program scores](figures/main/02_validation_scores.png)
3. [Figure 3 — Matched differences by tumor region](figures/main/03_paired_differences.png)
4. [Figure 4 — Sensitivity by selection variant](figures/main/04_sensitivity_analysis.png)
5. [Figure 5 — Three-program patient coverage](figures/main/05_patient_coverage.png)

See [the PNG-only figure index](docs/FIGURE_INDEX.md), `docs/FIGURE_LEGENDS.md` and `docs/RESULTS_EN.md`. An Egyptian Arabic summary is in `docs/RESULTS_AR.md`.

## Analysis selection and reproducibility

The project user **explicitly approved retaining the previously scored `primary_07e` masks**, with donor imbalance/cross-lineage flags and sensitivity analyses. Detailed policy: `docs/ANALYSIS_SET_APPROVAL.md`. The total primary selected pools in GSE144735 were 1,745 candidate macrophages (clusters 0/1/2/5), 4,230 candidate fibroblasts (0/2/3/4/9/11/12), and 1,264 author-labeled CD8 cells from T clusters 2/4/10 (723+381+160). These are *analytical populations*, not experimentally verified identities.

The current release was built from existing CSV outputs **without rerunning RDS-based scoring**. Exact per-cell ID masks, which only exist on the user's local R workstation in Step 07E/07F outputs, are not claimed to be embedded here. To freeze/checksum them (no UCell rerun): place this repo beside or inside the original R project and run the standalone script `scripts/00_export_and_verify_primary_07e_masks.R` from RStudio. By default it looks in `C:/Users/Hossa/Documents/CRC-ImmuneMap-project/CRC-ImmuneMap`, or supply `Sys.setenv(CRC_PROJECT_DIR="...")` first. It writes into the original project's `results/tables/final_annotation_review/primary_07e_LOCKED/`.

Recreate the five main PNG plots using Python after installing packages in `requirements-figures.txt`:

```bash
python -m pip install -r requirements-figures.txt
python scripts/01_reproduce_figures.py
```

Figure source tables are in `data/derived/`; pre-existing project R scripts (Steps 04 onward, where provided) in `scripts/pipeline_archive/`. **This is not a complete raw-to-results turnkey pipeline:** earlier stages 01–03 and large raw data/RDS are not packaged, and some archived scripts use a local Windows project path. For full rerun, acquire the public GEO series, retain original preprocessing scripts and environment, and recreate checkpoints.

## Repository structure

```text
CRC-ImmuneMap_GitHub_Release/
  README.md
  data/derived/              # small aggregate tables (CSV)
  data/manifests/            # frozen signature
  docs/                      # results, methods, approval, figure legends
  evidence/                  # supporting audit/session/review tables
  figures/main/              # 5 result PNGs (no SVG)
  figures/supplementary/     # 21 single-cell PNGs
  scripts/01_reproduce_figures.py
  scripts/00_export_and_verify_primary_07e_masks.R
  scripts/pipeline_archive/  # historic R checkpoints provided in conversation
  checksums/                 # SHA256 source/release manifest
  .gitignore
```

## Replace the existing GitHub repository files (keep the same repo URL)

The existing private `CRC-ImmuneMap` repository can be updated using GitHub Desktop. Clone the existing repository locally, replace its working-tree contents with this unpacked v5 release, then **Commit to main** and **Push origin**. Removed old SVGs and renamed paths are recorded as deletions; the new `figures/main/` and `figures/supplementary/` are additions. Review the change list before committing. Do not delete the repository or upload the ZIP itself. See [Arabic replacement instructions](docs/GITHUB_REPLACE_GUIDE_AR.md).

## Verified source study and manuscript-ready documentation (2026-10-14)

Original source and method reference: **Lee et al. (2020), Nature Genetics 52:594–603**, DOI [10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z), PMID 32451460. UCell method: **Andreatta & Carmona (2021)**, DOI [10.1016/j.csbj.2021.06.043](https://doi.org/10.1016/j.csbj.2021.06.043).

- [Source sample mapping verification](docs/SOURCE_METADATA_VERIFICATION.md) (51 public GEO sample–accession IDs cross-checked with project-derived metadata; **not full individual clinical case verification**)
- [GEO verified 51-sample accession map](data/derived/GEO_verified_51_sample_accession_map.csv)
- [QC, annotation and doublet limitations](docs/ANNOTATION_QC_DOUBLET_ASSESSMENT.md)
- [Manuscript-ready Methods / Results / Limitations draft](docs/MANUSCRIPT_READY_METHODS_LIMITATIONS.md)
- [الملخص النهائي للمراجعة العلمية](docs/FINAL_SCIENTIFIC_REVIEW_AR.md)
- [Verified primary references (BibTeX)](references/verified_references.bib)

GEO sample labels support the T/B/N pairing design in KUL3. The authors' upstream quality filtering was documented; no new independent doublet call was produced by this repository. The exact 07E/07F per-cell masks still require a **local export** if claimed as shared reproducibility artifacts.

## Publication cautions

There is **no BH-significant Discovery association** among the three hypotheses. A recurring *within-patient tumor-vs-normal elevation* is a different question from *association across patients*. The validation cohort is within the same published study and is small. Cross-lineage transcript detection can reflect several biological/technical causes; donor dominance and QC issues remain. No CellChat, ligand–receptor, spatial validation, treatment response or causal demonstration is provided.

Please see `docs/METHODS_AND_LIMITATIONS.md` and `docs/RELEASE_CHECKLIST.md` for remaining source-metadata and biological-review steps before a manuscript submission. No author names, license, or DOI has been invented.

## Native R/Seurat QC and cell-annotation figures (optional, local RDS required)

The manuscript-style study **workflow illustration generated separately is excluded**; it does not represent measured single-cell plots. To generate real QC, UMAP, broad annotation, subclustering, marker DotPlots and analytical-mask overlays **directly from original Seurat RDS checkpoints**, use `scripts/02_generate_singlecell_QC_annotation_figures.R` in RStudio:

```r
source("scripts/02_generate_singlecell_QC_annotation_figures.R")
```

The original R plotting script produces figures in `results/figures/R_singlecell_QC_annotation/` on the local workstation. **Twenty researcher-supplied images were curated into 21 distinct supplementary PNG views**, because two faint UMAPs were replaced by high-contrast R re-renders and a split tumor/normal view was added. Original source mappings are retained in `docs/`. Two are copied original pre-filtering QC plots and 18 are generated from RDS checkpoints. The RDS themselves are not packaged. The script reads checkpoints serially without reclustering or re-scoring, and creates `FIGURE_SOURCE_MANIFEST.csv` listing sources and skipped panels. Discovery post-filter QC and validation authors-processed QC are labeled as such; original pre-filter figures, if present, are copied rather than invented. Plots are supplementary candidates until inspected and incorporated into a curated figure set.


## Start here: curated PNG-only R figures (v5)

- [Arabic upload and figure guide](docs/START_HERE_R_FIGURES_AR.md)
- [Native R figures (PNG)](figures/supplementary/) — **21 curated PNG views**; [Figure index](docs/FIGURE_INDEX.md).
- The initial 20 plots were supplied by the researcher. Two pale discovery UMAPs were replaced with user-supplied R high-contrast re-renders, and a split view was included; no plotted data were invented. The images marked preliminary/unreviewed are not to be labeled as confirmed final cell types.
- The **five results figures** under `figures/main/` were rendered from R-exported numeric CSVs using **Python** (`scripts/01_reproduce_figures.py`). They remain Python-rendered from R-derived tables. This distinction matters if you want an R-only figure-generation pipeline.
- For original R figure scripts, use `scripts/02_generate_singlecell_QC_annotation_figures.R` and retained project RDS.
