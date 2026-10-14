# Biological annotation, QC and doublet assessment — transparent limitations

## QC verified from the original publication/GEO

The GEO GSE144735 and GSE132465 records state that the authors performed quality filtering using UMI counts, detected genes and mitochondrial fraction; combined Seurat multiCCA and RCA for concordant initial identities; removed discordant assignments and high-feature-count outliers as candidate doublets, leaving **27,414** and **63,689** cells respectively. KUL3 GEO sample metadata (e.g. GSM4294474 and GSM4294485) specify thresholds: **UMIs >1,000; detected genes >200 and <6,000; mitochondrial fraction <20%**. These are the *authors' upstream filters*, not novel thresholds tested by CRC-ImmuneMap. The present validation analysis imported the already-filtered GEO count matrix and did not apply a second automatic QC exclusion.

## User-approved analysis selection and warnings

- Primary `primary_07e`: macrophage candidates **1,745**, fibroblast candidates **4,230**, author-subtype-gated CD8 cells **1,264**. The latter consist of 723 (cluster 2), 381 (cluster 4), and 160 (cluster 10) cells. The analysis selection was approved by the project user on **2026-10-14**; this was not blind independent pathologist or flow-cytometry validation.
- **Donor-dominated examples:** macrophage cluster 0 has **94.7% KUL01**; cluster 1 **77.4% KUL21**; CD8 candidate cluster 10 **92.8% KUL19**. Donor dominance is a sensitivity/review flag, **not an automatic cell-exclusion rule**.
- **Cross-lineage co-detection:** macrophage cluster 0: **54.7%** show the configured two-marker macrophage and two-marker epithelial co-detection; fibroblast clusters 0 and 12: **45.3%** and **56.5%**. In targeted audits median epithelial-panel UMI counts are small (myeloid 0: **3**, fibroblast 0: **2**, fibroblast 12: **3**) despite frequent transcript detection. This does **not distinguish** ambient epithelial RNA, bona fide transcription, technical multiplets, or cell-bound contamination.
- Aggregate Step 07F primary-pool cross-lineage **flags**: 511/1,745 macrophage candidates, 817/4,230 fibroblast candidates, and 168/1,264 CD8 candidates. These flags do **not** measure a validated doublet rate.
- Published CD8 labels are combined with lineage checks for selection. A whole T-cell cluster was **not** automatically treated as CD8; mixed γδ/CD8 clusters were gated at the cell level using the author-provided CD8 subtype. Mere detection of a gamma-delta transcript is not a definitive γδ lineage call.

## What was NOT verified by this review

- No new unsupervised **doublet detection** (Scrublet/DoubletFinder/scDblFinder) on original barcodes, no ambient RNA correction (e.g. SoupX/CellBender), no orthogonal protein/histology or spatial colocalization.
- The authors' upstream filtering **reduces but cannot prove elimination of all doublets**. Therefore it would be incorrect to label the existing flagged cells as `confirmed doublets` or `confirmed ambient contamination`.
- Marker-panel co-detection is not a direct cell-type purity estimate; the adopted populations are analytic candidate compartments with recorded uncertainty.
- Cross-lineage exclusion and high-risk cluster omission are **sensitivity analyses**, not proof that one specific retained mask is physically pure.

## How to report in Methods/Limitations

Report the original QC as *performed by the source authors*, and the present work as *secondary analysis of author-QC-filtered processed UMI matrices*. Report cluster selection and author-subtype-gated CD8, donor dominance and cross-lineage flags; do not claim independently performed quantitative doublet filtering or definitive cellular identities.

**Source URLs:** [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465); [GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735); [GSM4294474](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4294474); [GSM4294485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4294485).
