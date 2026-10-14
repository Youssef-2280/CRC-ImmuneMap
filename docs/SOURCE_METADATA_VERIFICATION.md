# Original data and sample-pairing verification (2026-10-14)

## Evidence sources

- **Original published study:** Lee H-O, Hong Y, Etlioglu HE, et al. *Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer*. Nature Genetics **52**, 594–603 (2020). DOI: [10.1038/s41588-020-0636-z](https://doi.org/10.1038/s41588-020-0636-z), PubMed [PMID 32451460](https://pubmed.ncbi.nlm.nih.gov/32451460/).
- **Discovery source:** GEO [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465), 23 Korean CRC patients, 23 primary tumor and 10 matched normal mucosa samples; 63,689 author-QC-retained cells.
- **Validation source:** GEO [GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735), 6 Belgian CRC patients, 6 each tumor core, tumor border and matched normal mucosa; 27,414 author-QC-retained cells.
- The two series are **different patient cohorts within the same Nature Genetics publication**, not independently published datasets. KUL3 is a distinct patient/cohort check but not a separate-study, confirmatory external validation.

## Sample-to-patient checks completed

1. Compiled the **51 GEO sample IDs and GSM accessions** in `data/derived/GEO_verified_51_sample_accession_map.csv`: SMC 33 (23 T, 10 N); KUL3 18 (6 T, 6 B, 6 N). Accession IDs were transcribed from the relevant **GEO series sample listings**.
2. Compared all **33 SMC sample/patient labels** to the project's archived `data/derived/sample_mapping_check.csv`: one patient/condition per sample; 33/33 labels agree with GEO; normal suffix N only in SMC01–SMC10.
3. Compared the **18 KUL3 sample–patient–region combinations** represented in the user's Step 11B uploaded patient/sample/subtype evidence tables against GEO sample labels: all 18/18 combinations agree; patient ID matches KUL prefix; T = Tumor_Core, B = Tumor_Border, N = Normal.
4. The original Step 07B import code explicitly checked all 27,414-cell annotation `Patient`, `Sample` and `Class` columns against these suffixes, and the user completed the downstream analysis. That is corroborating computational evidence; the original 27,414-row author annotation file was **not independently re-read in this review session**.
5. Patient is the statistical unit. Tumor Core and Border are repeated within-patient regions, not two biological replicates. Core+Border is one aggregate per patient.

**Verification boundary:** This establishes agreement of public GEO sample IDs, project summaries, and author-derived sample class conventions. It **does not** independently review all 51 original sample-specific clinical records, clinical stages, molecular subtypes or paired tissue procurement logs. No molecular pathology or histological match beyond GEO's stated 'matched normal' design is claimed.

## Exact GEO annotations

- KUL3 suffix `-T` is *Core region*, `-B` is *Border region*, `-N` is *Normal mucosa*. Specific GEO samples [GSM4294474](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4294474) (KUL01-T) and [GSM4294485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4294485) (KUL31-B) separately corroborate the tissue labels.
- Do not call KUL3 **18 independent patients** or call SMC normals matched for all 23 subjects.
- Source count matrices in GEO are **author-processed, QC-retained raw UMI count values**; raw FASTQ availability differs by cohort (SMC controlled-access EGA, KUL3 ArrayExpress). Do not describe them as raw, unfiltered droplet libraries.
