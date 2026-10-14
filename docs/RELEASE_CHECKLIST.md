# Release checklist — honest scope

- [x] Step 09 technical audit in original R session: 26/26 checks PASS.
- [x] Frozen 29 program–gene entries (`maxRank = 1500`) checked.
- [x] `primary_07e` retained as the **final primary analysis set by explicit user approval**; cluster inclusion and CD8 cell-level gates documented.
- [x] Figure generation from saved aggregate tables; discovery inference and matched validation contrasts separately labeled.
- [x] Sensitivity and low patient-coverage caveats preserved.
- [x] README, results, figure legends, environment, methods and archive of available analysis scripts packaged.
- [ ] **On local Windows R project:** Run `scripts/00_export_and_verify_primary_07e_masks.R` to save and checksum exact per-cell masks from Step 07E/07F, if releasing those as processed identifiers is desired.
- [x] Verify public GEO series sample identifiers and project sample-to-patient/region conventions: 51/51 sample IDs consistent (SMC 33 + KUL3 18); see `docs/SOURCE_METADATA_VERIFICATION.md`.
- [ ] Independently validate the **full original clinical/clinicopathological metadata** (stage/MSI etc.) and, if needed, physical matched-tissue provenance; not possible from uploaded summary tables alone.
- [x] Document original GEO authors' QC, donor dominance and cross-lineage/doublet uncertainty; see `docs/ANNOTATION_QC_DOUBLET_ASSESSMENT.md`.
- [ ] Orthogonal confirmation or independent doublet analysis is **not completed**; required only if claiming definitive identities rather than candid exploratory candidate pools.
- [x] Verify original dataset study citation and UCell method citation; BibTeX saved under `references/`.
- [ ] Insert *this project's actual* authors/affiliations and decide repository license; none fabricated here.
- [x] The user created a private GitHub repository and uploaded the prior v4 files.
- [ ] Replace those old files with this v5 PNG-only release using GitHub Desktop. The assistant has not modified the remote repository.
