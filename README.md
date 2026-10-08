# Carotid injury and ferroptosis associated transcription

Retained R analysis scripts accompanying an exploratory study of rat carotid balloon injury, single-cell transcription and bulk transcriptomics.

## Study inputs and scope

- GSE174098: single-cell evidence used for the displayed atlas and score-associated analyses.
- GSE164050: bulk evidence used for the displayed differential-expression and enrichment analyses.
- Some retained branches refer to GSE126627 or unrelated template analyses. They are included as historical variants, not validated additional cohorts.

## Repository layout

`scripts/single_cell` contains input, QC, annotation, VSMC/score, enrichment, GSVA and communication script variants.

`scripts/bulk_transcriptome` contains preprocessing, limma, heatmap/volcano, enrichment, RF, alternative SVM/LASSO, intersections and exploratory validation/clustering scripts.

`docs/methods_and_code_audit.md` records parameters recovered from the scripts and the limits on their interpretation. `docs/script_inventory.json` maps each published script to its original filename and SHA256. Published copies use UTF-8; private absolute local paths are replaced by LOCAL_PROJECT_PATH. Source originals are retained separately by the authors.

## Reproduction status

This repository provides the retained source code, not a certified end-to-end reproduction. Scripts were inspected but not executed for this revision. Required raw/processed matrices, Seurat objects, final gene lists, session versions and some external helper scripts are absent. Script numbering indicates a historical workflow, not a verified execution order. Do not source all scripts sequentially: several clear the workspace, require objects from earlier interactive sessions, contain alternate branches or exploratory plotting fragments.

Before running, reconcile sample metadata and animal pairing, restore the required inputs, identify the figure-generating version, and record package versions. No lockfile or sessionInfo is fabricated. packages_detected.txt lists statically detected packages only.

## Important interpretation checks

- Rat mitochondrial/hemoglobin gene matching needs correction and recalculation; archived zero-valued QC panels are unverified.
- Clustering variants differ in variable-feature counts and resolution. No integration call is evident in the inspected branches.
- The final FerrScore execution and conflicting plotted scales are not established; symbol capitalization is not ortholog mapping.
- Bulk limma scripts omit animal blocking despite paired sampling; thresholds differ among branches.
- 05-05venn韦恩图.R explicitly adds Plin2, Atg13 and Cdo1 to three lists. This overlap cannot establish independent algorithmic discovery.
- Human Hallmark/CellChat references require validated rat-to-human mapping. Immune deconvolution and tissue image generation are not fully covered by the supplied scripts.

Do not interpret a repository publication as confirmation that these issues were corrected. Figures were not recalculated. The code archive preserves historical evidence for transparent reporting.

## Data availability

Public accession pages: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174098 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164050 . Data and manuscript files are not hosted here.

## Attribution and licensing

Existing script comments and third-party helper references are retained. No blanket license is assigned because ownership and licensing of embedded/referenced third-party fragments have not been verified.
