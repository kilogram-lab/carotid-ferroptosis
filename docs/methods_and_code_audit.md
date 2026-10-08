# Supplementary methods and retained code audit

The supplied R scripts are retained variants inspected during revision. They were not executed, and their correspondence to the final figures has not been established. Original input matrices, saved objects, gene-list CSVs, session versions and execution logs are not included. The scripts contain template fragments, implicit session dependencies and alternative branches. This archive is evidence for reporting choices and limitations, not a certified reproducible pipeline.

| Topic | Recovered setting | Interpretation boundary |
|---|---|---|
| QC | `02_01_scRNA_Mine_now.R`: `^MT-`; uppercase human-style HB list; >200 and <3000 detected genes; mt <20 and HB <3 | Gene names must be verified in the rat object. Existing zero panels were not recalculated. |
| Normalization and clustering | LogNormalize 10000; vst 6000, PCs 1:20, resolution 0.5, algorithm 1, seed 2021 in the now branch; 5000/resolution 0.2 in the earlier branch | No integration call found in inspected branches. Final branch and rationale are not verified. |
| VSMC analysis | `02_02_scRNA_Mine_now_ferr.R` subsets VSMCs and reuses the stored PC/UMAP representations; reclusters with PCs 1:20 and resolution 0.5 | Re-clustering is not evidence of a freshly recomputed VSMC embedding. Membership needs the saved object. |
| Score | Title-case symbol conversion, intersection with FindMarkers output; AddModuleScore call commented out in the VSMC script | Not an ortholog map. The final score input list, execution and scale discrepancy remain unresolved. |
| GSVA | `02_03_GSVA_Ferroptosis.R`: human Hallmark sets; AverageExpression by ferr_group; method ssgsea; row scaling | Species mapping and final matrix absent. Script clears the session but depends on existing objects. Two-group row scaling can produce opposing colors without showing absolute pathway activity. |
| CellChat | Below Q1 / above Q3 / middle; uppercase symbols; human secreted-signaling database and PPI; minimum ten cells; node size by group cell count | Uppercasing is not validated orthology; versions and final execution are missing. |
| Bulk preprocessing | normalizeBetweenArrays; distribution-triggered log2 | Final count matrix and execution absent; suitability for RNA sequencing requires reassessment. |
| Bulk limma | model.matrix(~group); Con reference, Neo coefficient; eBayes; no animal blocking | Sampling is paired but the retained design is unpaired. |
| Gene thresholds | Nominal P <0.05 or <0.03 with abs(logFC)>0.3; a heatmap branch uses adjusted P <0.05 and abs(logFC)>1 | Do not describe all lists as FDR-filtered; final lists need reconciliation. |
| Bulk ranking | Export of descending logFC from allDiff | Final enrichment engine, gene-set versions/permutations and correspondence to GSEA panels unavailable. |
| Candidate box plots | stat_compare_means method t.test, no paired argument | Nominal annotations are retained; they are not redefined as paired tests. |
| Feature ranking | classif.ranger, impurity importance; alternative SVM branches | Small internal cohort, no verified external validation. |
| Candidate overlap | `05-05venn韦恩图.R` manually adds Plin2, Atg13 and Cdo1 to ScRNADEG, RF and SVM lists | The displayed intersection must not be claimed as independent algorithmic discovery. Final source lists and execution require validation. |
| Other cohort branches | GSE126627 appears in exploratory intersection/clustering scripts | These do not establish a validated additional cohort or its role in displayed findings. |
| Immune / tissue | Complete immune generation code, animal age, antibody records, field selection and IHC/TEM test definitions not recovered | Those reporting gaps remain. |

No plot values, group labels, scores or micrographs were regenerated. Source R files were not modified. Repository publication status and the pinned commit are recorded separately in repository_status.json. The supplied unrelated PE-netMeta URL is not used as this study’s code link. Review source comments and local paths before external upload; scientific scripts retain their original provenance comments.
