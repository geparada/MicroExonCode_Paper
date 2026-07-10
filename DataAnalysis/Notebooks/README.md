# Data-analysis notebooks

This directory contains the Jupyter notebooks used for the data analysis in **Parada, Bretschneider, et al., "An expanded repertoire of brain microexons is directly impacted by autism-associated genetic variation."** The notebooks cover microexon discovery and characterization, development of the microexon code model, its interpretation and validation, and the analysis of autism-associated genetic variation.

This guide groups the notebooks by the part of the study they support and links them to the corresponding figures, so that a reader interested in a particular result can go directly to the source code behind it.

## How the notebooks fit into the wider repository

The notebooks sit on top of the pipelines and outputs provided elsewhere in this repository:

- **[`microexon-code/`](../../microexon-code)** — the microexon code model (gradient-boosted trees over splice-site, branch-point, RBP and k-mer features), including training, inference, and *in silico* saturation mutagenesis.
- **[`CNN_predictions_workflow/`](../../CNN_predictions_workflow)** — the Pangolin and DeltaSplice variant-scoring workflow used for benchmarking and complementary predictions.
- **[`OutputTables/`](../../OutputTables)** — model predictions for all microexon SNVs and ASD indels evaluated in the study.
- **[`Input/`](../Input)** — input tables used by several notebooks. Environment setup is described in the [parent README](../README.md).

Most notebooks read the model outputs and annotation tables produced by these pipelines and turn them into the statistics and figure panels of the paper. Inputs derived from the SSC, SPARK and MSSNG cohorts are controlled-access (SFARI Base and Autism Speaks MSSNG) and are referenced as `CONTROLLED_ACCESS/…` placeholders; access is governed by the respective data-use agreements.

---

## Microexon discovery and characterization
*Figure 1, Figure 2, Extended Data Figs. 1–5*

| Notebook | Contents |
|---|---|
| `human_stats.ipynb` | Human microexon discovery and characterization: splice-site and expression features, NMD proportions, PSI clustering across bulk and single-nucleus RNA-seq, and detection-saturation analysis (Figure 1; Extended Data Fig. 1). |
| `mouse_microexons.ipynb` | Mouse microexon discovery and quantitative characterization across bulk, TRAP and single-cell datasets. |
| `human_mouse_stats.ipynb` | Cross-species comparison of 5′/3′ splice-site strength and conservation for new vs. known vs. putative microexons (Figure 1b–d). |
| `Saturation Curve.ipynb` | Microexon detection as a function of sequencing depth (Extended Data Fig. 1a). |
| `Human novel microexons blacklist.ipynb` | Filtering of the newly detected human microexon set. |

## Nonsense-mediated decay (NMD)
*Figure 1f,i,j; Extended Data Figs. 2 and 4*

| Notebook | Contents |
|---|---|
| `NMD analysis.ipynb` | Prediction of NMD-inducing microexons in human and mouse from reading-frame consequences. |
| `NDM_II.ipynb` | Cycloheximide (CHX) validation of NMD-inducing microexons and the associated bar plot (Figure 1i,j; Extended Data Fig. 2). |
| `Rebutal_II_NMD.ipynb` | Correlation between CHX-responsive microexons and host-gene expression. |
| `heatmaps_NMD_not_NMD.ipynb` | PSI heatmaps contrasting NMD-inducing and frame-preserving microexons (Figure 1g; Extended Data Fig. 1i). |
| `Selecting_NMD_microexons_validation.ipynb` | Selection of candidate NMD-inducing microexons for experimental validation. |

## Cell-type and developmental regulation
*Figure 2a–d; Extended Data Figs. 3 and 4*

| Notebook | Contents |
|---|---|
| `celltype_specific_microexons_mouse_2025.ipynb` | Identification of neuronal-subtype-specific microexons (Figure 2a; Extended Data Fig. 3a). |
| `Mouse_microexon_browser.ipynb` | Developmental splicing trajectories and per-microexon browser plots (Figure 2b–d). |

## Human–mouse conservation

| Notebook | Contents |
|---|---|
| `Human_mouse_conservation_III.ipynb` | Conservation of new microexons between human and mouse by pairwise sequence alignment, including the fraction with 1:1 orthologs. |

## Functional enrichment and protein context
*Figure 2e; Extended Data Fig. 5; protein schematics in Figures 2, 4, 5*

| Notebook | Contents |
|---|---|
| `Enrichment_analyses_2025.ipynb` | Gene Ontology and Reactome functional-enrichment analysis of microexon host genes (Figure 2e; Extended Data Fig. 5). |
| `Domain analysis.ipynb` | Overlap of microexons with UniProt protein domains and intrinsically disordered regions (Extended Data Fig. 1f). |
| `Protein_draw.ipynb` | Protein ORF and domain schematics for individual microexon-containing genes. |

## Microexon code model: performance and interpretation
*Figure 3a,b; Extended Data Figs. 6 and 7*

| Notebook | Contents |
|---|---|
| `Information_leak_checkup.ipynb`, `Information_leak_checkup_II.ipynb` | Precision–recall/AUPRC performance and information-leakage controls (held-out chromosomes, SVM-BP and ortholog exclusions) (Extended Data Fig. 6d,e). |
| `Rebuttal_stats.ipynb` | SHAP feature-importance summaries, code-score comparisons across microexon classes and tissues, the ISM–conservation relationship, and definition of the "code-selected" microexon set (Extended Data Fig. 6b,c,f,g; Figure 3b). |
| `Novel microexons ISM.ipynb` | PCA and clustering of *in silico* mutagenesis profiles with silhouette scoring to define code-selected microexons (Extended Data Fig. 7c–f). |

## In silico mutagenesis and regulatory motifs
*Figure 3b–d; Extended Data Figs. 7 and 9*

| Notebook | Contents |
|---|---|
| `ISM hexamer analysis - known & novel microexons .ipynb` | Hexamer-resolved ISM impact heatmap clustered by motif similarity (Figure 3c). |
| `ISM motifs.ipynb` | Motif-cluster ISM analysis and positional profiles of regulatory hexamers, including QKI-related sites. |
| `Alignment logo.ipynb` | Sequence logos for regulatory motifs highlighted by the model (SRRM4, SRSF11, RBFOX, PTBP1, QKI) (Figure 3d). |
| `SIM_hexamer_heatmap.ipynb` | Hexamer heatmaps and Δ-score benchmarking of ISM outputs. |
| `Hexamer_ISM_extraction.ipynb` | Extraction of hexamer-level ISM scores and matched conservation. |
| `SSC_SPARK_var_hexamer_analysis.ipynb` | Hexamer/motif context of variants, including positional QKI-motif effects (Extended Data Fig. 9). |
| `Mouse ISM.ipynb` | Transformation of mouse ISM scores. |

## Model validation and benchmarking
*Extended Data Fig. 8*

| Notebook | Contents |
|---|---|
| `Variant_evaluation_Benchmark.ipynb` | Correlation of minigene ΔPSI measurements with code-model and Pangolin predictions (Extended Data Fig. 8a). |
| `Reporter PSI.ipynb` | Analysis of minigene reporter PSI, including validation constructs. |
| `Minigene sequences.ipynb` | Sequences of minigene reporter constructs. |
| `MPRA_analyses_II.ipynb` | Correlation of code scores with massively parallel splicing-assay (MaPSy) measurements across SRRM4 levels, and calibration of Δ code score to a ≥5% ΔPSI threshold (Extended Data Fig. 8c,d). |
| `PsychENCODE.ipynb`, `PsiENCODE_metadata.ipynb` | *In vivo* variant-dependent ΔPSI from matched PsychENCODE brain WGS/RNA-seq (Extended Data Fig. 8b). |
| `Pangolin_ISM_analysis.ipynb` | Per-tissue Pangolin ISM metaplots and conservation comparisons. |
| `ISM_DeltaSplice_Pangolin_code_predictions.ipynb` | Alignment of DeltaSplice/Pangolin predictions to microexon regions for comparison with the code model. |
| `Checking_pangolin_DeltaSplice_training.ipynb` | Assessment of overlap between the CNN models' training data and the microexon set. |

## Autism-associated genetic variation
*Figure 4a–e; Extended Data Fig. 10*

| Notebook | Contents |
|---|---|
| `Whole genome score TOTAL - known + novel.ipynb` | Core computation of sibling-pair Δ code scores from individualized genomes (Figure 4a,b). |
| `Minimal_sib_delta_top_delta_cm.ipynb` | Sibling-pair Δ code-score distributions and top-decile selection (Figure 4b). |
| `Variant plots III.ipynb` | Rare vs. common variant effects, ASD-gene enrichment, and candidate-variant tables (Figure 4c; Extended Data Fig. 10d,e). |
| `Integrating ASD SSC and SPARK.ipynb` | Integration of SSC and SPARK variants with conservation and population-frequency annotation (Figure 4d). |
| `Overlap SSC_SPARK with MSSNG.ipynb` | Three-cohort overlap of recurrently impacted microexons with permutation-based significance (Figure 4e; Extended Data Fig. 10g). |
| `Whole genome scores MSSNG - known + novel.ipynb` | MSSNG cohort analysis and cross-cohort Venn diagrams (Figure 4e). |
| `MSSNG variant analysis.ipynb`, `Variant characterization.ipynb` | Characterization of individual variants across cohorts. |

## Cumulative burden and ASD-status association
*Figure 4f–h; Extended Data Fig. 10*

| Notebook | Contents |
|---|---|
| `DeltaSplice_Pangolin_ASD_IV.ipynb` | Cumulative-impact analysis across models and thresholds, common-vs-rare model comparison, and conditional logistic regression of ASD status (Figure 4f; Extended Data Fig. 10f). |
| `RV_GDT.ipynb` | Preparation of pedigree/genotype inputs and Rare-Variant Generalized Disequilibrium Test (Figure 4g,h). |
| `RD_GDT_final_plots.ipynb` | RV-GDT threshold grids and final association plots (Figure 4g,h). |
| `Rebuttal_II_carrier_stats.ipynb` | Carrier statistics for predicted impactful variants across cohorts (Extended Data Fig. 10h). |

## Variant tables and resources
*Supplementary Tables; Extended Data Fig. 11*

| Notebook | Contents |
|---|---|
| `Wrapping_up_variant_tables.ipynb` | Assembly of the per-variant prediction tables and supplementary variant tables. |
| `Exome.ipynb` | Coverage of known and new microexons by commercial exome-capture panels (Extended Data Fig. 11). |

---

The protein–protein interaction network (Figure 5) and the functional-enrichment networks (Figure 2e, Extended Data Fig. 5) were assembled in Cytoscape (with STRING and EnrichmentMap) from gene and term lists generated by the notebooks above.

The directory also contains additional development notebooks accumulated over the course of the project.
