# Microexon Code Outputs

## Description

This repository provides data compiling code predictions for all Single Nucleotide Variants (SNVs) and insertions/deletions (indels) evaluated in our work, facilitating further exploration of regulatory mechanisms and genotype-phenotype relationships involving microexons.

## Data Files

### known_and_novel_microexons.known_microexon_model.hg38.csv

This file contains all the metadata for the known and new human neuronal microexons analyzed.

#### Column Description

- `event`: Microexon ID. Known microexons IDs were taken from VastDB.
- `chrom`: Chromosome.
- `strand`: Strand (+/-).
- `upIntStart`: Upstream intron start.
- `upIntEnd`: Upstream intron end.
- `dnIntStart`: Downstream intron start.
- `dnIntEnd`: Downstream intron end.
- `geneName`: Gene name.
- `lengthDiff`: Microexon length.
- `label`: Microexon label; with `hi` indicating known neuronal microexons and `novel` indicating new ones.
- `group`: Group of microexons (known/novel).
- `wt_neural_high_pred`: Wild-type prediction scores.
- `ME`: MicroExonator ID provided for novel microexons.
- `geneID`: ENSEMBL gene ID.
- `mes3`: MaxEntScan 3' splice site score.
- `mes5`: MaxEntScan 5' splice site score.
- `conservation`: phastCons conservation score.
- `neural_high_pred`: Wild-type prediction scores.

### ASD_indels_and_ISM_SNVs.deta_code_preditions.tsv

This file contains our model predictions for all possible SNVs across microexons and indels associated with Autism Spectrum Disorder (ASD) families analyzed in our study.

#### Column Description

- `variant`: Variant ID.
- `event`: Microexon ID.
- `region`: Microexon region impacted by the variant, with up/dn corresponding to upstream and downstream intronic regions respectively; ex microexon region; c1/c2 flanking exon regions.
- `rel_pos`: Distance from splice site microexon splice site. If region=c1/c2, then this corresponds to the distance to flanking exon splice sites.
- `var_type`: Variant type; indel or single nucleotide variant (SNV).
- `raw_delta`: Direct difference between wild-type and mutant prediction scores.
- `delta_logit_score`: Difference between wild-type and mutant logit prediction scores.
- `transformed_score`: Scaled delta logit scores used for all the analyses of our paper.
