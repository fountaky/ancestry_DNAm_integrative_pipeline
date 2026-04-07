# Ancestry DNA Methylation Integrative Pipeline

## Introduction
This repository provides code for integrative analysis of DNA methylation (DNAm), SNP genotyping array, and RNA-seq data to investigate the effects of population-related genetic variation on DNAm and its downstream effects on the transcriptome. Using preprocessed data from the TCGA breast cancer cohort, we present a pipeline for identifying ancestry-related methylation patterns, mapping methylation quantitative trait loci (meQTLs), and further uncovering expression quantitative trait methylation (eQTMs) associated with these patterns. We also provide code for the identification of genetic ancestry effects on DNAm in the NCI-Maryland Breast Cancer cohort, using raw whole exome sequencing, DNAm and RNAseq data.

The preprint associated with this pipeline can be found [here](https://doi.org/10.1101/2024.08.29.610316). 

<img width="6000" height="4200" alt="Founta-et-al-Schematic-revision-v2" src="https://github.com/user-attachments/assets/a830f5e7-9b02-46a7-9953-6f85c6757821" />


## Contents

### differential-methylation-analysis
Includes code to perform differential DNA methylation analysis.

### differential-expression-analysis
Includes code to perform differentially expression analysis.

### meQTL-analysis  
Includes code to:  
- Obtain methylation residuals, process them, and perform meQTL analysis on aDMSs (`ancestry_meQTL_analysis/`)  
- Perform control meQTL analysis with permutations (`control_perm_meQTL_analysis/`) 
- Obtain population-specific gnomAD allele frequencies for ancestry-related and control meSNPs (`gnomad_AF_analysis/`) 

### eQTM-analysis
Includes code to identify candidate eQTMs (map CpG sites to local genes), run eQTM analysis, identify eQTM communities, and perform Gene Ontology analysis on the eQTM target genes.

### nci-maryland-cohort-analysis
Includes code used to analyze the NCI-Maryland Breast Cancer cohort data. We have provided the necessary scripts to:
- Preprocess raw methylation (.idat files) and RNAseq (counts) data (`raw_data_preprocessing/`)  
- Perform genetic ancestry inference from WES data (`ancestry_inference_pipeline/`)  
- Conduct differential methylation analysis (`differential_methylation_analysis/`) 
- Conduct eQTM and eQTM community detection analysis (`eQTM_analysis/`) 

### paper-figures
Includes code to generate main paper figures.

## Citation  
If you use this code, please cite:

Founta & Chambwe (2024). Genetic ancestry-specific meQTLs control immune function regulation in a breast cancer cohort of African and European patients. *bioRxiv*.    
https://doi.org/10.1101/2024.08.29.610316  
