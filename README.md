# ancestry_DNAm_integrative_pipeline

## Introduction
This repository provides code for integrative analysis of DNA methylation (DNAm), SNP genotyping array, and RNA-seq data to investigate the effects of genetic ancestry on DNAm and its regulatory effects on the transcriptome. Using preprocessed data from the TCGA breast cancer cohort, we present a pipeline for identifying ancestry-specific methylation patterns, mapping methylation quantitative trait loci (meQTLs), and further uncovering expression quantitative trait methylation (eQTMs) associated with these patterns.

The preprint associated with this pipeline can be found [here](https://doi.org/10.1101/2024.08.29.610316)

## Contents

### differential-methylation-analysis
Includes code to perform differentially DNA methylation analysis

### differential-expression-analysis
Includes code to perform differentially expression analysis

### meQTL-analysis
Includes code to obtain methylation redisuals, process them, and perform meQTL analysis

### eQTM-analysis
Includes code to identify candidate eQTMs (map CpG sites to local genes), run eQTM analysis, identify eQTM communities, and perform Gene Ontology analysis on the eQTM target genes 

### paper-figures
Includes code to generate main figures

### extended-figures
Includes code to generate extended figures
