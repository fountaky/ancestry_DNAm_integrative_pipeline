#!/bin/bash
#$ -cwd
#$ -l m_mem_free=3G
#$ -pe threads 16

# Introduce global variable (Number of permutations to run)
n=1000 
export n

# Run meQTL control analysis with permutations
module load EBModules
module load  R/4.0.3-foss-2020a
 
## Script that subsets BRCA methylation array
Rscript 1_subset_BRCA_methylation_matrix.R

## Script to generate combined methylation df with aDMS locations (per tensorQTL requirements)
Rscript 2_process_meth_matrix_per_tensorQTL_requirements.R

## Script to preprocess methylation dataframe
Rscript 3_peer_file_format_meth_matrix.R

## Python script to run PEER and obtain methylation residuals
bash 4_run_peer.py

## Combine previously generated methylation residuals with methylation df containing aDMS locations 
Rscript 5_process_methylation_residuals.R

## Convert methylation residuals to sorted bed file (tensorQTL input)
bash 6_meth_residuals_to_sorted_bed.py

## Run meQTL analysis
python 7_run_meQTL_mapping_control_perms.py
