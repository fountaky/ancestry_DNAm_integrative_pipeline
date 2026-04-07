#!/bin/bash
#$ -cwd

# Load necessary modules
module load EB5
module load EB5Modules
module load EBModules
module load BCFtools/1.22-GCC-14.2.0
module load plink2/2.0.0

# Here we are about to subet African American and White samples from the multisample full cohort vcf and create 2 separate population-specific vcfs
WES_VCF="merged_maryland.vcf.gz"
OUTDIR="outpath"
SAMPLE_LIST_PATH="path_with_pop_samples"

# AA samples
bcftools view -S ${SAMPLE_LIST_PATH}/AA_WES_TAN_samples.txt ${WES_VCF} -O z -o ${OUTDIR}/merged_maryland_AA.vcf.gz

# White samples
bcftools view -S ${SAMPLE_LIST_PATH}/White_WES_TAN_samples.txt ${WES_VCF} -O z -o ${OUTDIR}/merged_maryland_White.vcf.gz