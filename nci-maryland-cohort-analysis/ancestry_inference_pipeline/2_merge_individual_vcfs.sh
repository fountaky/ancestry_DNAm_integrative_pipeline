#!/bin/bash
#$ -cwd
#$ -l m_mem_free=16G
#$ -pe threads 16

module load EB5
module load EB5Modules
module load BCFtools/1.22-GCC-14.2.0

# Combine individual vcf files to multisample vcf
bcftools merge -Oz -f PASS \
  -o merged_maryland.vcf.gz \
  --file-list maryland_ancestry_vcfs.txt # this .txt files is a list with all the individual vcf files to be combined into one
