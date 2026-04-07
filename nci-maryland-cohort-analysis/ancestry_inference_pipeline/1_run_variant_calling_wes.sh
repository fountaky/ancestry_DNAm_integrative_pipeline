#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --qos=cpu_snice  
#SBATCH --job-name=WES_sarek
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00 
#SBATCH --output=WES_sarek_%A_%a.out

conda activate nf_latest

# Set up paths
INPUT="sample_sheet_wes_sarek.csv" # .csv file with paths to fastq files. Headers shoud be: sample,fastq_1,fastq_2,status,lane,patient
OUTPUT="/nci_maryland/outpath/"
GENOME="Homo_sapiens_assembly38.fasta"
FASTA_DICT="Homo_sapiens_assembly38.dict"
FASTA_FAI="Homo_sapiens_assembly38.fasta.fai"
DBNSFP="dbNSFP5.3a_grch38.vcf.gz"
DBNSFP_TBI="dbNSFP5.3a_grch38.vcf.gz.tbi"
DBSNP="Homo_sapiens_assembly38.dbsnp138.vcf.gz"
DBSNP_TBI="Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"

# ------------------------
# Run nf-core/sarek
# ------------------------
nextflow run nf-core/sarek \
   -profile singularity \
   --input "$INPUT" \
   --outdir "$OUTPUT" \
   --wes \
   --tools haplotypecaller,vep \
   --trim_fastq \
   --aligner bwa-mem2 \
   --saved_mapped \
   --use_gatk_spark false \
   --vep_dbnsfp \
   --dbnsfp "$DBNSFP" \
   --dbnsfp_tbi "$DBNSFP_TBI" \
   --fasta "$GENOME" \
   --dict "$FASTA_DICT" \
   --fasta_fai "$FASTA_FAI" \
   --dbsnp "$DBSNP" \
   --dbsnp_tbi "$DBSNP_TBI" \
   -resume



   