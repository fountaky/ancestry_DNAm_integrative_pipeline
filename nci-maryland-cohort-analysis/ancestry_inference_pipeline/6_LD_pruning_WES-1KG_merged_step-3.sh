#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --job-name=LD_pruning_White-1KG
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=LD_pruning_%A_%a.out

# Load necessary modules
module load EB5
module load EB5Modules
module load EBModules
module load BCFtools/1.22-GCC-14.2.0
module load plink2/2.0.0

# Define the array of r² thresholds
R2_THRESHOLDS=(0.1 0.2 0.5) # run LD for various thresholds and then choose optimal one

# Get the r² threshold for this array task
# Array indices in SLURM start at 1, but array indices in bash start at 0
R2=${R2_THRESHOLDS[$SLURM_ARRAY_TASK_ID-1]}

# Define your output directory
OUTDIR="/ancestry_calling/pop_multi_sample_WES_vcf_qc/White_vcf" # Change this path to run LD prunning in the African American cohort

echo "Running LD pruning with r² threshold = $R2"

# -----------------------------
# Perform LD pruning on the combined WES - WGS dataset
# -----------------------------

plink2 \
  --bfile ${OUTDIR}/merged.WES_1KG_common_QCed \
  --indep-pairwise 200 50 $R2 \
  --out ${OUTDIR}/merged.LDprune_r2_${R2}

plink2 \
  --bfile ${OUTDIR}/merged.WES_1KG_common_QCed \
  --extract ${OUTDIR}/merged.LDprune_r2_${R2}.prune.in \
  --make-bed \
  --out ${OUTDIR}/merged.WES_WGS_pruned_r2_${R2}
  
# -----------------------------
# Perform PCA on merged LD pruned dataset
# -----------------------------

plink2 \
  --bfile ${OUTDIR}/merged.WES_WGS_pruned_r2_${R2} \
  --pca 20 approx \
  --out ${OUTDIR}/merged.WES_WGS_pca_r2_${R2}

echo "Completed LD pruning and PCA with r² threshold = $R2"
