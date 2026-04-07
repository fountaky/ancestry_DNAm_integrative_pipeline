#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --qos=cpu_snice  
#SBATCH --job-name=ADMIXTURE_Run
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00 
#SBATCH --output=ADMIXTURE_%A_%a.out

# Load necessary modules
module load EB5
module load EB5Modules
module load EBModules
module load BCFtools/1.22-GCC-14.2.0
module load plink2/2.0.0
module load ADMIXTURE/1.3.0

# Define the array of pop groups / Will be running admixture separately for African American and White samples
declare -a POP_ARRAY=("AA" "White")

# Get the current population from the array index
POP=${POP_ARRAY[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing population: $POP"

# -----------------------------
# Perform admixture on merged dataset
# -----------------------------
INPUT_PATH="/pop_multi_sample_WES_vcf_qc/${POP}_vcf"

# Range of K values to run ADMIXTURE for HOW DO YOU DECIDE ON THE K?
START_K=2
END_K=10

# Number of threads to use
THREADS=8

# Optional: output directory
ADMIXTURE_OUT="/pop_multi_sample_WES_vcf_qc/${POP}_vcf/admixture"
mkdir -p "$ADMIXTURE_OUT"

# Change to output directory to avoid file conflicts
cd "$ADMIXTURE_OUT"

for K in $(seq $START_K $END_K); do
  echo "Running ADMIXTURE for K=${K} with ${THREADS} threads..."

  # Run ADMIXTURE with multithreading enabled
  admixture --cv -j${THREADS} "${INPUT_PATH}/merged.WES_WGS_pruned_r2_0.2.bed" $K | tee "$ADMIXTURE_OUT/log${K}.out"

done

