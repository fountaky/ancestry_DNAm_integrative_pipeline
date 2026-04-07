#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --job-name=POP_QC_step1
#SBATCH --array=1-2
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=POP_QC_step1_%A_%a.out

# Load necessary modules
module load EB5
module load EB5Modules
module load EBModules
module load BCFtools/1.22-GCC-14.2.0
module load plink2/2.0.0


# Define the array of pop groups
declare -a POP_ARRAY=("AA" "White")

# Define MAF values to use
declare -a MAF_VALUES=(0.01 0.05) # perform MAF filtering for 2 different threshold and check which results make more sense

# Get the current population from the array index
POP=${POP_ARRAY[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing population: $POP"

WES_VCF="path_to_vcf"
OUTDIR="/pop_multi_sample_WES_vcf_qc/${POP}_vcf"

THREADS=8

# Convert vcf to bcf, include only bialelic regions
bcftools view "${WES_VCF}/merged_maryland_${POP}.vcf.gz" -Oz -Ob -o "${OUTDIR}/WES.pass.bcf"
bcftools index "${OUTDIR}/WES.pass.bcf"

bcftools view -m2 -M2 -v snps "${OUTDIR}/WES.pass.bcf" -Ob -o "${OUTDIR}/WES.pass.biallelic.bcf"
bcftools index -f "$OUTDIR/WES.pass.biallelic.bcf"

plink2 --bcf "${OUTDIR}/WES.pass.biallelic.bcf" --make-bed --not-chr X --out "${OUTDIR}/WES.pass.plink" --threads 8

# -----------------------------
# STEP 1: RUN VARIOUS QC PARAMETERS TO DETERMINE THRESHOLD FOR FILTERING (OPTIONAL)
# -----------------------------

#######################################
# 1a. Sample QC (missingness, F, freq)
#######################################
echo "[$(date)] Running sample QC metrics..."

plink2 \
  --bfile "${OUTDIR}/WES.pass.plink" \
  --missing \
  --het \
  --freq \
  --out "${OUTDIR}/WES.pass.sample_qc" \
  --threads ${THREADS}

#######################################
# 1b. Variant QC (missingness, MAF, HWE)
#    Note: HWE meaningful for WGS; skip or
#    interpret carefully for WES.
#######################################
echo "[$(date)] Running variant QC metrics..."

plink2 \
  --bfile "${OUTDIR}/WES.pass.plink" \
  --missing \
  --freq \
  --hardy midp \
  --out "${OUTDIR}/WES.pass.variant_qc" \
  --threads ${THREADS}


#######################################
# 1c. Filter variants

#######################################
# Loop through each MAF threshold
for maf in "${MAF_VALUES[@]}"; do
  echo "[$(date)] Processing with MAF threshold: $maf"
  
  plink2 --bfile ${OUTDIR}/WES.pass.plink \
    --geno 0.1 \
    --maf $maf \
    --make-bed \
    --out ${OUTDIR}/WES.pass.varQC.maf${maf}

  ##### Sample and variant QC for the variant QCed WES
  echo "[$(date)] Running sample QC metrics for MAF $maf..."

  plink2 \
    --bfile "${OUTDIR}/WES.pass.varQC.maf${maf}" \
    --missing \
    --het \
    --freq \
    --out "${OUTDIR}/WES.pass.varQC.maf${maf}.sample_qc" \
    --threads ${THREADS}

  echo "[$(date)] Running variant QC metrics for MAF $maf..."

  plink2 \
    --bfile "${OUTDIR}/WES.pass.varQC.maf${maf}" \
    --missing \
    --freq \
    --hardy midp \
    --out "${OUTDIR}/WES.pass.varQC.maf${maf}.variant_qc" \
    --threads ${THREADS}
done


