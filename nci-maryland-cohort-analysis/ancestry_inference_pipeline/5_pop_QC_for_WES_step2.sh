#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --job-name=POP_QC_step2
#SBATCH --array=1-2
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=POP_QC_ste21_%A_%a.out

# Load necessary modules
module load EB5
module load EB5Modules
module load EBModules
module load BCFtools/1.22-GCC-14.2.0
module load plink2/2.0.0

# Define the array of pop groups
declare -a POP_ARRAY=("AA" "White")

# Get the current population from the array index
POP=${POP_ARRAY[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing population: $POP"

OUTDIR="/pop_multi_sample_WES_vcf_qc/${POP}_vcf"
KG_QC2="/ancestry_inference/1kgp_ref/qc_1kg/qc2" # path to reference 1000 Genomes daya   

THREADS=8

# Remove samples that seem like outliers after QC (no outliers in our case)
#plink2 \
#  --bfile ${OUTDIR}/WES.pass.varQC \
#  --remove /NCI-Maryland/processed/WES_raw_remove_samples.txt \
#  --make-bed \
#  --out ${OUTDIR}/WES.pass.varQC.samplesQC

# Export the final WES qc vcf for intersection with 1KG common data 
plink2 --bfile ${OUTDIR}/WES.pass.varQC.maf0.01 --export vcf bgz --out ${OUTDIR}/WES.pass.varQC.samplesQC
bcftools index ${OUTDIR}/WES.pass.varQC.samplesQC.vcf.gz

# Add "chr" prefix to WES file using the existing renaming file (swapping columns) (we want to match the 1KGP data chromosome format)
bcftools annotate \
  --rename-chrs <(awk '{print $2, $1}' /1kgp_ref/qc_1kg/qc2/chr_rename.txt) \
  -Oz \
  -o ${OUTDIR}/WES.pass.varQC.samplesQC.withchr.vcf.gz \
  ${OUTDIR}/WES.pass.varQC.samplesQC.vcf.gz
bcftools index ${OUTDIR}/WES.pass.varQC.samplesQC.withchr.vcf.gz

# -----------------------------
# Step 1: finding the common snps between referemce WGS and QCed WES
# -----------------------------
for CHR in {1..22}; do
  bcftools isec -n=2 \
    -p ${OUTDIR}/chr${CHR} \
    ${OUTDIR}/WES.pass.varQC.samplesQC.withchr.vcf.gz \
    ${KG_QC2}/1KG.chr${CHR}.preIntersect.bcf

  bgzip -f ${OUTDIR}/chr${CHR}/0000.vcf
  bcftools index -f ${OUTDIR}/chr${CHR}/0000.vcf.gz

  bgzip -f ${OUTDIR}/chr${CHR}/0001.vcf
  bcftools index -f ${OUTDIR}/chr${CHR}/0001.vcf.gz

  echo "[INFO] chr${CHR} common sites: $(bcftools view -H ${OUTDIR}/chr${CHR}/0000.vcf.gz | wc -l)"
done

# -----------------------------
# STEP 2: CONCATENATE ALL CHROMOSOMES
# -----------------------------
echo "[STEP 2A] Concatenating all common WES chromosomes..."
# WES intersected with 1KG
bcftools concat -Ob --threads $THREADS -o "$OUTDIR/WES.pass.qced.common.bcf" \
  "$OUTDIR"/chr*/0000.vcf.gz
bcftools index -f "$OUTDIR/WES.pass.qced.common.bcf"

echo "[STEP 2B] Concatenating all common 1KG chromosomes..."
# 1KG intersected with WES
bcftools concat -Ob --threads $THREADS -o "$OUTDIR/1KG.common.bcf" \
  "$OUTDIR"/chr*/0001.vcf.gz
bcftools index -f "$OUTDIR/1KG.common.bcf"

# -----------------------------
# STEP 3: set the variant IDs within each dataset before merging
# -----------------------------
# WES common
plink2 \
  --bcf   ${OUTDIR}/WES.pass.qced.common.bcf \
  --set-all-var-ids @:#\$r,\$a \
  --make-bed \
  --out  ${OUTDIR}/WES.pass.qced.common2

# WGS common
plink2 \
  --bcf ${OUTDIR}/1KG.common.bcf \
  --set-all-var-ids @:#\$r,\$a \
  --make-bed \
  --out ${OUTDIR}/WGS.common2

# -----------------------------
# STEP 4: Merge PLINK datasets
# -----------------------------
echo "[STEP 4] Merging WES and 1KG datasets..."

# Use plink 1.9 for the merge operation
conda activate plink1_env

plink \
  --bfile "$OUTDIR/WGS.common2" \
  --bmerge "$OUTDIR/WES.pass.qced.common2.bed" "$OUTDIR/WES.pass.qced.common2.bim" "$OUTDIR/WES.pass.qced.common2.fam" \
  --make-bed \
  --out "$OUTDIR/merged.WES_1KG_common_QCed"

conda deactivate
