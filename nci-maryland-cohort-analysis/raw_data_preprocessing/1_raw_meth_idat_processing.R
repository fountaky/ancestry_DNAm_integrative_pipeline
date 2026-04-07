# Load libraries 
library("minfi")
library("ENmix")
library("IlluminaHumanMethylationEPICmanifest")  
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") 
library(ggplot2)
library(R.utils)

# ============================================
# Step 1: Load Data from gzipped IDATs
# ============================================
idat_dir <- "/raw/idat_files" # directory with raw .dat files

# Create temporary directory
temp_dir <- tempfile()
dir.create(temp_dir)

# Get all .gz files
gz_files <- list.files(idat_dir, pattern = "\\.idat\\.gz$", 
                       recursive = TRUE, full.names = TRUE)

# Copy and decompress to temp location
cat("Decompressing to temporary directory...\n")
for(gz_file in gz_files) {
  # Copy to temp
  temp_gz <- file.path(temp_dir, basename(gz_file))
  file.copy(gz_file, temp_gz)
  # Decompress
  gunzip(temp_gz, remove = TRUE)
}

# Read idat files with minfi
RGSet <- read.metharray.exp(base = temp_dir, 
                            recursive = TRUE, 
                            force = TRUE)


### Perform minimum QC to run the ENmix step. QC plots were generated using the old script ###

# ============================================
# Step 2: Minimal QC
# ============================================

# Get detection p-values with minfi
detP <- detectionP(RGSet)
cat("Samples with >5% failed probes:\n")
failed_samples <- colMeans(detP > 0.01) > 0.05
print(names(which(failed_samples)))


# ============================================
# Step 3: Background Correction & Dye Bias Correction with ENmix
# ============================================
cat("\n=== Background Correction & Dye Bias Correction ===\n")
mdat <- preprocessENmix(RGSet, 
                        bgParaEst = "oob",      # Out-of-band background estimation
                        dyeCorr = "RELIC",      # Dye bias correction method
                        QCinfo = NULL,            # Remove low quality samples/probes --> QC NULL here since ENmix QC function doesn't work with minfi generated RGSet
                        nCores = 6)             # Parallel processing

cat("After preprocessENmix:", dim(mdat), "\n")

# ============================================
# Step 4: Between-Array Normalization
# ============================================
cat("\n=== Between-Array Normalization ===\n")
mdat <- norm.quantile(mdat, method = "quantile1")

cat("After normalization:", dim(mdat), "\n")

# ============================================
# Step 5: Probe-Type Bias Adjustment (RCP)
# ============================================
cat("\n=== Probe-Type Bias Adjustment (RCP) ===\n")
beta <- rcp(mdat) #, qcscore = qc

cat("Beta values generated:", dim(beta), "\n")
cat("Beta value range:", range(beta, na.rm = TRUE), "\n")

# ============================================
# Step 6: QC Filtering & Optional Imputation
# ============================================
cat("\n=== QC Filtering ===\n")

# Filter low quality data (USING MINFI QC) points and samples/probes with >5% missing
beta <- qcfilter(beta,
                 detectionP = detP,
                 rthre = 0.05,        # Remove probes with >5% failed samples
                 cthre = 0.2,        # Remove samples with >20% failed probes
                 rmcr = TRUE,         # Remove cross-reactive probes
                 impute = FALSE)      # I DON'T WANT IMPUTATION

cat("After QC filtering:", dim(beta), "\n")

# Save ENmix generated betas as csv files
write.csv(beta, 
          "/processed_methylation/maryland_enmix_betas_20%_SM.csv")


