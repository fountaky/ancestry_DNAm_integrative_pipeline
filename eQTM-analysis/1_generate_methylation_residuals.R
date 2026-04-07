### Get methylation residuals for eQTM analysis ###
library(readr)
library(dplyr)
library(limma)

# Load methylation beta values for aDMS for our cohort
methylation <- read_csv("BRCA_betas_selected_probes.csv")

# Preprocess and format expression data
colnames(methylation) <- c("probe",substr(colnames(methylation[,-1]), start=1, stop = 12))

# Create appropriate format for beta matrix
methylation <- as.data.frame(methylation)
row.names(methylation) <- methylation[,1]
methylation <- methylation[,-1]

# Load demographics
FinalFullTCGAmeta <- read_csv("FinalFullTCGAmeta2.csv") 
colnames(FinalFullTCGAmeta)[3] <- "sample"

# Keep only BRCA samples
BRCA_demographics <- FinalFullTCGAmeta[which(FinalFullTCGAmeta$DISEASE == "BRCA"),]

# Select only demographics of interest 
sub <- unique(BRCA_demographics[,c("patient","purity","consensus_ancestry")])

# Load genetic ancestry admixture proportions
tcga_gdan_aim_ancestry_calls <- read_excel("tcga_gdan_aim_ancestry_calls.xlsx")
tcga_ancestry <- tcga_gdan_aim_ancestry_calls

# Subset to columns representing continental admixture proportions
suben <- merge(sub,tcga_ancestry[,c(1,9:13)])
non_na_suben <- na.omit(suben) # remove missing values

# For samples with duplicated purity value calculate the mean purity
df_mean <- non_na_suben %>%
  group_by(patient,consensus_ancestry,`Admixture % AFR`,`Admixture % AMR`,`Admixture % EAS`,`Admixture % EUR`,`Admixture % SAS`) %>%
  summarise(purity = mean(purity, na.rm = TRUE), .groups = "drop")

df_mean <- as.data.frame(df_mean)
row.names(df_mean) <- df_mean$patient
df_mean <- df_mean[,-1]

# Only keep entries for samples with methylation data available
final <- df_mean[intersect(row.names(df_mean),colnames(methylation)),]

## Get methylation residuals with limma ##
mval <- as.matrix(methylation[,row.names(final)]) # keep samples with available covariates of interest to control for

# Design limma model with admixture proportions and purity
design_full <- model.matrix(~`Admixture % AFR`+`Admixture % AMR`+`Admixture % EAS`+`Admixture % EUR` +`Admixture % SAS` +purity, data=final)
design_full

# Fit the model
fit <- lmFit(mval, design_full)

# Get methylation residuals (adjusted methylation values for selected covariates)
ResidualsMatrix <- residuals(fit, mval)

# Save methylation residuals
write.csv(ResidualsMatrix,"residuals_methylation_admix_proportions.csv",row.names=T)
