### Generate methylation residuals for eQTM analysis ###
library(readr)
library(dplyr)
library(limma)

# Load methylation file
betas <- read_csv("tumors_enmix_betas.csv")

colnames(betas) <- betas[1,]
betas <- betas[-1,]

betas <- as.data.frame(betas)

row.names(betas) <- betas[,1]
betas <- betas[,-1]

# Load and process EPIC metadata file
epic <- read_csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)

epic2 <- as.data.frame(epic)

row.names(epic2) <- epic2[,1]
epic2 <- epic2[,-c(1:2)]

epic2$Methyl450_Loci[which(is.na(epic2$Methyl450_Loci))] <- "FALSE"

# Subset epic probes to those shared with 450k array 
shared <- epic2[which(epic2$Methyl450_Loci == "TRUE"),]

# Keep only 450k array shared probes
mval <- betas[which(row.names(betas) %in% row.names(shared)),]

# Load purity
purity <- read.csv("normal-mast_purity_ALL_InfPurify.csv")
View(purity)

colnames(purity) <- c("accession","purity")

# Load admixture proportions
ancestry <- read.csv("final_ancestry_maryland.csv")
View(ancestry)

# Combine data
merged <- merge(purity, ancestry, by = "accession")
row.names(merged) <- merged$accession

final <- merged[intersect(row.names(merged),colnames(mval)),]

## Get residuals with limma
mval2 <- as.matrix(mval)
mval2 <- mval2[,row.names(final)]
all.equal(colnames(mval2),row.names(final))

# Design model
design_full <- model.matrix(~ AFR + AMR + EAS + EUR + SAS + purity, data=final)

# Train the model
fit <- lmFit(mval2, design_full)

# Get residuals
ResidualsMatrix <- residuals(fit, mval2)

# Save file
write.csv(ResidualsMatrix,
          "/processed_expression/residuals_methylation_tumors.csv",
          row.names=T)
