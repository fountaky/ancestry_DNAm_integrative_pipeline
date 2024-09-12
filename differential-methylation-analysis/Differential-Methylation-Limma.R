### BRCA analysis to identify ancestry-associated probes
library(readr)
library(limma)
library(dplyr)
library(stringr)
library(data.table)

# Load demographics for all TCGA samples
FinalFullTCGAmeta <- read_csv("FinalFullTCGAmeta2.csv") # Load demographics for TCGA samples
colnames(FinalFullTCGAmeta)[3] <- "sample"

# Remove samples with NA age
FinalFullTCGAmeta <- FinalFullTCGAmeta[-which(is.na(FinalFullTCGAmeta$age_at_initial_pathologic_diagnosis)),] 
dim(FinalFullTCGAmeta)
demographics <- FinalFullTCGAmeta

# Subset demographics dataframe to only BRCA samples
BRCA_demographics <- demographics[which(demographics$DISEASE == "BRCA"),]
BRCA_demographics <- na.omit(BRCA_demographics) 

# Preprocess BRCA demographics
BRCA_demographics <- BRCA_demographics[,-c(1:2)]
BRCA_demographics <- data.frame(BRCA_demographics)
row.names(BRCA_demographics) <- as.character(BRCA_demographics[,1])
BRCA_demographics <- BRCA_demographics[,-1]

# Kepp only AFR and EUR samples
BRCA_demographics <- BRCA_demographics[BRCA_demographics$consensus_ancestry == "eur" | BRCA_demographics$consensus_ancestry == "afr",]

# Load methylation matrix for the BRCA cohort
dat <- fread("BRCA_methylation.csv")
methylation <- as.data.frame(dat)
methylation <- methylation[,-1]
colnames(methylation)[1] <- "probe"

# Subset the methylation matrix - Keep samples with demographic data available
patients <- which(colnames(methylation) %in% row.names(BRCA_demographics))
methylation <- methylation[,c(1 , patients)]

# Use limma to fit the linear regression models
rownames(methylation) <- as.character(methylation[,1])
mval <- methylation[,-1]
mval <- mval[,order(colnames(mval))]

# Make sure there is the same order of samples between demographics and methylation array have
BRCA_demographics <- BRCA_demographics[order(row.names(BRCA_demographics)),]
all.equal(colnames(mval),row.names(BRCA_demographics))

# Ancestry and subtype information as factors
BRCA_demographics$consensus_ancestry <- as.factor(BRCA_demographics$consensus_ancestry)
BRCA_demographics$SUBTYPE <- as.factor(BRCA_demographics$SUBTYPE)

##### Create the model ####
design_full <- model.matrix(~0 + consensus_ancestry+purity+age_at_initial_pathologic_diagnosis + SUBTYPE, data=BRCA_demographics)

fit <- lmFit(mval, design_full) # Fit the model

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts("consensus_ancestryeur-consensus_ancestryafr",
                            levels = design_full)

contMatrix

# Further model fitting steps
fit2 <- contrasts.fit(fit, contMatrix) # Fit contrasts of interest
fit2 <- eBayes(fit2) # Bayesrian estimation on original fit

# Get statistics across all probes 
DMPs_subtype <- topTable(fit2, num=Inf, coef=1,  adjust.method = "BH") # If coef=NULL 38,017 significant probes
head(DMPs_subtype)
sig.probes <- DMPs_subtype[DMPs_subtype$adj.P.Val < 0.05,] # significant pval probes

effective <- sig.probes[abs(sig.probes$logFC) >= 0.1,] # effect size filtering
dim(effective)

write.csv(DMPs_subtype, "BRCATopTable.csv")
write.csv(sig.probes, "BeforeBRCAAncestryProbes.csv")
write.csv(effective, "AfterBRCAAncestryProbes.csv")
