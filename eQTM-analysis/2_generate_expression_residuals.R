### Get expression residuals for eQTM analysis ###
library(readr)
library(dplyr)
library(limma)

# Gene expression matrix for our BRCA cohort
expression <- read_csv("expression_brca_entrez_ids.csv")

# Preprocess and format expression data
expression <- as.data.frame(expression)
expression <- aggregate(expression[,-1], list(expression$gene_id),data=expression,FUN=mean) # for genes with duplicated values, get the mean value

row.names(expression) <- expression[,1]
expression <- expression[,-1]

# Load demographics
FinalFullTCGAmeta <- read_csv("FinalFullTCGAmeta2.csv") # Load demographics for your samples
colnames(FinalFullTCGAmeta)[3] <- "sample"

# Keep only BRCA samples
BRCA_demographics <- FinalFullTCGAmeta[which(FinalFullTCGAmeta$DISEASE == "BRCA"),]

# Select only demographics of interest 
sub <- unique(BRCA_demographics[,c("patient","purity","consensus_ancestry")])

# Load admixture proportions
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

final <- df_mean[intersect(row.names(df_mean),colnames(expression)),] # only keep demographics for samples with expression data available

## Get expression residuals with limma ##
mval <- as.matrix(expression)
mval <- mval[,row.names(final)] 
all.equal(colnames(mval),row.names(final)) # sanity check for sample order across the two files

# Design limma model with admixture proportions and purity
design_full <- model.matrix(~`Admixture % AFR`+`Admixture % AMR`+`Admixture % EAS`+`Admixture % EUR` +`Admixture % SAS` +purity, data=final)
design_full

# Fit the model
fit <- lmFit(mval, design_full)

# Get expression residuals (adjusted expression values for selected covariates)
ResidualsMatrix <- residuals(fit, mval)

# Save expression residuals
write.csv(ResidualsMatrix,"residuals_expression_admix_proportions.csv",row.names=T)
