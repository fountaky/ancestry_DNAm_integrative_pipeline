### Generate expression residuals for eQTM analysis ###
library(readr)
library(dplyr)
library(limma)

# Gene expression file
expression <- read_csv("deseq2_processed_tumors.csv")

# Preprocess and format expression data
expression <- as.data.frame(expression)
row.names(expression) <- expression[,1]
expression <- expression[,-1]

# Load purity
purity <- read.csv("normal-mast_purity_ALL_InfPurify.csv")
colnames(purity) <- c("accession","purity")

# Load admixture proportions
ancestry <- read.csv("final_ancestry_maryland.csv")

# Combine data
merged <- merge(purity, ancestry, by = "accession")
row.names(merged) <- merged$accession

# Load expression metadata to get accession number (for further sample matching with methylation samples)
expr_metadata <- read_csv("expr_metadata.csv")
expr_metadata$title<- gsub("S_","", expr_metadata$title)

expr_metadata <- as.data.frame(expr_metadata)
row.names(expr_metadata) <- expr_metadata$title

expr_metadata$accession <- gsub("accession: ", "",expr_metadata$`accession:ch1`)

expr <- expr_metadata[colnames(expression),]

# Add accession number to expression matrix
expression2 <- expression
colnames(expression2) <- expr$accession

final <- merged[intersect(row.names(merged),colnames(expression2)),]

## Get expression residuals with limma
mval <- as.matrix(expression2)
mval <- mval[,row.names(final)]
all.equal(colnames(mval),row.names(final))

# Design model
design_full <- model.matrix(~ AFR + AMR + EAS + EUR + SAS + purity, data=final)

# Train the model
fit <- lmFit(mval, design_full)

# Get residuals
ResidualsMatrix <- residuals(fit, mval)

# Save file
write.csv(ResidualsMatrix,
          "residuals_DESeq2expression_tumors.csv",
          row.names=T)

