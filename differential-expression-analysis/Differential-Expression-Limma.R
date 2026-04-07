# Load libraries
## ----------------------------------------------------------------------------
library(readr)
library(limma)
library(dplyr)
library(stringr)
library(data.table)
#' 
#' ## Load data
#' Load demographics for the whole TCGA cohort
## ----------------------------------------------------------------------------
FinalFullTCGAmeta <- read_csv("FinalFullTCGAmeta.csv") # Load metadata for TCGA samples
colnames(FinalFullTCGAmeta)[3] <- "sample"

FinalFullTCGAmeta <- FinalFullTCGAmeta[-which(is.na(FinalFullTCGAmeta$age_at_initial_pathologic_diagnosis)),] # Remove samples with NA age

#' Create dataframe with BRCA patients demographics
## ----------------------------------------------------------------------------
BRCA_demographics <- FinalFullTCGAmeta[which(FinalFullTCGAmeta$DISEASE == "BRCA"),]
BRCA_demographics <- na.omit(BRCA_demographics) 

# Keep only sample barcode
BRCA_demographics <- BRCA_demographics[,-c(2:3)]
BRCA_demographics <- data.frame(BRCA_demographics)

# Additional preprocessing
BRCA_demographics <- unique(BRCA_demographics)

# For duplicated patients, get mean average purity value
BRCA_demographics2 <- aggregate(BRCA_demographics$purity,by=list(patient=BRCA_demographics$patient,consensus_ancestry=BRCA_demographics$consensus_ancestry,age_at_initial_pathologic_diagnosis=BRCA_demographics$age_at_initial_pathologic_diagnosis,gender=BRCA_demographics$gender,race=BRCA_demographics$race,DISEASE=BRCA_demographics$DISEASE,SUBTYPE=BRCA_demographics$SUBTYPE),data=df,FUN=mean)

colnames(BRCA_demographics2)[8] <- "purity"

# Keep only individuals of Afr or Eur ancestry
BRCA_demographics2 <- BRCA_demographics2[BRCA_demographics2$consensus_ancestry == "eur" | BRCA_demographics2$consensus_ancestry == "afr",]

#' 
#' Load gene expression for the BRCA cohort
## ----------------------------------------------------------------------------
expression <- read_csv("expression_brca_tcga.csv.csv")

## ----------------------------------------------------------------------------
# Obtain batch info for every patient from their TCGA barcode
batch <- substr(colnames(expression)[-1], start = 22, stop = 25) 

# Patient id
patient <- substr(colnames(expression)[-1], start = 1, stop = 12) 

# combine patient id and batch dataframe
df <- data.frame(patient=patient,batch=batch)
#'
#' Merged df with brca demographics
## ----------------------------------------------------------------------------
demo <- merge(df,BRCA_demographics2[,c("patient","SUBTYPE","consensus_ancestry")],by="patient",all.x=T)

#'
#' Further preprocessing of the expression matrix
## ----------------------------------------------------------------------------
# Get mean expression for duplicated genes
expr <- expression
expr <- aggregate(expr[,-1], list(expr$gene_id),data=expr,FUN=mean)

row.names(expr) <- expr$Group.1
expr <- expr[,-1] # 572 BRCA samples
    
mval <- expr

#' 
#' Log2 transformation of gene expression data
## ----------------------------------------------------------------------------
# Log2 transformation of expression
trans <- log2(mval)
View(trans)
# For value == 0 replace with NA
trans[sapply(trans, is.infinite)] <- NA

#' Keep genes present in >80% of the samples
## ----------------------------------------------------------------------------
IndexMat <- sapply(trans, is.na)
trans <- subset(trans, rowSums(IndexMat) < 114) # 572*0.2=114.4

#' Identify samples without demographic info
## ----------------------------------------------------------------------------
rows_with_na <- which(apply(demo, 1, function(x) any(is.na(x))))
head(rows_with_na) # 39 such samples

demo <- demo[-rows_with_na,]
trans <- trans[,-rows_with_na]

#' 
#' Preprocesssing of demographic info that are gonna be used in the model
## ----------------------------------------------------------------------------
demo$consensus_ancestry <- as.factor(demo$consensus_ancestry)

demo$SUBTYPE <- as.factor(demo$SUBTYPE)

demo$batch <- as.factor(demo$batch)

#' 
#' ## Design and fit linear regression model
## ----------------------------------------------------------------------------
# Model design
design_full <- model.matrix(~ 0+ consensus_ancestry+ SUBTYPE+batch, data=demo)

# Model fitting
fit <- lmFit(trans, design_full)

#' 
#' Create a contrast matrix for specific comparisons
## ----------------------------------------------------------------------------
contMatrix <- makeContrasts("consensus_ancestryeur-consensus_ancestryafr", levels = design_full)

#' Fit model with contrasts
## ----------------------------------------------------------------------------
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

#' Results
## ----------------------------------------------------------------------------
DMPs <- topTable(fit2, num=Inf, coef=1,  adjust.method = "BH") 
head(DMPs)

#' Identify genes with significant pval
## ----------------------------------------------------------------------------
sig.probes <- DMPs[DMPs$adj.P.Val < 0.05,]
dim(sig.probes)

#' Identify genes with es >= 0.1
## ----------------------------------------------------------------------------
effective <- sig.probes[abs(sig.probes$logFC) >= 0.1,]
dim(effective)

#' Save all results
## ----------------------------------------------------------------------------
write.csv(DMPs,"diff-expr-results.csv",row.names = T)

