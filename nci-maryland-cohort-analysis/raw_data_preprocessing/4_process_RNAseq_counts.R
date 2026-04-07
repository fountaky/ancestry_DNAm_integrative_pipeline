library(readr)
library(data.table)
library(DESeq2)
library(ggplot2)
library(ggpubr)

# Import raw counts for tumor RNAseq samples 
data <- read_csv("/raw/expression/tumor_raw_reads.csv")
data <- as.data.frame(data)
row.names(data) <- data[,1]
data <- data[,-1]

# Create basic metadata
sample_names <- colnames(data)

### Remove low count genes
# Remove genes expressed in less than 80% percent of the samples
threshold <- 0.8 * ncol(data)  # 80% of the samples

# Identify genes (rows) where the number of nonzero counts exceeds the threshold
keep_genes <- rowSums(data != 0) > threshold

# Subset matrix to keep only those genes
filtered_mat <- data[keep_genes, ]
dim(filtered_mat)

# Load expression metadata obtained from GEO
expr_metadata <- read_csv("expr_metadata.csv")
expr_metadata$title<- gsub("S_","", expr_metadata$title)

expr_metadata <- as.data.frame(expr_metadata)
row.names(expr_metadata) <- expr_metadata$title

expr <- expr_metadata[colnames(filtered_mat),]
all.equal(colnames(filtered_mat), row.names(expr))

colnames(expr)[52] <- "race"

# Optional: add race information in the metadata dataframe
metadata <- data.frame(sample = sample_names, race=expr$race)
rownames(metadata) <- sample_names

# Round decimals --> they should be integers since the values represent raw counts
filtered_mat <- round(filtered_mat)

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_mat,
                              colData = metadata,
                              design = ~ race)  # No actual comparison is being performed here

# Normalize data
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Generate log2 transformed data
vst <- vst(dds, blind=FALSE)

# Extract the data
vst_matrix <- assay(vst)

# Diagnostic PCA plot 
plotPCA(vst,intgroup = "race")+ geom_text(aes(label = name))+theme(text = element_text(size = 15),axis.text.x = element_text(angle=0, hjust=1),
                                                                         legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(color = "black"),
                                                                         axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.7))                          

# Save log transformed data
write.csv(vst_matrix,"deseq2_processed_tumors.csv",
          row.names = T)
