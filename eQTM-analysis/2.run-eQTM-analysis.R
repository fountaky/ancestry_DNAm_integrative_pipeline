library(readr)
library(stats)

# Load probe-gene pairs (candidate eQTMs)
probe_gene_pairs <- read_csv("probe-gene-pairs.csv")
View(probe_gene_pairs)

# Load TCGA gene expression data
expression <- read_csv("expression_brca_tcga.csv")

# Load BRCA methylation beta matrix
methylation <- read_csv("BRCA_betas_selected_probes.csv")
colnames(methylation) <- c("probe",substr(colnames(methylation[,-1]), start=1, stop = 12)) # Process TCGA barcode

methylation <- methylation[,c("probe",intersect(colnames(methylation), colnames(expression)))]

methylation <- as.data.frame(methylation)
row.names(methylation) <- methylation[,1]
methylation <- methylation[,-1]

expression <- as.data.frame(expression)
expression <- aggregate(expression[,-1], list(expression$gene_id),data=expression,FUN=mean)

row.names(expression) <- expression[,1]
expression <- expression[,-1]

pairs <- probe_gene_pairs[,c(6,12)]
colnames(pairs) <- c("probe", "gene_symbol")

overall_results <- data.frame(probe = pairs$probe, gene=pairs$gene_symbol)

# Run association analysis for every probe-gene pair
for (i in 1:dim(pairs)[1]) { 
  temp_meth <- t(methylation[pairs$probe[i],])
  temp_exp <- t(expression[pairs$gene_symbol[i],])
  
  correlation <- cor.test(temp_meth[,1], temp_exp[,1], method = c("spearman"),adjust.method="fdr",exact=FALSE)
  
  overall_results$p.val[i] <- correlation$p.value
  overall_results$rho[i] <- correlation$estimate
  overall_results$statistic[i] <- correlation$statistic
  
}

# Calculate adjusted p-values
overall_results$p.adjusted <- p.adjust(overall_results$p.val, method = "BH")

eqtms <-overall_results[which(overall_results$p.adjusted <= 0.05),]
dim(eqtms)

# Save significant results
write.csv(eqtms, "sig_eQTMs.csv", row.names = F)

# Save all results
write.csv(overall_results, "overall_eQTMs_revised.csv", row.names = F)

