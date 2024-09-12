#' 
## --------------------------------------------------------------------------------
library(readr)
library("org.Hs.eg.db", character.only = TRUE)
library(GO.db)
library(clusterProfiler)

#' 
#' Load gene list with number of eQTM hits per gene
## --------------------------------------------------------------------------------
gene_list_gsea <- read_csv("eQTM_gene_list_gsea.csv")

#' Transform dataframe to gene list
## --------------------------------------------------------------------------------
gene_list <-gene_list_gsea$x
names(gene_list) <- gene_list_gsea$...1

#' Download most recent BP GO
## --------------------------------------------------------------------------------
recent_go_BP <- gson_GO(org.Hs.eg.db, keytype = "ENTREZID", ont = "BP")

#' Run GSEA
## --------------------------------------------------------------------------------
set.seed(1) 
gse <- GSEA(
       gene_list,
        exponent = 1,
        minGSSize = 3,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        gson = recent_go_BP, 
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")

View(gse@result)

#' Save results
## --------------------------------------------------------------------------------
write.csv(gse@result,"eQTM-genes-GSEA-results.csv", row.names = T)

