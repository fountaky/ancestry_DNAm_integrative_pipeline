## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(readr)
library("org.Hs.eg.db", character.only = TRUE)
require(DOSE)
library(GO.db)
library(clusterProfiler)
library(enrichplot)

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
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


#' Circular plot
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
edox <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID') ## convert gene ID to Symbol
p3 <- cnetplot(edox, foldChange=gene_list, circular = TRUE,node_label="all", colorEdge = TRUE,showCategory=10,color_category = "gray48",color_gene="red",cex_label_gene=1.5) 

pdf("genes-gsea-circular-10-groups.pdf",25,15)
p3
dev.off()

#' Sugiyama plot
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf("sugiyama_gsea.pdf", 15,9)
goplot(gse,layout = "sugiyama")
dev.off()

