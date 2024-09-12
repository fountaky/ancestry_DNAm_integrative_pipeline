#' # Generate plots for Figure 4 of the paper
#' ## Load libraries and data
#' Libraries
## ------------------------------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggridges)
library("kSamples")

#' Load metadata (mean AFR and EUR probe methylation and gene expression value for each eQTM probe-gene pair)
## ------------------------------------------------------------------
data <- read_csv("eqtms-means-ancestry.csv")

data$dif_expr <- data$mean_expr_afr - data$mean_expr_eur  
data$dif_meth <- data$mean_meth_afr - data$mean_meth_eur 

#' 
#' Get number of eQTM hits per gene
## ------------------------------------------------------------------
hits <- as.data.frame(table(data$gene))
hits <- hits[order(hits$Freq, decreasing = T),]

## ------------------------------------------------------------------
colnames(hits) <- c("gene","hits")

#' Add expression and methylation differences to eqtms hits dataframe
## ------------------------------------------------------------------
merged <- merge(hits, data[,c("gene","dif_expr","dif_meth")],by="gene")

#' Expression as separate dataframe
## ------------------------------------------------------------------
expression <- unique(merged[,-4])

#' Methylation as separate dataframe
## ------------------------------------------------------------------
methylation <- unique(merged[,-3])

#' Create bining to group eQTM hits into different categories
## ------------------------------------------------------------------
# 5 hits per group
for (i in 1:dim(expression)[1]) {
  if (expression$hits[i] > 5 & expression$hits[i] < 11) {
    expression$group[i] <- "6-10"
  } else if (expression$hits[i] >= 11 & expression$hits[i] < 16){
    expression$group[i] <- "11-15"
  } else if (expression$hits[i] >= 16 & expression$hits[i] < 25){
    expression$group[i] <- "16-24"
  } else {
    expression$group[i] <- "1-5"
  }
}

#' Perform Anderson-Darling statistical test to examine whether difference in expression is associated with number of eQTM hits
## ------------------------------------------------------------------
ad.diff <- ad.test(expression$dif_expr[expression$group == "1-5"],expression$dif_expr[expression$group == "6-10"],expression$dif_expr[expression$group == "11-15"],expression$dif_expr[expression$group == "16-24"])
ad.diff

#' # Generate Plot
#' Ridgeplot
## ------------------------------------------------------------------
library(RColorBrewer)

a <- c(brewer.pal(9, "Set1"),brewer.pal(3,"Set2"),brewer.pal(4,"Accent"),"cyan2","#B3CDE3","darkred") # Set up distribution colots

# Create final dataframe and plot the distributions
prov <- expression
prov$group <- factor(prov$group, levels=unique(prov$group))

p1 <- ggplot(prov, 
             aes(x=log2(abs(dif_expr)),y=group,fill=group)) +  geom_density_ridges(alpha = 0.65)+  ylab("")+xlab("Expression Difference Per Gene [AFR-EUR]")+
  scale_fill_manual(values=a)  +
  theme(text = element_text(size=23),axis.text=element_text(size=21),plot.title =   element_text(size = 12,  face = "bold"),
        axis.title = element_text(),
        axis.text.x=element_text(color="black") ,
        axis.text.y=element_text(color="black") ,
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position="bottom") 

p2 <- ggpar(p1, legend.title = "No of eQTM Hits")

#' Save plot
## ------------------------------------------------------------------
p2
ggsave("eQTM_hits_ridge_plot.jpeg",width = 19, height = 22, device='jpeg', dpi=700,units = c("cm"))

