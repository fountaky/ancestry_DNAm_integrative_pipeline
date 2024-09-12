#' ## Visualization of GSEA eQTM results using ggplot2 barplots
#' 
#' Load libraries
## ---------------------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)

#' Load GSEA results
## ---------------------------------------------------------------------------
data <- read_csv("eQTM-genes-GSEA-results.csv")

#' Data preprocessing
## ---------------------------------------------------------------------------
to_plot <- data[,c("ID", "Description","setSize","p.adjust")] # select field of interest
to_plot <- to_plot[order(to_plot$p.adjust,decreasing = F,to_plot$setSize),] 
to_plot <- to_plot[c(1:10),] # select first 10 most significantly eQTM enriched ontologies

to_plot <- to_plot[order(to_plot$p.adjust,to_plot$setSize,decreasing = T,to_plot$setSize),] # change ordering for plotting purposes
to_plot$Description <- factor(to_plot$Description,levels=to_plot$Description)

#' Generate plot
## ---------------------------------------------------------------------------
plt <- ggplot(to_plot) +
  geom_col(aes(-log10(p.adjust), Description, fill = setSize), width = 0.6) +theme(text = element_text(size=15),axis.text=element_text(size=15),plot.title = element_text(size = 15,  face = "bold"),
        
         element_text(size=20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black") ,
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1))+ylab("")+scale_fill_viridis(option = "rocket",direction = -1)+ xlab("Significance (-log10)")

p2 <- ggpar(plt, legend.title = "Set Size",legend = "bottom")

#' Save plot
## ---------------------------------------------------------------------------
p2
ggsave("all-eQTM-genes-GSEA-barplots.jpeg", width = 32, height = 22, device='jpeg', dpi=700,units = c("cm"))

