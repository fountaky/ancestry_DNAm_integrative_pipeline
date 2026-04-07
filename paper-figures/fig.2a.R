# Load libraries
library(readr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Load statistics for all illumina array probes
paired_test <- read_csv("BRCATopTable_afr_eur.csv")
paired_test <- as.data.frame(paired_test)
row.names(paired_test) <- paired_test[,1]
paired_test <- paired_test[,-1]

# Load statistics for probes with pval < 0.05
sig.probes <- read_csv("BeforeBRCAAncestryProbes.csv")
sig.probes <- as.data.frame(sig.probes)
row.names(sig.probes) <- sig.probes[,1]
sig.probes <- sig.probes[,-1]

# Load statistics for probes with pval < 0.05 & effect size >= 0.1
effective <- read_csv("AfterBRCAAncestryProbes.csv")
effective <- as.data.frame(effective)
row.names(effective) <- effective[,1]
effective <- effective[,-1]

hyper_eur <- effective[effective$logFC > 0,]
hyper_afr <- effective[effective$logFC < 0,]

# Remove "rs" probes
paired_test <- paired_test[-which(grepl("rs",row.names(paired_test))),]

# Create Volcano plot
## Add column to be used for color annotation
color <- c()

for (i in 1:nrow(paired_test)){
  
  if(row.names(paired_test)[i] %in% row.names(sig.probes) & row.names(paired_test)[i] %in% row.names(hyper_eur)){
    color[i] <- "Significant_Effective_eur"
  } else if (row.names(paired_test)[i] %in% row.names(sig.probes) & row.names(paired_test)[i] %in% row.names(hyper_afr)) {
    color[i] <- "Significant_Effective_afr"
  } else if(row.names(paired_test)[i] %in% row.names(sig.probes) & !(row.names(paired_test)[i] %in% row.names(effective))){
    color[i] <- "Significant_Not_Effective" 
  } else{
    color[i] <- "Non_Significant"
  }
}

paired_test$color <- color
length(which(paired_test$color == "Significant_Effective_eur"))
length(which(paired_test$color == "Significant_Effective_afr"))
length(which(paired_test$color == "Significant_Not_Effective"))

## Generate volcano plot
paired_test$color <- as.factor(paired_test$color)

v1 <- ggplot(data=paired_test, aes(x=logFC, y=-log10(adj.P.Val), col=color, label=color)) + 
  geom_point() + 
  theme_minimal() +
  geom_text(fontface = "bold",position=position_jitter(width=0.1,height=0.05))

v2 <- ggplot(data=paired_test, aes(x=logFC, y=-log10(adj.P.Val), col=color)) + #, label=delabel
  geom_point() + 
  theme(text = element_text(size=23),
        axis.title = element_text(),
        axis.text.x=element_text(color="black") ,
        axis.text.y=element_text(color="black") ,
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black") ,
        legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, size=0.7))+
  scale_color_manual(values=c("black","orange","royalblue", "snow3")) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype='dashed', col = 'black') +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', col = 'black') + xlab("Fold Change (log)") + ylab("Significance (-log)")

v2 <- ggpar(v2, legend.title = "")
v2
ggsave("Volcano_diff_probes.jpeg",width = 16, height = 16, device='jpeg', dpi=700,units = c("cm"))
