#' # Load libraries and data
#' Load libraries
## -------------------------------------------------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)

#' 
#' Load methylation data for the genes with the highest number of DMS (along with gene annotation information)
## -------------------------------------------------------------------------------------
dat1 <- read_csv("meth-levels-selected-genes.csv")
dat2 <- dat1

#' 
#' # Preprocess gene structure annotations
#' For multiple annotations we follow the priority order: TSS200 > TSS1500 > 5′UTR > 1st Exon > Body > 3′ UTR > Intergenic
## -------------------------------------------------------------------------------------
for (i in 1:dim(dat2)[1]) {
    if (grepl("TSS200", dat2$UCSC_RefGene_Group[i])){
        dat2$UCSC_RefGene_Group[i] <- "TSS200"
    } else if (grepl("TSS1500",dat2$UCSC_RefGene_Group[i])) { 
        dat2$UCSC_RefGene_Group[i] <- "TSS1500"
    } else if (grepl("5'UTR",dat2$UCSC_RefGene_Group[i])) {
        dat2$UCSC_RefGene_Group[i] <- "5'UTR"
    } else if (grepl("1stExon",dat2$UCSC_RefGene_Group[i])) {
        dat2$UCSC_RefGene_Group[i] <- "1stExon"
    } else if (grepl("Body",dat2$UCSC_RefGene_Group[i])) {
        dat2$UCSC_RefGene_Group[i] <- "Body"
    } else if (grepl("3'UTR",dat2$UCSC_RefGene_Group[i])) {
        dat2$UCSC_RefGene_Group[i] <- "3'UTR"
    } else {
        dat2$UCSC_RefGene_Group[i] <- "Intergenic"
    }
}

#' 
#' # Load mean methylation values (along with std) per population group for differentially methylated probes
## -------------------------------------------------------------------------------------
std <- read_csv("probes-means-meth-ancestry.csv")
dat3 <- merge(dat2,std[,c("probe","mean_afr","mean_eur","std_afr","std_eur")],by=c("probe","mean_afr","mean_eur"))

#' # Generate stacked barplot -Mean methylation per probe per gene
#' Appropriate formatting of data to be visualized
## -------------------------------------------------------------------------------------
to_plot<- data.frame(probe=rep(dat3$probe,2),gene=rep(dat3$UCSC_RefGene_Name,2),meth_levels=c(dat3$mean_afr,dat3$mean_eur),ancestry=c(rep("AFR",dim(dat3)[1]),
                     rep("EUR",dim(dat3)[1])),std=c(dat3$std_afr,dat3$std_eur),annotation1=rep(dat3$UCSC_RefGene_Group,2),annotation2=rep(dat3$Relation_to_UCSC_CpG_Island,2))

to_plot <- to_plot[order(to_plot$gene,to_plot$annotation1,to_plot$annotation2),]
to_plot$probe <- factor(to_plot$probe,levels=unique(to_plot$probe))

#' Use alternative gene name for C21orf56
## -------------------------------------------------------------------------------------
a <- which(to_plot$gene == "C21orf56")
to_plot$gene[a] <- "SPATC1L"

#' Generate methylation barplots per gene
## -------------------------------------------------------------------------------------
g = "TNXB" # gene to be plotted

p1 <- ggplot(to_plot[which(to_plot$gene == g),], aes(factor(ancestry,levels=c("AFR","EUR")), meth_levels,fill=ancestry)) + 
  geom_bar(stat='identity') +
  theme(text = element_text(size=20),strip.text.x = element_blank(),panel.background = element_blank(),axis.line = element_line(color = "black"),
        axis.text.x=element_blank() ,
        axis.text.y=element_text(color="black"),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.7),legend.position="none") + scale_fill_manual(values = c("orange","royalblue"))+
  facet_grid(~ probe)+xlab("") + ylab("Methylation Levels")+
  geom_errorbar(aes(ymin = meth_levels - std, ymax = meth_levels + std), color = "Black", size = .2, width = .3) +ylab("")


p1
ggsave(paste0("C:/Users/KFounta/OneDrive - Northwell Health/Founta-Chambwe-Shared-Folder/Manuscripts&Abstracts/manuscripts/Founta-et-al/figures/paper figures/Figure-2/R-figures/diff-meth-",g,".jpeg"), width = 8, height = 4, device='jpeg', dpi=700)
