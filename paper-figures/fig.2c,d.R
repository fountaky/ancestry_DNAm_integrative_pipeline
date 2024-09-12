#Load libraries
library(readr)
library(ggplot2)
library(ggpubr)

# Load odds ratios for each ancestry group of interest
ancestry="eur"
odds_ratios <- read_csv(paste0("C:/Users/KFounta/OneDrive - Northwell Health/Founta-Chambwe-Shared-Folder/Projects/Thesis/Results/odds_ratios/odds_rations_",ancestry,".csv"))
odds_ratios <- as.data.frame(odds_ratios)

# Transform significance to asterisk symbol
# Asterisk anntoations for odds ratios p-val
anns <- c()
for (i in 1:dim(odds_ratios)[1]) {
  if (odds_ratios$`p-value`[i] >= 0.05) {
    anns <- c(anns,"")
  } else if (odds_ratios$`p-value`[i] < 0.05 & odds_ratios$`p-value`[i] >= 0.03) {
    anns <- c(anns,"*")
  } else if (odds_ratios$`p-value`[i] < 0.03 & odds_ratios$`p-value`[i]>= 0.01) {
    anns <- c(anns,"**")
  }  else {
    anns <- c(anns,"***")
  }
}  


# Fix column names
colnames(odds_ratios)[c(1)] <- c("Annotation")

# Plot CGI regions
cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
       "#FFDB6D", "#C4961A", "grey", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

# Generate and save plots
odds_ratios$l =odds_ratios$Annotation

cgi <- c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea")
to_plot <- odds_ratios[which(odds_ratios$Annotation %in% cgi),]

plot1 <- ggplot(to_plot, aes(y= odds-1, x = l,fill=l)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_pointrange(aes(ymin=odds_conf_low-1, ymax=odds_conf_high-1), width=0.5,size=0.5) +
  geom_hline(yintercept = 0, linetype=2) +
  coord_flip() +
  labs(x="", y = "Odds Ratio") +
  scale_fill_manual(values = cbp1[1:6])+
  theme(plot.title = element_text(size = 23,  face = "bold"),
        text = element_text(size = 23),
        axis.title = element_text(),
        axis.text.x=element_text(color="black") ,
        axis.text.y=element_text(color="black") ,
        legend.position = "bottom",
        plot.margin = margin(2, 1, 1, 1, "cm"),panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.7)) +
  scale_y_continuous(labels = function(y) y + 1) + 
  geom_text(aes(x = c(1:dim(to_plot)[1]),  y = rep(2.8,dim(to_plot)[1]), label = anns[c(8:13)]), size = 7)


plot1
ggsave(paste0("figure_",ancestry,"_odds_ratios_CGI.jpeg"), width = 8, height = 6, device='jpeg', dpi=700)


gene_group <- c("TSS200","1stExon","5'UTR","3'UTR","TSS1500","Body","Intergenic")
to_plot <- odds_ratios[which(odds_ratios$Annotation %in% gene_group),]

plot2 <-ggplot(to_plot, aes(y= odds-1, x = l,fill=l)) +
  geom_bar(stat = "identity") +
  geom_pointrange(aes(ymin=odds_conf_low-1, ymax=odds_conf_high-1), width=0.5,size=0.5) +
  geom_hline(yintercept = 0, linetype=2) +
  coord_flip() +
  labs(x="", y = "Odds Ratio") +
  scale_fill_manual(values = cbp1[7:15])+
  theme(plot.title = element_text(size = 23,  face = "bold"),
        text = element_text(size = 23),
        axis.title = element_text(),
        axis.text.x=element_text(size = 15,color="black") ,
        axis.text.y=element_text(color="black") ,
        legend.position = "bottom",
        plot.margin = margin(2, 1, 1, 1, "cm"),panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=0.7)) +
  scale_y_continuous(labels = function(y) y + 1) + 
  geom_text(aes(x = c(1:dim(to_plot)[1]),  y = rep(4.6,dim(to_plot)[1]), label = anns[c(1:7)]), size = 7)

plot2
ggsave(paste0("figure_",ancestry,"_odds_ratios_RefGeneGroup.jpeg"), width = 8, height = 6, device='jpeg', dpi=700)
