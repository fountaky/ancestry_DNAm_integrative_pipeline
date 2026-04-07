### Generate barplots to visualize pop specific frequencies of meSNPs ###
library(readr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Load dataframe with number of SNPs for different ancestry bins per population for ancestry-related and control permutation meSNPs
dat <- read_csv("POP_MAF_meQTLs_Revised.csv")

dat$freq <- 100*as.numeric(dat$freq)
dat$category <- factor(dat$category,levels=c("≤1%","<5%", "≥5%"))

dat$type <- c(rep("Ancestry \nmeSNPs",6),rep("Control \nmeSNPs",6))

# Now generate plot
p <- ggplot(data=dat, aes(x=category, y=freq, fill=type)) +
  geom_bar(stat="identity", width=0.8,position=position_dodge(width=0.9))+
 scale_fill_manual(values=c("darksalmon","lavenderblush4"))+ylab("%meSNPs")+
  scale_x_discrete(labels=c("≤1%", "<5%", "≥5%"))+
  xlab("gnomAD Allele Frequency")+
  theme(text=element_text(size=20),plot.title = element_text(size = 17,  face = "bold"),
         axis.title.y = element_text(size = 25),
         axis.title.x = element_text(size = 25),
         axis.text.x=element_text(color="black"),
         axis.text.y=element_text(color="black") ,
         panel.background = element_blank(),axis.line = element_line(colour = "black"),
         axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "right",
        legend.title=element_blank(),
         legend.key = element_blank())+
  facet_grid(~pop, labeller = labeller(type = c(
   AFR = "AFR Frequency",
   EUR = "EUR Frequency")))

p

# Save plot
ggsave("POP_FREQ_plots_meQTLs.jpeg", width = 19, height = 14, device='jpeg', dpi=700,units = c("cm"))
