# Required libraries
library(ggplot2)
library(readr)
library(data.table)
library(gridExtra)

# Load methylation beta value matrix
matrix <- read.csv("expression.csv")
matrix <- matrix[,-c(1:3)]
colnames(matrix)[1] <- "id"    

# Load table with genotypes across all chrosomes for every samples in our cohort
genotypes <- fread("genotypes_all_chr.txt")
genotypes <- genotypes[,-1]

# Load demographic information for our cohort
demo <- read.csv("demographics_legacy_ids.csv")
colnames(demo) <- demo[1,]
demo <- demo[-1,]
demo <- t(demo)
demo <- as.data.frame(demo)
colnames(demo) <- demo[1,]
demo <- demo[-1,]

# Load trans-meQTLs (along with integrated AIM SNP annotations)
tensor_trans <- read.csv("trans_3_annotations_es_0.25.csv")
tensor_trans <- tensor_trans[,-c(1:2)]

# Load cis-meQTLs (along with integrated AIM SNP annotations)
tensor_cis <- read.csv("cis_3_annotations_es_0.25.csv")
tensor_cis <- tensor_cis[,-c(1:2)]

all_meqtl_snps <- rbind(tensor_cis[,c("snps","rs")],tensor_trans[,c("snps","rs")])

# Create Function to plot meQTLs
## meQTL visualization function
plot_meqtl <- function(snp_id, probe_id){

    GT <- genotypes[genotypes$id == snp_id,]
    b_value <- matrix[matrix$id == probe_id,]

    to_plot <- as.data.frame(rbind(GT, b_value))
    row.names(to_plot) <- to_plot[,1]
    to_plot <- to_plot[,-1]
    to_plot[1,] <- as.factor(to_plot[1,])

    df <- t(to_plot)
    df <- as.data.frame(df)
    colnames(df) <- c("GT","beta_value")
    df <- na.omit(df)

    if (length(which(df$GT == "NA")) != 0){ 
      df <- df[-which(df$GT == "NA"),]
    }

    df$beta_value <- as.numeric(df$beta_value)
    df$GT <- as.factor(df$GT)
    
    a <- which(all_meqtl_snps$snps == snp_id)
    if (!is.na(all_meqtl_snps$rs[a]))  {

        p<-ggplot(df, aes(x=GT, y=beta_value,fill=GT)) + 
          geom_boxplot()+
          geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + ggtitle(paste0(all_meqtl_snps$rs[a],", ",probe_id)) +
          labs(y= "DNA Methylation Level", x="Genotype")+theme(axis.text=element_text(size=20),plot.title = element_text(size = 16,  face = "bold"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black") ,
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "none")+
        scale_fill_manual(values=c("#7ba2fc","#9cff75","#f20060"))
        
        p
        
    } else {
        p<-ggplot(df, aes(x=GT, y=beta_value,fill=GT)) + 
          geom_boxplot()+
          geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + ggtitle(paste0(snp_id,", ",probe_id)) +
          labs(y= "DNA Methylation Level", x="Genotype")+theme(axis.text=element_text(size=20),plot.title = element_text(size = 16,  face = "bold"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black") ,
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(color = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position = "none")+
        scale_fill_manual(values=c("#7ba2fc","#9cff75","#f20060"))
        
        p
    }
    }