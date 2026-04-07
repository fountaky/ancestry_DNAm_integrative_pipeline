#' Load Libraries
## -------------------------------------------------------------------------------------
library(ggplot2)
library(readr)
library(ggpubr)
library(dplyr)

#' Load cis- and trans-meQTLs (along with AIM SNP annotations)
## -------------------------------------------------------------------------------------
cis <- read_csv("cis_3_annotations_es_0.25.csv")
trans <- read_csv("trans_3_annotations_es_0.25.csv")

#' Keep only columns of interest and combine data
## -------------------------------------------------------------------------------------
cis <- cis[,c("probe","snps","rs","Significant")]
cis$type <- "cis"
trans <- trans[,c("probe","snps","rs","Significant")]
trans$type <- "trans"

merged <- rbind(cis,trans)

#' ## Load SNP locations
## -------------------------------------------------------------------------------------
locations <- read_csv("all_snps_locations_meQTLs_0.25.csv")
colnames(locations)[1] <- "snps"

#' Add SNP location information to meQTL dataframe
## -------------------------------------------------------------------------------------
final <- merge(merged,locations[,c("snps","Chr","Pos")],by="snps")

#' Order based on SNP location
## -------------------------------------------------------------------------------------
library(dplyr)
final <- final %>% mutate_at(c('Chr', 'Pos'), as.numeric)

final <- final[(order(final$Chr,final$Pos)),]

## -------------------------------------------------------------------------------------
type <- final[,c("probe","type")]
type <- unique(type)

#' Identify probes regulated by both cis- and trans-meQTLs
## -------------------------------------------------------------------------------------
a <- type$probe[duplicated(type$probe)]
type[which(type$probe %in% a),"type"] <- "cis&trans"

# Now remove duplicated rows
type <- unique(type)

#' Load probe annotations
## -------------------------------------------------------------------------------------
annotations <- read_csv("diff-probes-annotations.csv")
colnames(annotations)[1] <- "probe"

toplot <- annotations

#' 
#' For multiple UCSC_RefGene_Group annotations we follow the priority order: TSS200 > TSS1500 > 5′UTR > 1st Exon > Body > 3′ UTR > Intergenic
## -------------------------------------------------------------------------------------
for (i in 1:dim(toplot)[1]) {
    if (grepl("TSS200", toplot$UCSC_RefGene_Group[i])){
        toplot$UCSC_RefGene_Group[i] <- "TSS200"
    } else if (grepl("TSS1500",toplot$UCSC_RefGene_Group[i])) { 
        toplot$UCSC_RefGene_Group[i] <- "TSS1500"
    } else if (grepl("5'UTR",toplot$UCSC_RefGene_Group[i])) {
        toplot$UCSC_RefGene_Group[i] <- "5'UTR"
    } else if (grepl("1stExon",toplot$UCSC_RefGene_Group[i])) {
        toplot$UCSC_RefGene_Group[i] <- "1stExon"
    } else if (grepl("Body",toplot$UCSC_RefGene_Group[i])) {
        toplot$UCSC_RefGene_Group[i] <- "Body"
    } else if (grepl("3'UTR",toplot$UCSC_RefGene_Group[i])) {
        toplot$UCSC_RefGene_Group[i] <- "3'UTR"
    } else {
        toplot$UCSC_RefGene_Group[i] <- "Intergenic"
    }
}

## -------------------------------------------------------------------------------------
type <- merge(type, toplot[,c("probe","Relation_to_UCSC_CpG_Island","UCSC_RefGene_Group")],by="probe") # Transfer all annotations to main dataframe
type[which(is.na(type$Relation_to_UCSC_CpG_Island)),"Relation_to_UCSC_CpG_Island"] <- "OpenSea"

#' Load chromatin states annotations
## -------------------------------------------------------------------------------------
states <- read_csv("probes-universal-chr-states.csv")
colnames(states) <- c("probe","chr_state")

#' Incorporate chromatin states annotations into main dataframe
## -------------------------------------------------------------------------------------
type <- merge(type, states,by="probe")

#' ## Generate appropriate format to plot meQTL events
## -------------------------------------------------------------------------------------
matrix <- table(final[,c("probe","snps")])
matrix <- (matrix > 0) + 0

#' Further format processing
## -------------------------------------------------------------------------------------
plot <- as.data.frame(t(matrix))
plot <- plot[,type$probe]
head(plot)

#' # Order based on decreasing order of meQTL hits per probe
## -------------------------------------------------------------------------------------
plot <- plot[,order(colSums(plot),decreasing = T)]

#' ## Add SNP significance info
## -------------------------------------------------------------------------------------
sig <- final[c("rs","Significant","snps","Chr","Pos")]

sig <- as.data.frame(sig)
sig <- unique(sig)

sig$Significant[which(is.na(sig$Significant))] <- "Significant_FALSE"

sig <- as.data.frame(sig)
row.names(sig) <- sig$snps

## -------------------------------------------------------------------------------------
final_sig <- merge(plot,sig,by=0)

row.names(final_sig) <- final_sig$Row.names
final_sig <- final_sig[,-1]
 
#' Order based on significance "TRUE" or "FALSE"
## -------------------------------------------------------------------------------------
final_sig <- final_sig[order(final_sig$Significant,final_sig$Chr,final_sig$Pos,decreasing = T),]
head(final_sig)

#' Replace number encoders and non suitable names with appropriate character variables
## -------------------------------------------------------------------------------------
rownames <- row.names(final_sig)

final_sig <- final_sig %>%
  mutate(across(everything(), as.character))

final_sig[,c(1:198)]<- data.frame(lapply(final_sig[,c(1:198)], function(x) {
 gsub("0", "NO", x)
}))

final_sig[,c(1:198)]<- data.frame(lapply(final_sig[,c(1:198)], function(x) {
 gsub("1", "YES", x)
}))

final_sig<- data.frame(lapply(final_sig, function(x) {
 gsub("Significant_TRUE", "TRUE", x)
}))

final_sig<- data.frame(lapply(final_sig, function(x) {
 gsub("Significant_FALSE", "FALSE", x)
}))

row.names(final_sig) <- rownames

#' ## Add probe meQTL type info
#' First, you need to order "type" dataframe based on probe order in the main "final_sig" dataframe
## -------------------------------------------------------------------------------------
type <- as.data.frame(type)
row.names(type) <- type$probe

type$chr_state <- gsub(".*_","",type$chr_state)

type <- type[colnames(final_sig)[1:198],]

#' Order "final_sig" based on significance, chr and position
## -------------------------------------------------------------------------------------
final_sig <- final_sig %>% mutate_at(c('Chr', 'Pos'), as.numeric)
final_sig <- final_sig[order(final_sig$Significant,final_sig$Chr,final_sig$Pos),]

#' ## Create the heatmap using ComplexHeatmap
#' Load plot specific libraries
## -------------------------------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(randomcoloR)

#' Generate annotations
## -------------------------------------------------------------------------------------
### Relation to CGI
type$Relation_to_UCSC_CpG_Island <- factor(type$Relation_to_UCSC_CpG_Island,levels=c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea"))
cgi_col = c("#E69F00", "#0072B2","#56B4E9", "#009E73",  "#D55E00", "#F0E442")
names(cgi_col) <- levels(type$Relation_to_UCSC_CpG_Island)

### Relation to Gene regulatory element
type$UCSC_RefGene_Group <- factor(type$UCSC_RefGene_Group,levels=c("1stExon","TSS200","5'UTR","3'UTR","TSS1500","Body","Intergenic"))
genegroup_col = c("#CC79A7", "#52854C","#C4961A", "#FFDB6D", "#C3D7A4", "grey", "#D16103")
names(genegroup_col) <- levels(type$UCSC_RefGene_Group)

# Heatmap Column split
col_split = type$type
split = final_sig$Significant

# Row annotations
left_annot <- rowAnnotation('AIM SNP' = split,col = list('AIM SNP' = c("TRUE" = "lightblue2", "FALSE" = "cornsilk4")))

# Column annotation on top
ha_column=HeatmapAnnotation(border = TRUE,'meQTL regulation' = type$type,col=list('meQTL regulation'=c("cis"="lavender","cis&trans"="burlywood1","trans"="indianred1")))

# Column annotation bottom
ha_column_bottom = HeatmapAnnotation(border = TRUE,"CGI Annotation" = type$Relation_to_UCSC_CpG_Island,"GeneStructure Annotation"=type$UCSC_RefGene_Group,col=list("CGI Annotation"=cgi_col,"GeneStructure Annotation"=genegroup_col))

colors = c("gray97","gray7") # complex heatmap colors

#' 
#' Now generate heatmap
## -------------------------------------------------------------------------------------
set.seed(12)

heat <- Heatmap(as.matrix(final_sig[,c(1:198)]), name = "meQTL Regulation", cluster_rows = F, cluster_columns = T,cluster_column_slices=T, left_annotation = left_annot,show_row_names = F,
                show_parent_dend_line = FALSE,row_split = split,border = TRUE, row_title = NULL, show_column_names = F, cluster_row_slices = TRUE, top_annotation = ha_column,
                row_gap = unit(2, "mm"),col = colors,column_split = col_split,bottom_annotation =ha_column_bottom) 

#' Save heatmap
## -------------------------------------------------------------------------------------
pdf("meQTLs-heatmap-0.25-legend.pdf",24,12)
draw(heat, annotation_legend_side = "bottom",heatmap_legend_side="bottom",merge_legend=T)
dev.off()

