#Load libraries
library(readr)
library(readxl)
library(dplyr)
library(RColorBrewer)

# Load methylation matrix with beta values of for AFR and EUR patients across the differentially methylated probes
betas <- read_csv("brca_meth_afr_eur_only.csv")
data <- as.data.frame(betas)
row.names(data) <- data[,1]
data <- data[,-1]

data <- data[-which(grepl("rs",row.names(data))),] # remove "rs" probes

# Load and process demographic information for TCGA patients
demo <- read_csv("FinalFullTCGAmeta2.csv")
demo <- as.data.frame(demo)

colnames(demo)[3] <- "sample"
row.names(demo) <- demo$sample

# Demographics only for cohort of interest
brca_demo <- demo[colnames(data),]

# Sanity check
all.equal(row.names(brca_demo),colnames(data))

brca_demo$consensus_ancestry <- toupper(brca_demo$consensus_ancestry) # ancestry categories uppercase letters

# Create annotations and order betas based on ancestry
df1 <- as.data.frame(brca_demo[,c("SUBTYPE","consensus_ancestry")])
df1 <- df1[order(df$consensus_ancestry),]

df1$patient <- substr(row.names(df1),1,12) # truncate tcga barcode to corespond to patient id

# Ancestry annotations for heatmap
ancestry <- as.data.frame(df1$consensus_ancestry)
row.names(ancestry) <- row.names(df1)
colnames(ancestry) <- "Genetic Ancestry"

# order samples in betas matrix based on ancestry category
to_plot <- data[,row.names(ancestry)]

# Load probe ancestry specific hypermethylation (hypermethylated in AFR or EUR)
probes_type <- read.csv("diff_probes_type.csv")
type <- probes_type

# Load annotations for differentially metylated probes
annotations <- read_csv("diff-probes-annotations.csv")
colnames(annotations)[1] <- "probe"
toplot <- annotations

# For multiple UCSC_RefGene_Group annotations 
# We follow the priority order: TSS200 > TSS1500 > 5′UTR > 1st Exon > Body > 3′ UTR > Intergenic
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

# Merge annotations with hypermethylation type
type <- merge(type, toplot[,c("probe","Relation_to_UCSC_CpG_Island","UCSC_RefGene_Group")],by="probe")

# OpenSea annotation for NA CGIs
type[which(is.na(type$Relation_to_UCSC_CpG_Island)),"Relation_to_UCSC_CpG_Island"] <- "OpenSea"

# Load and process chromatin states annotations
states <- read_csv("probes-universal-chr-states.csv")

colnames(states) <- c("probe","chr_state")
type <- merge(type, states,by="probe")

type$chr_state <- gsub(".*_","",type$chr_state)

# Order "type" dataframe based on ancestry sepcific hypermethylation
type <- type[order(type$type),]

# Order methylation matrix row.names based on "type" dataframe
to_plot <- to_plot[type$probe,]

###### Create the heatmap using ComplexHeatmap
# Generate annotations
library(tidyverse)
library(ComplexHeatmap)
library(randomcoloR)

# Color for subtypes
n = length(unique(df$SUBTYPE))
set.seed(12)
col_subtypes <- brewer.pal(n=5,"Set2")
names(col_subtypes) <- unique(df1$SUBTYPE)

### Chrom states - Annotate states into higher level categories
type[grepl("PromF",type$chr_state),"chr_state"] <- "PromF"
type[grepl("BivProm",type$chr_state),"chr_state"] <- "BivProm"
type[grepl("ReprPC",type$chr_state),"chr_state"] <- "ReprPC"
type[grepl("TSS",type$chr_state),"chr_state"] <- "TSS"
type[grepl("EnhA",type$chr_state),"chr_state"] <- "EnhA"
type[grepl("DNase",type$chr_state),"chr_state"] <- "DNase"
type[grepl("TxWk",type$chr_state),"chr_state"] <- "TxWk"
type[grepl("TxEnh",type$chr_state),"chr_state"] <- "TxEnh"
type[grepl("Quies",type$chr_state),"chr_state"] <- "Quies"
type[grepl("HET",type$chr_state),"chr_state"] <- "HET"
type[grepl("Acet",type$chr_state),"chr_state"] <- "Acet"
type[grepl("TxEx",type$chr_state),"chr_state"] <- "TxEx"
type[grepl("EnhWk",type$chr_state),"chr_state"] <- "EnhWk"
type[grepl("GapArtf",type$chr_state),"chr_state"] <- "GapArtf"
type[grepl("znf",type$chr_state),"chr_state"] <- "znf"
type[grepl("Tx2",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx6",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx5",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx8",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx1",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx7",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx3",type$chr_state),"chr_state"] <- "Tx"
type[grepl("Tx4",type$chr_state),"chr_state"] <- "Tx"

states <- unique(type$chr_state)
states <- states[order(states)]

type$chr_state <- factor(type$chr_state,levels=states)

states_col = c("lemonchiffon","purple","#fff44f","orange","yellow","#fff5ee","slateblue1",
               "#ff4500","azure2","dimgrey","red","darkgreen","#ADFF2F","#3cb371","forestgreen","#7fffd4")

names(states_col) <- levels(type$chr_state)

### Relation to CGI
type$Relation_to_UCSC_CpG_Island <- factor(type$Relation_to_UCSC_CpG_Island,levels=c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea"))
cgi_col = c("#E69F00", "#0072B2","#56B4E9", "#009E73",  "#D55E00", "#F0E442")
names(cgi_col) <- levels(type$Relation_to_UCSC_CpG_Island)

### Relation to Gene strucure
type$UCSC_RefGene_Group <- factor(type$UCSC_RefGene_Group,levels=c("1stExon","TSS200","5'UTR","3'UTR","TSS1500","Body","Intergenic"))
genegroup_col = c("#CC79A7", "#52854C","#C4961A", "#FFDB6D", "#C3D7A4", "grey", "#D16103")
names(genegroup_col) <- levels(type$UCSC_RefGene_Group)

# Column split (based on genetic ancestry)
col_split = ancestry$`Genetic Ancestry`
split = type$type

# Row annotations
left_annot <- rowAnnotation("Hypermethylation" = type$type,col=list("Hypermethylation"=c("AFR"="orange","EUR"="royalblue")))                                                                                                          

right_annot <- rowAnnotation("CpG Island" = type$Relation_to_UCSC_CpG_Island,"Gene Structure"=type$UCSC_RefGene_Group,col=list("Gene Structure"=genegroup_col,"CpG Island"=cgi_col)) 

# Column annotation on top
ha_column=HeatmapAnnotation(border = TRUE, "Subtype"=df1$SUBTYPE,col=list("Subtype"= col_subtypes))

###### Generate plot
library(circlize)
library(viridis)

heat_color=  colorRamp2(seq(0,1,0.01), viridis(101))

set.seed(12)
heat <- Heatmap(as.matrix(to_plot), name = "DNAm Level", cluster_rows = T, cluster_columns = T,cluster_column_slices=T, left_annotation = left_annot, right_annotation = right_annot,show_row_names = F,
                show_parent_dend_line = FALSE,row_split = split,border = TRUE, row_title = NULL, show_column_names = F, cluster_row_slices = TRUE, top_annotation = ha_column,
                row_gap = unit(4, "mm"),column_split = col_split,col = heat_color,column_gap=unit(4, "mm")) 


#Save heatmap
pdf("Methylation-Heatmap.pdf",18,18)
draw(heat, annotation_legend_side = "right",heatmap_legend_side="bottom",merge_legend=T)
dev.off()

