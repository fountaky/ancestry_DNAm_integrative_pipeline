library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(readxl)

########## Remove ancestry probes from all probes
control <- read.csv("control_probes_betas.csv")

row.names(control) <- control[,1]
control <- control[,-1]
head(control)

#### Identify and remove X and Y chromosome probes
# Load Illumina annotations
humanmethylation450_15017482_v1_2_1_ <- read_excel("humanmethylation450_15017482_v1-2.xlsx")
annotations <- humanmethylation450_15017482_v1_2_1_[-c(1:6),]
colnames(annotations) <- annotations[1,]
annotations <- annotations[-1,]

locations <- as.data.frame(annotations[,c("IlmnID","Genome_Build","MAPINFO","CHR")])
row.names(locations) <- locations$IlmnID
locations <- locations[,-1]

# Remove sex chromosome probes
sex_chrom <- locations[which(locations$CHR %in% unique(locations$CHR)[c(1:2)]),]

######### Sample 757 random probes with no ancestry associations ###################
control <- control[-which(row.names(control) %in% row.names(sex_chrom)),] # remove sex chromosome probes
control <- control[-which(grepl("rs",row.names(control))),] # remove all "rs" probes

# Load probes with ancestry associated pvalue < 0.5 
ancestry <- read.csv("Ancestry_Effect_Probes.csv")

# Remove ancestry associated probes (pvalue < 0.5) from control probe pool 
control2 <- control[-which(row.names(control) %in% as.character(ancestry[,1])),]

n <- as.numeric(Sys.getenv("n"))
for (i in 1:n){
    
    final <- control2[sample(row.names(control2), 757),]
    loc <- locations[row.names(final),]

    # Save control methylation matrix
    write.csv(final, paste0("control_meqtls_generated_files/control_meth_",i,".csv"), row.names = TRUE)
    write.csv(loc, paste0("control_locations_",i,".csv"), row.names = TRUE)
    }
