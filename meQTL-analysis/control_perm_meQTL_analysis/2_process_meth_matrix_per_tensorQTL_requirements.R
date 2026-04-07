library(readr)
library(GenomicRanges)
library(rtracklayer)

#Load probe location and beta matrix / Add aDMS location to methylation matrix (per tensorQTL requirements)
n <- as.numeric(Sys.getenv("n"))
for (i in 1:n) {

    betas <- read_csv(paste0("control_meth_",i,".csv"))
    head(betas)
    location <- read.csv(paste0("control_locations_",i,".csv"))
    head(location)

    location <- location[,-2]
    location <- location[,c('X','CHR','MAPINFO')]
    location$end <- location$MAPINFO

    colnames(betas)[1] <- "geneid" # in this framework, geneid represents probe id
    colnames(location)[1] <- "geneid"

    merged <- merge(location, betas, by="geneid")
    colnames(merged)[c(1:4)] <- c("phenotype_id",'chr', 'start', 'end')
    merged$start <- merged$start - 1
    merged[,c(1:4)] <- merged[,c("chr","start","end","phenotype_id")]
    colnames(merged)[c(1:4)] <- c("chr","start","end","phenotype_id")
    merged$chr <- paste("chr",merged$chr,sep="")

    write.csv(merged,paste0("control_expression_v2_",i,".csv"), row.names=F) # expression file represents methylation file
    
    }