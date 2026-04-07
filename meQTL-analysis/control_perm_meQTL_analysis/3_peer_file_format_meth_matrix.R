library(readr)
library(dplyr)

# Further process methylation matrix to obtain input for PEER software
n <- as.numeric(Sys.getenv("n"))
for (i in 1:n) {
    meth <- read.table(paste0("control_expression_v2_",i,".csv"),sep=",")
    head(meth)
    meth <- meth[,-c(1:3)]

    colnames(meth) <- meth[1,]
    meth <- meth[-1,]

    row.names(meth) <- meth[,1]
    meth <- meth[,-1]

    meth <- na.omit(meth)
    dim(meth)

    set.seed(6543)
    meth <- meth[sample(row.names(meth),567),] # After removing probes with missing values (peer requirement), the number of aDMS included in the ancestry meQTL analysis dropped to 567
    dim(meth)

    write.csv(meth, paste0("control_expression_v3_",i,".csv"), row.names=TRUE)

    peer_meth <- t(meth)
    head(peer_meth)
    dim(peer_meth)

    row.names(peer_meth) <- NULL
    colnames(peer_meth) <- NULL

    peer_meth <- as.data.frame(peer_meth)
    peer_meth <- mutate_all(peer_meth, function(x) as.numeric(as.character(x)))
    head(peer_meth)

    dim(peer_meth)

    write.table(peer_meth,paste0("methylation_to_peer_",i,".csv"), row.names = F, col.names = F)
    
    }