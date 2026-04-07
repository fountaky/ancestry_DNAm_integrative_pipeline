library(readr)
library(dplyr)
library(preprocessCore, lib="~/R.libs")

n <- as.numeric(Sys.getenv("n"))
for (i in 1:n) {

    meth <- read.table(paste0("control_expression_v2_",i,".csv"),sep=",")
    probes <- read.table(paste0("control_expression_v3_",i,".csv"),sep=",")
    
    colnames(meth) <- meth[1,]
    meth <- meth[-1,]
    row.names(meth) <- meth$phenotype_id

    colnames(probes) <- probes[1,]
    probes <- probes[-1,]
    
    row.names(probes) <- probes[,1]
    
    meth <- meth[row.names(probes),]
    dim(meth)
    
    # Load methylation residuals
    residuals <- read.table(paste0('tensor_meth_residuals_3_',i,'.txt'))
    dim(residuals)
    
    row.names(residuals) <- colnames(meth)[c(5:dim(meth)[2])]
    head(residuals)
    
    peer_matrix <- t(residuals)
    peer_matrix <- as.data.frame(peer_matrix)
    
    peer_matrix <- as.data.frame(normalize.quantiles(as.matrix(peer_matrix)))
    
    colnames(peer_matrix) <- colnames(meth)[c(5:dim(meth)[2])]
    head(peer_matrix)
    dim(peer_matrix)
    
    meth_residuals <- cbind(meth[,c(1:4)], peer_matrix)
    dim(meth_residuals)
    
    write.csv(meth_residuals, paste0("methylation_residuals_3_v2_",i,".csv"), row.names=FALSE)
    
    }