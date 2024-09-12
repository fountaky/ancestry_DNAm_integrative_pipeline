# Load libraries
library(readr)
library(ggplot2)

# Import methylation matrix with beta values across differentialy methylated probes for AFR and EUR patients
data <- read_csv("brca_afr_eur_only.csv")
data <- as.data.frame(data)
row.names(data) <- data[,1]
data <- data[,-1]

data <- data[-which(grepl("rs",row.names(data))),]
na_omit <- na.omit(data)

# Create appropriate data format
t_data <- as.data.frame(t(na_omit))

# Run PCA
pc <- prcomp(t_data,
             center = TRUE,
             scale. = TRUE)


attributes(pc)
summary(pc)

# Import demographics
FinalFullTCGAmeta <- read_csv("FinalFullTCGAmeta2.csv") 
FinalFullTCGAmeta <- as.data.frame(FinalFullTCGAmeta)
colnames(FinalFullTCGAmeta)[3] <- "sample"

row.names(FinalFullTCGAmeta) <- FinalFullTCGAmeta$sample

# Create ancestry for all samples
to_plot <- merge(t_data, FinalFullTCGAmeta[,c("sample","consensus_ancestry","SUBTYPE")], by=0)
to_plot <- to_plot[,-c(dim(to_plot)[2]-2)]
row.names(to_plot) <- to_plot[,1]
to_plot <- to_plot[,-1]

# Generate scatter plot matrix with multiple
multi <- merge(pc$x,to_plot[,c("SUBTYPE","consensus_ancestry")],by=0)

row.names(multi) <- multi[,1]
multi <- multi[,-1]

my_cols <- c("orange","royalblue")

multi$consensus_ancestry <- as.factor(multi$consensus_ancestry)

pdf("PC-Scatterplot-Matrix-bAFR-probes.pdf",width=10, height=8)
pairs(multi[,1:4], pch = 19,  cex = 0.5,
      col = my_cols[multi$consensus_ancestry], upper.panel = NULL)
par(xpd = TRUE)
legend("top",legend=c(levels(multi$consensus_ancestry)),
       fill=my_cols,horiz=T)
dev.off()
