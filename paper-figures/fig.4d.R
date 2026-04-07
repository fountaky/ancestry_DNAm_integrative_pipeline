library(readr)
library(dplyr)
library(network)
library(ggplot2)
library(ggnetwork)

# Load all significant eQTMs (probe-gene pairs)
eqtms <- read_csv("sig_eQTMs.csv")

# Load a list with all probes and genes in the HLA eQTM community
comm <- read_csv("HLA_eQTM_community_networkx.txt",
                 col_names = FALSE)

# Across all significant eQTMs, select only genes and probes that are part of the largest HLA eQTM community
data1 <- eqtms[which(eqtms$gene %in% comm),]
data2 <- data1[which(data1$probe %in% comm),]

#  Generate a matrix indicating eQTM probe-gene pairs (0 --> no eQTM, 1 --> eQTM)
matrix <- table(unique(data2[,c("probe","gene")]))
matrix <- (matrix > 0) + 0

plot <- as.matrix(t(matrix))

# Create R network for the largest community 
net = network(plot, directed = FALSE)

# Color probes and genes differently
net %v% "variable" = ifelse(c(row.names(plot),colnames(plot)) %in% data2$probe, "Probe", "Gene")

# Plot R network
n = ggnetwork(net)

p2 <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend))+
  geom_nodes(aes(x, y,color = variable), size = 5) +
  geom_edges(color = "gray35", alpha = 0.15) +
  theme_blank()+
  geom_nodelabel_repel(aes(label = vertex.names),size = 5)


p2
ggsave("HLA-eQTM-community-R-network.jpeg", width = 24, height = 12, device='jpeg', dpi=700)
