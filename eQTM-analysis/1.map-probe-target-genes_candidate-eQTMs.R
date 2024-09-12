
# Load libraries
library(readr)
library(rtracklayer)

# Load location info for differentially meth. probes 
probe_location <- read.csv("ancestry-probe_location.txt",sep="")

# Load gene positions
transcript_positions <- read_csv("gene_positions.csv")
dim(transcript_positions)

# Add  genomic window for candidate local eQTM associations 
probe_location$upstream_window <- probe_location$s1 -1500000
probe_location$downstream_window <- probe_location$s1 + 1500000

# Format probe and gene locations as bed files
probe_location <- na.omit(probe_location)
probes_bed <- GRanges(probe_location$chr, IRanges(probe_location$upstream_window, probe_location$downstream_window), names= probe_location$geneid)

genes_bed <- GRanges(transcript_positions$chromosome_name,IRanges(transcript_positions$start_position,transcript_positions$end_position),names=transcript_positions$hgnc_symbol)


# Find overlaps between the 2 files (probe-local gene mapping)
hits <- findOverlaps(genes_bed,probes_bed)
genes_bed[queryHits(hits)]
probes_bed[subjectHits(hits)]

genes_hits <- as.data.frame(genes_bed[queryHits(hits)])
probe_hits <- as.data.frame(probes_bed[subjectHits(hits)])

# Create a dataframe with probe-gene pairs (candidate eQTMs)
merged <- cbind(probe_hits,genes_hits)

# Save file
write.csv(merged,"probe-gene-pairs.csv", row.names = F)


