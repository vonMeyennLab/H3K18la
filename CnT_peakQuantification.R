
# This R script has been used to generate the following figures:
### Figures 1C, 4A

library(chromVAR)
library(GenomicRanges)
library(tidyverse)

# Load sample metadata
# Sample metadata table has the following columns:
# ID  tissue  histone replicate bamFilePath peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Select only the metadata for merged peaks
# One can also additionally filter for tissue type and/or histone mark
metaComb <- subset(metadata, replicate == "merged")

# Select only the replicate specific metadata
# One can also additionally filter for tissue type and/or histone mark
metaRep <- subset(metadata, replicate != "merged")

# List merged peaks to be used
mergedPeaks <- metaComb$peakFile

# Combine all selected merged peaks to create a master peak list
mPeak <- GRanges()

for(mp in mergedPeaks)
{
  peakRes <- read.table(mp, header = TRUE, fill = TRUE)
  mPeak   <- GRanges(seqnames = peakRes$chr, IRanges(start = peakRes$start, end = peakRes$end), strand = "*") %>% append(mPeak, .)
}

# Reduce genomic ranges (only for genome wide quantification)
masterPeak <- GenomicRanges::reduce(mPeak)

# List sample replicates
samples <- metaRep$ID

# Initialize count matrix
countMat <- matrix(NA, length(masterPeak), length(samples))
dim(countMat)

# Peak overlap with bam files to get count table
for(i in (1:nrow(metaRep)))
{
  bam <- metaRep[i, "bamFile"]
  fragment_counts <- chromVAR::getCounts(bam, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i]   <- counts(fragment_counts)[,1]
}

# Rename the count matrix with sample replicates
colnames(countMat) <- samples

# Save the count matrix as a data frame
# Convert the genomic regions into unique row identifiers
df <- data.frame(masterPeak)
rownames(countMat) <- paste0(df$seqnames, ":", df$start, "-", df$end)

