
# This R script has been used to generate the following figures:
### Figures 1D, 4B
### Supp. Figures 1F, 1G, 6B, 6C

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load sample metadata
# Sample metadata table has the following columns:
# ID  tissue  histone replicate bamFilePath peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Initialize empty list to store the peak annotations
sampleL <- list()
samples <- metadata$ID

# Convert the peak data frames into GRanges
for(sample in samples)
{
  
  hist <- subset(peakSummary, ID == sample)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")
  
  histGR <- renameSeqlevels(histGR, mapSeqlevels(seqlevels(histGR), "UCSC"))
  sampleL[[sample]] <- histGR
}

# Annotate the peaks
peakAnno <- lapply(sampleL, function(x){annotatePeak(x, tssRegion=c(-2000, 2000),
                                                     TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")})

# Plot distance of called peaks from TSS
plotDistToTSS(peakAnno) 

# Plot peak overlap with different genomic regions
plotAnnoBar(peakAnno) 

