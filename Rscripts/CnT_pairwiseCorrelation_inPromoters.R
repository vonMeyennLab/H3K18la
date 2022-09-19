
# This R script has been used to generate the following figures:
### Figures 1G, 4E
### Supp. Figures 2A, 2B, 6H

library(biomaRt)
library(chromVAR)
library(dplyr)
library(ezRun)
library(GenomicRanges)
library(ggpubr)
library(tidyverse)

# Load sample metadata
# Sample metadata table has the following columns:
# ID  tissue  histone replicate bamFilePath peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Select only the replicate specific metadata
# One can also additionally filter for tissue type and/or histone mark
metaRep <- subset(metadata, replicate != "merged")

# Load promoter regions
# Similarly, CGI promoter regions can also be used for quantification
promoter <- read.table("../mm10_promoters.bed", header = TRUE, fill = TRUE)

# Convert promoter regions into GRanges
promoterGR   <- GRanges(seqnames = promoter$chr, IRanges(start = promoter$start, end = promoter$end), strand = "*")

# List sample replicates
samples <- metaRep$ID

# Initialize count matrix
countMat <- matrix(NA, length(promoterGR), length(samples))
dim(countMat)

# Promoter overlap with bam files to get count table
for(i in (1:nrow(metaRep)))
{
  bam <- metaRep[i, "bamFile"]
  fragment_counts <- chromVAR::getCounts(bam, promoterGR, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i]   <- counts(fragment_counts)[,1]
}

# Rename the count matrix with sample replicates
colnames(countMat) <- samples

# Save the count matrix as a data frame
# Convert the promoter regions into unique row identifiers
df <- data.frame(masterPeak)
rownames(countMat) <- paste0(df$seqnames, ":", df$start, "-", df$end)

# List all histone marks
histMarks <- unique(gsub("_[1|2]", "", colnames(countMat)))

# Create unique loci for each promoter region
promoter$loci <- paste0(promoter$seqnames, ":", promoter$start, "-", promoter$end)

# Combine the promoter annotations with the raw counts
rawCounts <- data.frame(countMat[promoter$loci,])
rawCounts_prom <- merge(promoter, countMat, by.x="loci", by.y="row.names")

# Summarize promoter counts at gene level
rawCounts_prom <- rawCounts_prom[,c("gene_id", samples)]

prom_rawGeneCounts <- rawCounts_prom %>% 
  group_by(gene_id) %>% 
  summarise(across(everything(), sum))

# For each tissue, average the raw counts across replicates 
for (i in 1:length(histMarks))
{
  prom_rawGeneCounts[,histMarks[i]] <- rowMeans(prom_rawGeneCounts[,grep(histMarks[i], colnames(prom_rawGeneCounts))])
}

prom_rawGeneCounts <- column_to_rownames(prom_rawGeneCounts, "gene_id")

# Select only the histone mark wise average raw promoter counts
rawPromCountsbyHist <- prom_rawGeneCounts[histMarks]

# Log normalize the counts
prom_normGeneCounts <- data.frame(cpm(rawPromCountsbyHist, log=TRUE))

# Customize the graphics panels
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

# Generate pairwise scatter plots for all histone marks
pairs(peakPr_normGeneCounts, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

########################################################################################################################

# Load normalized (RPKM) RNAseq data
normCounts_RNA <- read.table("../normCounts_RPKM.txt", header = TRUE, sep = "\t")
normCounts_RNA <- column_to_rownames(normCounts_RNA, "Gene_ID")

# Average the RPKM values across replicates
normCounts_RNA$RNA <- rowMeans(normCounts_RNA)

# Remove genes with no counts
normCounts_RNA <- normCounts_RNA[rowSums(normCounts_RNA) > 0, ]

# log normalization
normGeneCounts <- data.frame(row.names = rownames(normCounts_RNA), 
                             normCounts = log2(normCounts_RNA$RNA))

# Add gene names
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", values=rownames(normCounts_RNA), mart=ensembl)

# Combine normalized promoter counts with normalized RNAseq counts
peakPr_normGeneCounts <- merge(geneAnnotation, peakPr_normGeneCounts, by.x="entrezgene_id", by.y="row.names")
normCounts_comb <- merge(peakPr_normGeneCounts, normGeneCounts, by.x="ensembl_gene_id", by.y="row.names") 

# Annotate top 2000 low expressed genes
normCounts_comb$GeneExp <- ""
normCounts_comb <- normCounts_comb[order(normCounts_comb$RNA),]
normCounts_comb[1:2000, "GeneExp"] <- "Low"

# Annotate top 2000 high expressed genes
normCounts_comb <- normCounts_comb[order(-normCounts_comb$RNA),]
normCounts_comb[1:2000, "GeneExp"] <- "High"

# Select only top 2000 high and low expressed genes
normCounts_sub <- subset(normCounts_comb, GeneExp %in% c("High", "Low"))

# Make scatter plot for pairwise histone marks
# Only top 2000 high and low expressed genes are considered
ggscatter(normCounts_sub, x = "H3K18la", y = "H3K4me3", 
                cor.coef = TRUE, cor.method = "pearson", main="ESC.ser") +
  labs(x="\n H3K18la promoters \n (log normalized counts)", 
       y="H3K4me3 promoters \n (log normalized counts) \n") +
  geom_point(alpha = 0.05) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 
