
# This R script has been used to generate the following figures:
### Supp. Figures 4C, 7A

library(biomaRt)
library(ChIPpeakAnno)
library(ChIPseeker)
library(chromVAR)
library(dplyr)
library(ezRun)
library(GenomicRanges)
library(ggpubr)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load tissue and histone mark specific count matrix 
# Count matrix was generated using CnT_peakQuantification.R script 
rawCounts <- readRDS("../rawPeakCounts.rds")

# List sample replicates
samples <- colnames(rawCounts)

# List all histone marks
histMarks <- unique(gsub("_[1|2]", "", samples))

# Make GRanges object
rawPeakCounts <- data.frame(rawCounts)

# Extract the genomic region info
rawPeakCounts$loci <- rownames(rawPeakCounts)

rawPeakCounts_GR <- tidyr::separate(rawPeakCounts, "loci", c("chr", "pos"), sep=":")
rawPeakCounts_GR <- separate(rawPeakCounts_GR, "pos", c("start", "end"), sep="-")

rawPeakCounts_GR <- makeGRangesFromDataFrame(rawPeakCounts_GR,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("chromosome", "chrom",
                                                              "chr", "chromosome_name"),
                                             start.field="start",
                                             end.field=c("end", "stop"),
                                             strand.field="strand",
                                             starts.in.df.are.0based=FALSE)

rawPeakCounts_GR <- renameSeqlevels(rawPeakCounts_GR, mapSeqlevels(seqlevels(rawPeakCounts_GR), "UCSC"))

# Load distal enhancer like sequences (dELS) tracks from ENCODE
dELS <- read.table("dELS.bed", header = TRUE, sep = "\t") 

# Convert dELS regions into GRanges object
dELS_GR <- makeGRangesFromDataFrame(dELS,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3")

# Find peak regions overlapping with dELS
dELS_overlapsGR <- GenomicRanges::findOverlaps(rawPeakCounts_GR, dELS_GR, minoverlap = 75L)
dELS_ol         <- data.frame(rawPeakCounts_GR[unique(dELS_overlapsGR@from)])  
dELS_ol$loci    <- paste0(dELS_ol$seqnames, ":", dELS_ol$start, "-", dELS_ol$end)

# Convert the peak overlapping dELS regions into GRanges object
dELS_ol_GR <- GRanges(seqnames = dELS_ol$seqnames, IRanges(start = dELS_ol$start, end = dELS_ol$end), strand = "*", keep.extra.columns=TRUE)

# Annotate the peak overlapping dELS region
peakAnno <- annotatePeak(dELS_ol_GR, tssRegion=c(-2000, 2000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert the overlapping peak annotations into data frame
peakAnnoTab <- data.frame(peakAnno)
peakAnnoTab$loci <- paste0(peakAnnoTab$seqnames, ":", peakAnnoTab$start, "-", peakAnnoTab$end)

peakAnnoTab <- merge(peakAnnoTab, dELS_ol, by="loci")

# Select only dELS overlapping peak regions which can be linked to the closest genes
peak_dELS <- subset(peakAnnoTab, abs(distanceToTSS) > 2000 & abs(distanceToTSS) < 100000)

# Summarize peak counts (overlapping with promoters) at gene level
peak_dELS <- peak_dELS[,c("geneId", samples)]

peak_dELS_rawGeneCounts <- peak_dELS %>% 
  group_by(geneId) %>% 
  summarise(across(everything(), sum))

# Average the raw counts across replicates 
for (i in 1:length(histMarks))
{
  peak_dELS_rawGeneCounts[,histMarks[i]] <- rowMeans(peak_dELS_rawGeneCounts[,grep(histMarks[i], colnames(peak_dELS_rawGeneCounts))])
}

peak_dELS_rawGeneCounts <- column_to_rownames(peak_dELS_rawGeneCounts, "geneId")

# Select only the histone mark wise average raw counts
rawPeakCountsbyHist <- peak_dELS_rawGeneCounts[,histMarks]

# Log normalization
peak_dELS_normGeneCounts <- data.frame(cpm(rawPeakCountsbyHist, log=TRUE))

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

# Combine normalized peak (overlapping with promoters) counts with normalized RNAseq counts
peak_dELS_normGeneCounts <- merge(geneAnnotation, peak_dELS_normGeneCounts, by.x="entrezgene_id", by.y="row.names")
normCounts_comb <- merge(peak_dELS_normGeneCounts, normGeneCounts, by.x="ensembl_gene_id", by.y="row.names") 

# Correlation scatter plot
ggscatter(normCounts_comb, x = "normCounts", y = "H3K18la", 
          cor.coef = TRUE, cor.method = "spearman", color="#2E8B57") +
  labs(x="\n Gene expression \n (log normalized counts)", 
       y="H3K18la peaks in promoters \n (log normalized counts) \n") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 