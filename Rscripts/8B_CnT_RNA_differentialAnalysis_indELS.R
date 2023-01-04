# This R script has been used to generate the following figures:
### Figure 3E
### Supp. Figure 5E

library(biomaRt)
library(dplyr)
library(edgeR)
library(ezRun)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(pheatmap)
library(tibble)
library(tidyverse)
library(ChIPseeker)
library(ChIPpeakAnno)
library(RColorBrewer)
library(scales)
library(openxlsx)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

### CnT: Differential expression analysis based on peak counts in promoters
# Load raw peak counts
rawCounts <- readRDS("~/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/MB_MT_H3K18la_rawPeakCounts.rds")

# Load sample metadata
metadata <- read.table("~/public/_Projects/MS_H3K18la/Mouse/metadata.txt", header = TRUE, sep = "\t")

# Select only info for the specific histone 
metadata <- subset(metadata, histone == "H3K18la" & replicate != "merged" & tissue %in% c("MB", "MT"))

# Extract sample info
samples <- colnames(rawCounts)
condition <- as.factor(metadata$tissue)
histone <- as.factor(metadata$histone)

# Add sample information
group <- data.frame(condition, row.names = samples)

# Load enhancer tracks
cCRE <- read.table("cCRE.bed", header = TRUE, sep = "\t") 

dELS <- subset(cCRE, V6 %in% c("dELS", "dELS,CTCF-bound"))

# Subset for proximal and distal enhancers
dELS_GR <- makeGRangesFromDataFrame(dELS,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3")

# Make GRanges of the peak regions from count table
rawCounts <- data.frame(rawCounts)
rawCounts$loci <- rownames(rawCounts)

peakRawCounts <- separate(rawCounts, "loci", c("chr", "pos"), sep=":")
peakRawCounts <- separate(peakRawCounts, "pos", c("start", "end"), sep="-")
peakRawCounts$chr <- gsub("chr", "", peakRawCounts$chr)

peakRawCounts <- makeGRangesFromDataFrame(peakRawCounts,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          seqinfo=NULL,
                                          seqnames.field=c("chromosome", "chrom",
                                                           "chr", "chromosome_name"),
                                          start.field="start",
                                          end.field=c("end", "stop"),
                                          strand="*",
                                          starts.in.df.are.0based=FALSE)

peakRawCounts <- renameSeqlevels(peakRawCounts, mapSeqlevels(seqlevels(peakRawCounts), "UCSC"))

# Find peak regions overlapping with dELS
table(!is.na(findOverlaps(peakRawCounts, dELS_GR, select="arbitrary", minoverlap = 75L)))

dELS_overlapsGR <- GenomicRanges::findOverlaps(peakRawCounts, dELS_GR, minoverlap = 75L)
dELS_ol <- data.frame(peakRawCounts[unique(dELS_overlapsGR@from)])  
dELS_ol$loci <- paste0(dELS_ol$seqnames, ":", dELS_ol$start, "-", dELS_ol$end)

dELS_ol_GR <- GRanges(seqnames = dELS_ol$seqnames, IRanges(start = dELS_ol$start, end = dELS_ol$end), strand = "*", keep.extra.columns=TRUE)

# Annotate peaks
peakAnno <- annotatePeak(dELS_ol_GR, tssRegion=c(-2000, 2000), 
                         TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

# Convert the peak annotations into data frame
peakAnno_df <- data.frame(peakAnno)
peakAnno_df$loci <- paste0(peakAnno_df$seqnames, ":", dELS_ol$start, "-", dELS_ol$end)

peakAnno_df <- merge(peakAnno_df, dELS_ol, by="loci")

# Summarize peak read counts in enhancer
peak_dELS <- subset(peakAnno_df, abs(distanceToTSS) > 2000 & abs(distanceToTSS) < 100000)
peak_dELS <- peak_dELS[-grep("predicted gene|microRNA|RIKEN|Riken|pseudogene", peak_dELS$GENENAME),]

peak_dELS <- peak_dELS[,c("geneId", samples)]

peakPr_rawCounts <- peak_dELS %>% 
  group_by(geneId) %>% 
  summarise(across(everything(), sum))

peakPr_rawCounts <- column_to_rownames(peakPr_rawCounts, "geneId")

# Create DGEList object
y <- DGEList(counts=peakPr_rawCounts, group=condition) 

# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

# Construct design matrix
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(MT - MB, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table                         

# Annotate genes based on log2FC / adjusted p-value thresholds
res$regulation <- ""
res[which(res$logFC >= 0.5 & res$FDR <= 0.05),"regulation"] <- "up"
res[which(res$logFC <= -0.5 & res$FDR <= 0.05),"regulation"] <- "down"
res[which(abs(res$logFC) < 0.5 | res$FDR > 0.05),"regulation"] <- "non-significant"

# Select significant DEGs based on FDR corrected p values
sigGenes <- subset(res, FDR < 0.05 & abs(logFC) > 0.5)

# Add gene names
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("entrezgene_id", "external_gene_name"), 
                        filters="entrezgene_id", values=rownames(sigGenes), mart=ensembl)

sigGenesPr <- merge(geneAnnotation, sigGenes, by.x="entrezgene_id", by.y="row.names")

# Combine with normalized counts 
sigGenesPr <- merge(sigGenesPr, normCounts, by.x="entrezgene_id", by.y="row.names")

# Split significant DEGs into up and down
upGenes <- subset(sigGenesPr, regulation == "up")
downGenes <- subset(sigGenesPr, regulation == "down")

######################################################################################################################################

### RNA: Differential expression analysis
# Load raw gene counts
rawCounts <- readRDS("~/public/_Projects/MS_H3K18la/Mouse/RNA/MB_MT_rawCounts.rds")

rawCounts <- rawCounts[,c(1,3,4,5)]

# Remove genes with no count across all samples
rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]

# Define sample and conditions
samples   <- colnames(rawCounts)
condition <- factor(c(rep("MT",2), rep("MB",2)))

# Add sample information
group <- data.frame(condition, row.names = samples)

# Create DGEList object
y <- DGEList(counts=rawCounts, group=condition)

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- cpm(y, log=TRUE)

# Construct design matrix
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                  

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(MT - MB, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract MT vs MB DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table               

# Gene mapping
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                        filters="ensembl_gene_id", 
                        values=rownames(res), mart=ensembl)

# Combine DE results with gene names
comb <- merge(geneAnnotation, data.frame(res), 
              by.x="ensembl_gene_id", 
              by.y="row.names")

# Rename columns
colnames(comb)[1:2] <- c("geneID", "geneName")

# Combine the DE results with normalized counts
comb <- merge(comb, normCounts, by.x="geneID", by.y="row.names")

# Select significant DE genes 
sigDEG <- subset(comb, PValue < 0.01 & abs(logFC) > 0.05)

######################################################################################################################################

### Combine CnT/RNA differential analyses
# Keep only the overlapping genes from the RNA and CnT differential analysis
DEG_overlaps <- sigDEG[which(sigDEG$geneName %in% sigGenesPr$external_gene_name),]

# Keep only relevant columns
multiDEG <- sigDEG[,c("geneName","PValue","FDR","logFC")]

# Combine log2FC from CnT and RNAseq
combDE <- merge(sigGenesPr, multiDEG, by.x="external_gene_name", by.y="geneName")

# Reorder columns
combDE <- combDE[,c("entrezgene_id", "external_gene_name", 
                    "logFC.x", "PValue.x", "FDR.x",
                    "logFC.y", "PValue.y", "FDR.y")]

# Update the column names
colnames(combDE) <- c("Entrez_geneID", "Gene",
                      "logFC_H3K18la", "pVal_H3K18la", "FDR_H3K18la",
                      "logFC_RNA", "pVal_RNA", "FDR_RNA")

# Categorize the log2FC trends 
combDE$regulation <- "opposite"
combDE[which(combDE$logFC_H3K18la > 0 & combDE$logFC_RNA > 0),"regulation"] <- "up"
combDE[which(combDE$logFC_H3K18la < 0 & combDE$logFC_RNA < 0),"regulation"] <- "down"

combDE <- subset(combDE, abs(combDE$logFC_RNA) > 0.5 & combDE$pVal_RNA < 0.05)

# Multiomics correlation scatter plot
ggscatter(combDE, x = "logFC_H3K18la", y = "logFC_RNA", 
          label = "Gene", font.label = c(10, "bold.italic"),
          color = "regulation", palette = c("#00AFBB", "#2E8B57", "#FC4E07")) +
  stat_cor(aes(label = ..r.label..)) +
  labs(x="log2FC H3K18la", y="log2FC RNAseq") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) +
  geom_vline(xintercept = 0, colour = "#323232") + geom_hline(yintercept = 0, colour = "#323232")

