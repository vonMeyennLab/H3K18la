library(chromVAR)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(ggVennDiagram)
library(ggplot2)
library(GenomicAlignments)
library(DESeq2)
library(edgeR)
library(ggcorrplot)
library(ggfortify)
library(clusterProfiler)

########################################################################################
###################### Quantify bams ##################################################
########################################################################################

#################
#### Fig 1B ##### 
#################

# Create reference peak list

# Load sample metadata
metadata <- read.table("~/metadata_mm.txt", header = TRUE, sep = "\t")

# List merged peak files
mergedPeaks <- metadata$Peaks

# Merge all called peaks for all samples
mPeak <- GRanges()

for(mp in mergedPeaks)
{
  peakRes <- read.table(mp, header = FALSE, fill = TRUE)
  mPeak   <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}

# Reduce genomic ranges
masterPeak <- GenomicRanges::reduce(mPeak)

# List samples
samples <- metadata$SampleID

# Initialize count matrix
countMat <- matrix(NA, length(masterPeak), length(samples))

# Overlap with bam files to get count table
for(i in (1:nrow(metadata)))
{
  bam <- metadata[i, "bamFiles"]
  fragment_counts <- chromVAR::getCounts(bam, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
  countMat[, i]   <- counts(fragment_counts)[,1]
}

# Add sample names
colnames(countMat) <- samples

# Add peak regions
df <- data.frame(mPeak)
rownames(countMat) <- paste0(df$seqnames, ":", df$start, "-", df$end)

# Remove peaks with no count across all samples
rawCounts <- countMat[rowSums(countMat) > 0, ]

# Add sample information
condition <- as.factor(metadata$Tissue)
histone <- as.factor(metadata$Factor)
group <- data.frame(condition, histone, row.names = samples)

# Create DGEList object
y <- DGEList(counts=rawCounts, group=condition)

# Remove peak regions with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Apply TMM normalization
y <- calcNormFactors(y)

# Extract normalized counts
normCounts <- edgeR::cpm(y, log=TRUE)

# Pearson correlation plot 
pheatmap(
  mat               = cor(normCounts, use="complete.obs", method = "pearson"),
  silent            = F,
  cutree_rows       = 4,
  cutree_cols       = 4,
  annotation_col    = group,
  color             = brewer.pal(n = 9, name = "YlOrRd"),
  fontsize_row      = 8, 
  fontsize_col      = 8,
  display_numbers   = FALSE)

#################
#### Fig S1A #### 
#################

# load binned genome
binned_genome <- read.delim("~/binned_genome/binned_genome.mm10.bed")
binned_genome<-binned_genome[-grep('MT',binned_genome$Chromosome),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$Chromosome,
                         IRanges(start = binned_genome$Start,
                                 end = binned_genome$End))

# read metadata with bamPaths for all files
metadata <- read.delim("~/metadata.txt")

# count reads within bins for all bams
mouse_bins_counts <- matrix(data = NA, nrow=length(binned_chromosome_ggranges),ncol=1)
for (i in c(1:length(metadata$bamPath))){
bam.oi<-metadata$bamPath[[i]]
fragment_counts_mouse <- summarizeOverlaps(binned_chromosome_ggranges,bam.oi,inter.feature=F,ignore.strand=T) %>% assays() %>% .[[1]]
mouse_bins_counts<- cbind(mouse_bins_counts,fragment_counts_mouse)
colnames(mouse_bins_counts)[dim(mouse_bins_counts)[2]]<-bam.oi
}
mouse_bins_counts<-mouse_bins_counts[,2:dim(mouse_bins_counts)[2]]

all_samples_names <- list.files(path="~/peakFiles",pattern=".bed")

colnames(mouse_bins_counts) <- gsub(".bed","",all_samples_names)[-grep('merged',all_samples_names)]

# normalise raw counts for MDS
libSize <- Matrix::colSums(mouse_bins_counts)
rpm <- sweep(mouse_bins_counts * 1e6, MARGIN=2, STATS=libSize, FUN="/")
normCounts <- log2(rpm+1)

histone<-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))]
condition <-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))-1]

#plot MDS
plotMDS(normCounts, col=c("blue","orange","darkred","darkgreen")[as.numeric(as.factor(histone[-grep('merged',all_samples_names)]))],  
        cex = 1, pch = c(1,2,3,4,5,6)[as.numeric(as.factor(condition[-grep('merged',all_samples_names)]))])


########################################################################################
###################### From here: work with peaks ######################################
########################################################################################

#################
#### Fig 1C #####
#### Fig S1E ####
#### Fig S1F ####
#################

# Read all peaks used in this project
peak.list <- list.files(path="~/peakFiles",pattern=".bed")
peaks<-list()
for (i in 1:length(peak.list)){
peaks[[i]] <- read.delim(peak.list[i])
colnames(peaks[[i]])<-c('peak.chr','peak.start','peak.end','peak.total.value','peak.max.value','peak.max.value.region')}

# For table S1 : QC metrics: Calculate peak number and peak width for each sample
peakN       <- c()    # peak number
peakSummary <- c()    # peak summary

for (i in 1:length(peak.list)){
sample=peak.list[i]
  sampleInfo <- strsplit(sample, "_")[[1]]
  peakInfo   <- peaks[[i]]
  
  # Calculate peak width
  peakInfo$width      <- abs(peakInfo$peak.end - peakInfo$peak.start)
  peakInfo$histone    <- sampleInfo[2]
  peakInfo$biology  <- sampleInfo[1]
  peakInfo$rep  <- sampleInfo[3]
  
  peakSummary <- rbind(peakSummary, peakInfo)
  peakN       <- data.frame(peakN = nrow(peakInfo), 
                            peakTotalBases = sum(peakInfo$width), 
                            histone = sampleInfo[2], 
                            biology = sampleInfo[1] , 
                            rep = sampleInfo[3]) %>% rbind(peakN, .)
}

# Annotate peaks to genomic regions

# Get TSS regions
promoter_mm <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=3000, downstream=3000)

#convert peak dataframes to ggranges objects
peaks.ggranges<-list()
for(i in 1:length(peaks)){
  if(length(grep('chrUn',peaks[[i]]$peak.chr))>0){peaks[[i]]<-peaks[[i]][-grep('chrUn',peaks[[i]]$peak.chr),]}
  if(length(grep('random',peaks[[i]]$peak.chr))>0){peaks[[i]]<-peaks[[i]][-grep('random',peaks[[i]]$peak.chr),]}
  peaks.ggranges[[i]]<- makeGRangesFromDataFrame(peaks[[i]],
                                                 keep.extra.columns=TRUE,
                                                 ignore.strand=FALSE,
                                                 seqinfo=NULL,
                                                 seqnames.field=c("chromosome", "chrom",
                                                                  "chr", "chromosome_name"),
                                                 start.field="start",
                                                 end.field=c("end", "stop"),
                                                 strand.field="strand",
                                                 starts.in.df.are.0based=FALSE) }
names(peaks.ggranges) <- names(peaks)

# annotate peaks to genomic elements with ChIPseeker
peak_anno.mouse<-list()
for (i in 1:length(peaks.ggranges)){
  peak_anno.mouse[[i]] <- annotatePeak(peaks.ggranges[[i]],
                                  tssRegion=c(-2000, 2000),
                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                  annoDb='org.Mm.eg.db')
}

names(peak_anno.mouse) <- names(peaks)

#plot genomic distribution of peaks
subset_peak_anno.mouse <- peak_anno.mouse[-grep('merged',names(peak_anno.mouse))]
PTM.order <- unlist(strsplit(names(subset_peak_anno.mouse),'_'))
PTM.order <- PTM.order[grep('H3',PTM.order)]
subset_peak_anno.mouse<-subset_peak_anno.mouse[order(PTM.order)]

plotAnnoBar(subset_peak_anno.mouse[c(2,3,4,5,9,10,6,7,8,1,11,22,23,24,25,26,27,13,14,15,16,17,18,19,12)])
plotDistToTSS(subset_peak_anno.mouse[c(2,3,4,5,9,10,6,7,8,1,11,22,23,24,25,26,27,13,14,15,16,17,18,19,12)],
              title = 'Distribution relative to TSS')

subset_peak_anno.mouse <- peak_anno.mouse[grep('merged',names(peak_anno.mouse))]
PTM.order <- unlist(strsplit(names(subset_peak_anno.mouse),'_'))
PTM.order <- PTM.order[grep('H3',PTM.order)]
subset_peak_anno.mouse<-subset_peak_anno.mouse[order(PTM.order)]

plotAnnoBar(subset_peak_anno.mouse)
plotDistToTSS(subset_peak_anno.mouse, title = 'Distribution relative to TSS')

#################
#### Fig 1D ##### 
#################

# Initialize empty list
sampleFE <- list()

# Calculating genome size based on gene regions and intergenic regions
genome_size <- 
  sum(
    GenomicRanges::reduce(genome_annotation$`Gene bodies`, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum,
    GenomicRanges::reduce(genome_annotation$`Intergenic regions`, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
  )

calc_fold_enrichment <- 
  function(gr1, gr2, genome_size){
    
    A <- 
      GenomicRanges::intersect(
        GenomicRanges::reduce(gr2, ignore.strand = TRUE),
        GenomicRanges::reduce(gr1, ignore.strand = TRUE)
      ) %>% GenomicRanges::width() %>% sum
    
    B <- 
      GenomicRanges::reduce(gr1, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
    
    C <- 
      GenomicRanges::reduce(gr2, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum
    
    log2((A/B)/(C/genome_size))
  }

# Compute peak fold enrichment for different genomic features
for (sample in samples)
{
  hist <- subset(peakSummary, sampleID == sample)
  hist$chr <- gsub("chr", "", hist$chr)
  histGR <- GRanges(hist$chr, IRanges(start = hist$start, end = hist$end), strand = "*")
  
  histFE <- lapply(genome_annotation, function(x){calc_fold_enrichment(histGR, x, genome_size = genome_size)})
  sampleFE[[sample]] <- histFE
}

sampleFE <- lapply(sampleFE, function(x){x <- do.call("rbind", x)})

# Combine the fold enrichment values
comb_FE <- do.call("cbind", sampleFE)
colnames(comb_FE) <- names(sampleFE)

breaksList <- seq(-2, 2, by = 0.1)

# Peak fold enrichment heat map
pheatmap::pheatmap(t(comb_FE[c(1:2,5:6,3,7),tissueOrder]),
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "row",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(length(breaksList)),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 40, 
                   cellwidth         = 40,
                   border_color      = 'black',
                   breaks            = breaksList)


#################
#### Fig S1G #### 
#################

# Quantify the promoters

# Load promoter regions
promoter <- genome_annotation$promoters
  
# Initialize count matrix
countMatPr <- matrix(NA, length(promoters), length(samples))

# Overlap with bam files to get count table
for(i in (1:nrow(metadata)))
{
  bam <- metadata[i, "bamFiles"]
  fragment_counts <- chromVAR::getCounts(bam, promoters, paired = TRUE, by_rg = FALSE, format = "bam")
  countMatPr[, i]   <- counts(fragment_counts)[,1]
}

# Add sample names
colnames(countMatPr) <- samples

# Add promoter regions
dfPr <- data.frame(promoters)
rownames(countMatPr) <- paste0(dfPr$seqnames, ":", dfPr$start, "-", dfPr$end)

# Perform MDS on promoter peaks specifically:

# Remove promoters with no count across all samples
rawCountsPr <- countMatPr[rowSums(countMatPr) > 0, ]

# normalise raw counts for MDS
libSize <- Matrix::colSums(rawCountsPr)
rpm <- sweep(rawCountsPr * 1e6, MARGIN=2, STATS=libSize, FUN="/")
normCountsPr <- log2(rpm+1)

plotMDS(normCountsPr,col=c("blue","orange","darkred","darkgreen")[as.numeric(as.factor(histone[-grep('merged',all_samples_names)]))],  
        cex = 1, pch = c(1,2,3,4,5,6)[as.numeric(as.factor(condition[-grep('merged',all_samples_names)]))])


#################
#### Fig 1E #####
#### Fig 1F ##### 
#################

# Promoter peak correlation to gene expression

# List all histone marks
histMarks <- unique(gsub("_[1|2]", "", colnames(rawCountsPr)))

# Add loci for promoter regions
promoters$loci <- paste0(promoters$seqnames, ":", promoters$start, "-", promoters$end)

rawPeakCounts <- data.frame(rawCountsPr[promoters$loci,])

rawPeakCountsPr <- merge(promoters, rawPeakCounts, by.x="loci", by.y="row.names")

# Summarize peak read counts in promoter
# Subset by tissue
peakPr <- rawPeakCountsPr[,c("gene_id", "MT_H3K18la_1", "MT_H3K18la_2")]

peakPr_rawGeneCounts <- peakPr %>% 
  group_by(gene_id) %>% 
  summarise(across(everything(), sum))

# Average raw counts across biological replicates
for (i in 1:length(histMarks))
{
  peakPr_rawGeneCounts[,histMarks[i]] <- rowMeans(peakPr_rawGeneCounts[,grep(histMarks[i], colnames(peakPr_rawGeneCounts))])
}

peakPr_rawGeneCounts <- column_to_rownames(peakPr_rawGeneCounts, "gene_id")

# Select only the mark wise peak raw counts
rawPeakCountsbyHist <- peakPr_rawGeneCounts[histMarks]

# Log normalization
peakPr_normGeneCounts <- data.frame(cpm(rawPeakCountsbyHist, log=TRUE))

# Load raw RNAseq data
rawCounts_RNAseq <- read.table("~/rawCounts.txt", header = TRUE, sep = "\t", check.names = FALSE) 

# Average raw counts across biologicl replicates
rawCounts_RNAseq$RNA <- rowMeans(rawCounts_RNAseq) 


# Remove genes with no count across all samples
rawCounts_RNAseq <- rawCounts_RNAseq[rowSums(rawCounts_RNAseq) > 0, ]

# log normalization
normCounts_RNA <- data.frame(row.names = rownames(rawCounts_RNAseq), 
                             normCounts = cpm(rawCounts_RNAseq$RNA, log=TRUE))

# Add gene names
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("entrezgene_id", "ensembl_gene_id", "external_gene_name"), 
                        filters="entrezgene_id", values=rownames(peakPr_normGeneCounts), mart=ensembl)

peakPr_normGeneCounts <- merge(geneAnnotation, peakPr_normGeneCounts, by.x="entrezgene_id", by.y="row.names")

# Correlation scatter plot
normCounts_comb <- merge(peakPr_normGeneCounts, normCounts_RNA, by.x="ensembl_gene_id", by.y="row.names")

# Extract lists of high and low expressed genes
normCounts_sub <- normCounts_comb[,c(1,4:7)]

normCounts_sub <- normCounts_sub %>% 
  group_by(EnsemblID) %>% 
  summarise(across(everything(), mean))
normCounts_sub <- column_to_rownames(normCounts_sub, "EnsemblID")

# Annotate low expressed genes
normCounts_sub$geneExpression <- ""
normCounts_sub <- normCounts_sub[order(normCounts_sub$RNA),]
normCounts_sub[1:2000, "geneExpression"] <- "Low"

# Annotate high expressed genes
normCounts_sub <- normCounts_sub[order(-normCounts_sub$RNA),]
normCounts_sub[1:2000, "geneExpression"] <- "High"

# Filer only for high and low expressed genes
normCounts_sub <- subset(normCounts_sub, geneExpression %in% c("Low", "High"))
normCounts_sub$geneExpression <- factor(normCounts_sub$geneExpression)

normCounts_sub <- rownames_to_column(normCounts_sub, "EnsemblID")
normCounts_sub <- merge(normCounts_sub, geneAnnotation, by.x= "EnsemblID", by.y="ensembl_gene_id")
normCounts_sub <- normCounts_sub[,c(1,7:8,2:6)]
colnames(normCounts_sub) <- c("EnsemblID", "EntrezID", "GeneName", 
                              "H3K18la", "H3K27me3", "H3K4me3",
                              "RNA", "GeneExp")

# Add colors
my_cols <- c("#FC4E07", "#00AFBB")

# Pair wise peak correlation scatter plot
ggscatter(normCounts_sub, x = "H3K4me3", y = "H3K27me3", 
          cor.coef = TRUE, cor.method = "pearson",  color = "GeneExp", palette = my_cols) +
  labs(x="\n H3K4me3 peaks in promoters \n (log normalized counts)", 
       y="H3K27me3 peaks in promoters \n (log normalized counts) \n") +
  geom_point(alpha = 0.05) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        axis.text=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.title=element_text(size=11,face="bold"),
        legend.title = element_text(size=10, face="italic")) 

write.table(normCounts_comb, "PIM_hPTM_RNA_Pr_normCounts.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(normCounts_sub, "PIM_hPTM_RNA_Pr_normCounts_top2000.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

ESC <- read.delim("~/ESC_hPTM_RNA_Pr_normCounts.txt")
GAS <- read.delim("~/GAS_hPTM_RNA_Pr_normCounts.txt")
PIM <- read.delim("~/PIM_hPTM_RNA_Pr_normCounts.txt")

ggplot(ESC,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('ESC')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw()
ggplot(GAS,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('GAS')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)")+ theme_bw()
ggplot(PIM,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('PIM')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)")+ theme_bw()

cor.test(PIM$H3K18la,PIM$RNA, method='spearman')
cor.test(GAS$H3K18la,GAS$RNA, method='spearman')
cor.test(ESC$H3K18la,ESC$RNA, method='spearman')

# perform GO analysis on highest expressed genes and genes associated to highest lactylated promoters
#   RNA

# load gene expression + promoter lactylation data with annotation of highest 2000 genes in each category
x <- list.files("/home/egalle/public/_Projects/MS_H3K18la/Mouse/GO",pattern='comb')
genes <- list()
for (i in 1:length(x)){
    genes[[i]] <- read.delim(x[i])
}
names(genes)<-x
genes.RNA <- genes[grep('RNA',x)]
genes.H3K18la <- genes[grep('H3K18la',x)]

genes<-list(unlist(genes.RNA[[1]]$EntrezID[genes.RNA[[1]]$GeneExp=='High']),
            unlist(genes.RNA[[2]]$EntrezID[genes.RNA[[2]]$GeneExp=='High']),
            unlist(genes.RNA[[3]]$EntrezID[genes.RNA[[3]]$GeneExp=='High']))

GO.RNA<-list()
for (i in 1:3){
    GO.RNA[[i]] <- enrichGO(unlist(genes[i]),pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             OrgDb='org.Mm.eg.db',
                             ont="ALL")}
names(GO.RNA) <- c('ESC','GAS','PIM')
temp.GO.RNA<-list(GO.RNA[[1]]$Description,GO.RNA[[2]]$Description,GO.RNA[[3]]$Description)
names(temp.GO.RNA) <- c('ESC','GAS','PIM')

dotplot(GO.RNA[[1]],showCategory=15,title='ESC')
dotplot(GO.RNA[[2]],showCategory=15,title='GAS')
dotplot(GO.RNA[[3]],showCategory=15,title='PIM')

genes<-list(unlist(genes.H3K18la[[1]]$EntrezID),unlist(genes.H3K18la[[2]]$EntrezID),unlist(genes.H3K18la[[3]]$EntrezID))
genes<-list(ESC[order(ESC$H3K18la,decreasing=T),]$EntrezID[1:2000],GAS[order(GAS$H3K18la,decreasing=T),]$EntrezID[1:2000],PIM[order(PIM$H3K18la,decreasing=T),]$EntrezID[1:2000])

GO.H3K18la<-list()
for (i in 1:3){
GO.H3K18la[[i]] <- enrichGO(unlist(genes[i]),pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH",
                     OrgDb='org.Mm.eg.db',
                     ont="ALL")}

names(GO.H3K18la) <- c('ESC','GAS','PIM')
dotplot(GO.H3K18la[[1]],showCategory=15,title='ESC')
dotplot(GO.H3K18la[[2]],showCategory=15,title='GAS')
dotplot(GO.H3K18la[[3]],showCategory=15,title='PIM')

# CGI Promoter peak correlation to gene expression
ESC <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_CGIpromoter/ESC_H3K18la_RNA_CGIpr_normCounts.txt")
GAS <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_CGIpromoter/GAS_hPTM_RNA_CGIpr_normCounts.txt")
PIM <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_CGIpromoter/PIM_hPTM_RNA_CGIpr_normCounts.txt")

ggplot(ESC,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('ESC')+ylab("H3K18la CGI promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw()
ggplot(GAS,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('GAS')+ylab("H3K18la CGI promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)")+ theme_bw()
ggplot(PIM,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('PIM')+ylab("H3K18la CGI promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)")+ theme_bw()

cor.test(PIM$H3K18la,PIM$RNA, method='spearman')
cor.test(GAS$H3K18la,GAS$RNA, method='spearman')
cor.test(ESC$H3K18la,ESC$RNA, method='spearman')
