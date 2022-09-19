# This R script has been used to generate the following figures:
### Figure 1B
### Supp. Figures 1D, 6A

library(edgeR)
library(GenomicAlignments)
library(GenomicRanges)
library(Matrix)
library(matrixStats)

# load binned mouse genome mm10 coordinates
binned_genome <- read.delim("mm10_genes.bed",header=F)
colnames(binned_genome) <- c('Chromosome','Start','End')
binned_genome<-binned_genome[-grep('CHR',binned_genome$Chromosome),]
binned_genome<-binned_genome[-grep('GL',binned_genome$Chromosome),]
binned_genome<-binned_genome[-grep('JH',binned_genome$Chromosome),]

binned_chromosome_ggranges <- GRanges(seqnames = binned_genome$Chromosome,
                                      IRanges(start = binned_genome$Start-2000,
                                              end = binned_genome$Start+2000))
binned_genome<-binned_genome[-grep('MT',binned_genome$Chromosome),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$Chromosome,
                                        IRanges(start = binned_genome$Start,
                                                end = binned_genome$End))

# Load sample metadata
# Sample metadata table has the following columns:
# ID  tissue  histone replicate bamFilePath peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# Count reads within bins for all bam files
mouse_bins_counts <- matrix(data = NA, nrow=length(binned_chromosome_ggranges),ncol=1)

for (i in c(1:length(metadata$bamPath)))
{
  bam.oi<-metadata$bamPath[[i]]
  fragment_counts_mouse <- summarizeOverlaps(binned_chromosome_ggranges,bam.oi,inter.feature=F,ignore.strand=T) %>% assays() %>% .[[1]]
  mouse_bins_counts<- cbind(mouse_bins_counts,fragment_counts_mouse)
  colnames(mouse_bins_counts)[dim(mouse_bins_counts)[2]]<-bam.oi
}

mouse_bins_counts<-mouse_bins_counts[,2:dim(mouse_bins_counts)[2]]
colnames(mouse_bins_counts) <- metadata$SampleID
all_samples_names<-colnames(mouse_bins_counts) 

# Normalise raw counts for multi-dimensional scaling (MDS)
libSize <- Matrix::colSums(mouse_bins_counts)
rpm <- sweep(mouse_bins_counts * 1e6, MARGIN=2, STATS=libSize, FUN="/")
normCounts <- log2(rpm+1)

histone<-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))]
condition <-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))-1]

# Perform MDS
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
plotMDS(normCounts, col=c("blue","orange","darkred","darkgreen")[as.numeric(as.factor(histone))],  cex = 1, pch = c(1,2,3,4,5,6,7,8)[as.numeric(as.factor(condition))])
legend(par("usr")[2], par("usr")[4], legend=levels(as.factor(condition)), pch=c(1,2,3,4,5,6,7,8), col="black", ncol=1, cex=0.6)


