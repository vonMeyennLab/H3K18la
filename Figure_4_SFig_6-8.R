library(chromVAR)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(ggVennDiagram)
library(ggplot2)
library(GenomicAlignments)
library(DESeq2)
library(edgeR)
library(ggcorrplot)
library(ggfortify)

# SFig 6A : 
########################################################################################
###################### Counting bams inside bins ######################################
########################################################################################

# load the binned genome and format 
binned_genome <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/binned_genome/binned_genome_hg38.bed.txt")
binned_genome <- data.frame(binned_genome$Chromosome,binned_genome$Start,binned_genome$End)
colnames(binned_genome)<-c('chr','start','end')
binned_genome$chr<-paste0('chr',binned_genome$chr)
#masterPeak_df <- bind_rows(peaks, .id = "column_label")
binned_genome<-binned_genome[-grep('MT',binned_genome$chr),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$chr, 
                         IRanges(start = binned_genome$start, 
                                 end = binned_genome$end))	
								 
# load the metadata for all human samples containing paths to bam files in column 'bamPath'
metadata <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Human/metadata.txt")
	
# count reads within bins for all samples	  
human_bins_counts <- matrix(data = NA, nrow=length(binned_chromosome_ggranges),ncol=1)
for (i in c(1:length(metadata$bamPath))){
bam.oi<-metadata$bamPath[[i]]
fragment_counts_mouse <- summarizeOverlaps(binned_chromosome_ggranges,bam.oi,inter.feature=F,ignore.strand=T) %>% assays() %>% .[[1]]
human_bins_counts<- cbind(human_bins_counts,fragment_counts_mouse)
colnames(human_bins_counts)[dim(human_bins_counts)[2]]<-bam.oi
}
human_bins_counts<-human_bins_counts[,2:dim(human_bins_counts)[2]]
colnames(human_bins_counts) <- paste0(c(rep('H3K18la',2),rep('H3K27ac',2),rep('H3K27me3',2),rep('H3K4me3',2),rep('H3K9me3',2)),rep(c('_1','_2'),5))

# normalize to log2CPM 
libSize <- Matrix::colSums(human_bins_counts)
rpm <- sweep(human_bins_counts * 1e6, MARGIN=2, STATS=libSize, FUN="/")
human_bins_norm.counts <- log2(rpm+1)

# define metadata for plots
samples <- colnames(human_bins_norm.counts)
histone<-unlist(strsplit(samples,"_"))[c(1:10)*2-1]
rep<-unlist(strsplit(samples,"_"))[c(1:10)*2]

# plot MDS
pdf('/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/SupplFig6A_human_binned_genome.pdf',width=4,height=4)
plotMDS(human_bins_norm.counts, col=c("blue","orange","darkred","darkgreen","grey")[as.numeric(as.factor(histone))], cex = 1,pch=20)
dev.off()


# Fig 4A : Adhi
# Fig 4B/SFig 6B:

# Load all human peak sets
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks_human/peaks")
peak.list <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks_human/peaks",pattern=".bed")
peaks <-list()
for (i in 1:length(peak.list)){
peaks[[i]] <- read.delim(peak.list[i])
colnames(peaks[[i]])<-c('peak.chr','peak.start','peak.end','peak.total.value','peak.max.value','peak.max.value.region')}
names(peaks) <- gsub("_0.01.stringent.bed","",peak.list)


# Annotate all human peak sets

promoter_hs <- getPromoters(TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, upstream=3000, downstream=3000)
promoter_mm <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.hg38.knownGene, upstream=3000, downstream=3000)

peaks.ggranges<-list()
for(i in 1:length(peaks)){
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

peak_anno.human<-list()
for (i in 1:length(peaks.ggranges)){
  peak_anno.human[[i]] <- annotatePeak(peaks.ggranges[[i]], 
                                  tssRegion=c(-3000, 3000), 
                                  TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                  annoDb='org.Hs.eg.db')
}
names(peak_anno.human) <- names(peaks)

# Figure 4B: Adhi
pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/Fig4B.peaks.genomic.annotation.human.pdf")
plotAnnoBar(peak_anno.human[grep('merged',names(peak_anno.human))])
plotDistToTSS(peak_anno.human[grep('merged',names(peak_anno.human))],
              title = 'Distribution relative to TSS')
dev.off()

pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/SupplFig6B-C.peaks.genomic.annotation.human.pdf")
plotAnnoBar(peak_anno.human[-grep('merged',names(peak_anno.human))])
plotDistToTSS(peak_anno.human[-grep('merged',names(peak_anno.human))],
                title = 'Distribution relative to TSS')
dev.off()

# Make 1 master peak file containing all samples

peaks<-bind_rows(peaks, .id = "column_label")
peaks <- peaks[,c(2,3,4,5,6,7,1)]
write.table(data.frame(peaks),'masterpeak.human.bed',sep='\t',quote=F, col.names=F,row.names=F)


# Fig 4C: Adhi











