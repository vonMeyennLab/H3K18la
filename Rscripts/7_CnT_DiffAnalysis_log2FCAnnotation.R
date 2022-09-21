# This R script has been used to generate the following figures:
### Figure 3A
### Supp. Figure 5A

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(ggpubr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(vioplot)

# Load raw peak counts based on the union of MB and MT peak sets (master peak)
MB_MT_H3K18la_rawPeakCounts<-readRDS("../MB_MT_H3K18la_rawPeakCounts.rds")

# Normalize the counts
libSize <- Matrix::colSums(MB_MT_H3K18la_rawPeakCounts)
MB_MT_H3K18la_rpmCounts <- sweep(MB_MT_H3K18la_rawPeakCounts * 1e6, MARGIN=2, STATS=libSize, FUN="/")

# Load raw peak counts based on the union of ESC-2i and ESC-ser peak sets (master peak)
ESC_H3K18la_rawPeakCounts<-readRDS("../ESC_H3K18la_rawPeakCounts.rds")

# Normalize the counts
libSize <- Matrix::colSums(ESC_H3K18la_rawPeakCounts)
ESC_H3K18la_rpmCounts <- sweep(ESC_H3K18la_rawPeakCounts * 1e6, MARGIN=2, STATS=libSize, FUN="/")

# Calculate fold-change from MB to MT and from 2i to ser
MB_MT_H3K18la_rpmFC <- (1+MB_MT_H3K18la_rpmCounts)/(MB_MT_H3K18la_rpmCounts[,1]+1)
MB_MT_H3K18la_rpmFC.mean<-rowMeans(MB_MT_H3K18la_rpmFC[,2:3])

ESC_H3K18la_rpmFC <- ESC_H3K18la_rpmCounts/(1+rowMeans(ESC_H3K18la_rpmCounts[,3:4]))
ESC_H3K18la_rpmFC.mean<-rowMeans(ESC_H3K18la_rpmFC[,1:2])

# Extract master peak coordinates
temp <- unlist(strsplit(split = '-',unlist(strsplit(split = ':',row.names(MB_MT_H3K18la_rpmCounts)))))
MB_MT_H3K18la_rpmCounts_coordinates <- data.frame(paste0('chr',c(temp[c(1:length(MB_MT_H3K18la_rpmCounts[,1]))*3-2])),
                                                  c(temp[c(1:length(MB_MT_H3K18la_rpmCounts[,1]))*3-1]),
                                                  c(temp[c(1:length(MB_MT_H3K18la_rpmCounts[,1]))*3]))
colnames(MB_MT_H3K18la_rpmCounts_coordinates) <- c('chr','start','end')

temp <- unlist(strsplit(split = '-',unlist(strsplit(split = ':',row.names(ESC_H3K18la_rpmCounts)))))
ESC_H3K18la_rpmCounts_coordinates <- data.frame(paste0('chr',c(temp[c(1:length(ESC_H3K18la_rpmCounts[,1]))*3-2])),
                                                c(temp[c(1:length(ESC_H3K18la_rpmCounts[,1]))*3-1]),
                                                c(temp[c(1:length(ESC_H3K18la_rpmCounts[,1]))*3]))
colnames(ESC_H3K18la_rpmCounts_coordinates) <- c('chr','start','end')

# Convert master peak sets to GRanges objects
MB_MT_peaks.ggranges<- makeGRangesFromDataFrame(MB_MT_H3K18la_rpmCounts_coordinates,
                                                keep.extra.columns=TRUE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c("chromosome", "chrom",
                                                                 "chr", "chromosome_name"),
                                                start.field="start",
                                                end.field=c("end", "stop"),
                                                strand.field="strand",
                                                starts.in.df.are.0based=FALSE)
ESC_peaks.ggranges<- makeGRangesFromDataFrame(ESC_H3K18la_rpmCounts_coordinates,
                                              keep.extra.columns=TRUE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field=c("chromosome", "chrom",
                                                               "chr", "chromosome_name"),
                                              start.field="start",
                                              end.field=c("end", "stop"),
                                              strand.field="strand",
                                              starts.in.df.are.0based=FALSE)

# Annotate master peaks to genomic features
promoter_mm <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=2000, downstream=2000)

peaks_anno_MB_MT <- annotatePeak(MB_MT_peaks.ggranges,
                                 tssRegion=c(-2000, 2000),
                                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                 annoDb='org.Mm.eg.db')
peaks_anno_ESC <- annotatePeak(ESC_peaks.ggranges,
                               tssRegion=c(-2000, 2000),
                               TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                               annoDb='org.Mm.eg.db')

# Define colors
color <- peaks_anno_MB_MT@anno@elementMetadata$annotation
color[grep('Intron',color)] <- 'Intron'
color[grep('Exon',color)] <- 'Exon'
color_MB_MT <- color

color <- peaks_anno_ESC@anno@elementMetadata$annotation
color[grep('Intron',color)] <- 'Intron'
color[grep('Exon',color)] <- 'Exon'
color_ESC <- color

# Box plots to visualize log2FC distribution across different genomic fetaures
par(mfrow=c(2,2),mar=c(6,10,5,5))
boxplot(abs(log2(rowMeans(MB_MT_H3K18la_rpmCounts[,2:3])/(1+MB_MT_H3K18la_rpmCounts[,1])))~color_MB_MT,
        horizontal=T,las=1, main = 'MB to MT',ylab='',xlab='abs(log2 FC)',notch=T,outpch=16,outcex=0.4)
abline(v=median(abs(log2(rowMeans(MB_MT_H3K18la_rpmCounts[,2:3])/(1+MB_MT_H3K18la_rpmCounts[,1])))[color_MB_MT=='Promoter (<=1kb)']),lty=2,col='red')

par(mfrow=c(2,2),mar=c(6,10,5,5))
boxplot(abs(log2(rowMeans(ESC_H3K18la_rpmCounts[,1:2])/(1+rowMeans(ESC_H3K18la_rpmCounts[,3:4]))))~color_ESC,
        horizontal=T,las=1, main = '2i to ser',ylab='',xlab='abs(log2 FC)',notch=T,outpch=16,outcex=0.4)
abline(v=median(abs(log2(rowMeans(ESC_H3K18la_rpmCounts[,3:4])/(1+ESC_H3K18la_rpmCounts[,1:2])))[color_ESC=='Promoter (<=1kb)']),lty=2,col='red')
