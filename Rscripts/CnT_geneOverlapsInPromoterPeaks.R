# This R script has been used to generate the following figures:
### Figures 1F, 4D
### Supp. Figures 3C

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomicRanges)
library(ggVennDiagram)
library(gridExtra)							
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

# Select peaks for specific histone mark (H3K27ac for example)  
subset_peak_anno.mouse <- peakAnno[grepl('GAS_H3K27ac|PIM_H3K27ac_1|ESC.ser_H3K27ac',names(peak_anno.mouse))]
subset_peaks <- samples[grepl('GAS_H3K27ac|PIM_H3K27ac_1|ESC.ser_H3K27ac',names(peak_anno.mouse))]

# Select peaks in promoters
subset_peaks.promoters <-list()
subset_peak_anno.mouse.promoters <- list()
subset_peaks.promoters.ggranges <- list()

for (i in 1:length(subset_peak_anno.mouse)){
  subset_peaks.promoters[[i]] <- subset_peaks[[i]][which(subset_peak_anno.mouse[[i]]@detailGenomicAnnotation$Promoter),]
  subset_peaks.promoters.ggranges[[i]] <- makeGRangesFromDataFrame(subset_peaks.promoters[[i]],
                                                                  keep.extra.columns=TRUE,
                                                                  ignore.strand=T,
                                                                  seqinfo=NULL,
                                                                  seqnames.field=c("chromosome", "chrom",
                                                                                   "chr", "chromosome_name"),
                                                                  start.field="start",
                                                                  end.field=c("end", "stop"),
                                                                  strand.field="strand",
                                                                  starts.in.df.are.0based=T) 
  subset_peak_anno.mouse.promoters[[i]] <- annotatePeak(subset_peaks.promoters.ggranges[[i]],
                                                        tssRegion=c(-2000, 2000),
                                                        TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                        annoDb='org.Mm.eg.db')
}

names(subset_peak_anno.mouse.promoters) <- names(subset_peaks)

# Plot Venn diagram showing gene overlaps based on peaks falling into promoter region
genes.w.promoter.peaks <- lapply(subset_peak_anno.mouse.promoters, function(i) as.data.frame(i)$geneId)
genes.w.promoter.peaks <- genes.w.promoter.peaks[order(names(genes.w.promoter.peaks))]
ggvennlist_promoters <- list()

ggvennlist_promoters[[1]]<-ggVennDiagram(genes.w.promoter.peaks[grep('GAS',names(genes.w.promoter.peaks))],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkred","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('GAS')+ scale_color_manual(values=c("blue","orange","darkred","darkgreen"))
ggvennlist_promoters[[2]]<-ggVennDiagram(genes.w.promoter.peaks[grep('PIM',names(genes.w.promoter.peaks))],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkred","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('PIM')+ scale_color_manual(values=c("blue","orange","darkred","darkgreen"))
ggvennlist_promoters[[3]]<-ggVennDiagram(genes.w.promoter.peaks[grep('ESC.ser',names(genes.w.promoter.peaks))],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkred","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('PIM')+ scale_color_manual(values=c("blue","orange","darkred","darkgreen"))

ggvennlist_promoters[[4]]<-ggVennDiagram(genes.w.promoter.peaks[grep('GAS',names(genes.w.promoter.peaks))][c(1,2,4)],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('GAS')+ scale_color_manual(values=c("blue","orange","darkgreen"))
ggvennlist_promoters[[5]]<-ggVennDiagram(genes.w.promoter.peaks[grep('PIM',names(genes.w.promoter.peaks))][c(1,2,4)],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('PIM')+ scale_color_manual(values=c("blue","orange","darkgreen"))
ggvennlist_promoters[[6]]<-ggVennDiagram(genes.w.promoter.peaks[grep('ESC.ser',names(genes.w.promoter.peaks))][c(1,2,4)],
                                         label_alpha = 0,set_size = 3,label_size = 3,set_color = c("blue","orange","darkgreen"))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('ESC-ser')+ scale_color_manual(values=c("blue","orange","darkgreen"))

do.call(grid.arrange,c(ggvennlist_promoters[1:6],nrow=2))
