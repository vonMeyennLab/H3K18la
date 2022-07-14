library(chromVAR)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library('org.Mm.eg.db')
library('TxDb.Mmusculus.UCSC.mm10.knownGene')
library(ChIPseeker)
library(ggVennDiagram)
library(ggplot2)
library(GenomicAlignments)
library(DESeq2)
library(edgeR)
library(ggcorrplot)
library(ggfortify)
library(clusterProfiler)

# Fig 1A: WB
# Fig S1A: Lactate
# Fig S1B: WB
# Fig S1C: Workflow

# Fig 1B+S1D:

########################################################################################
###################### Quantify bams inside bins ######################################
########################################################################################

# load binned genome
binned_genome <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/mm10_genes.bed",header=F)
colnames(binned_genome) <- c('Chromosome','Start','End')
binned_genome<-binned_genome[-grep('CHR',binned_genome$Chromosome),]
binned_genome<-binned_genome[-grep('GL',binned_genome$Chromosome),]
binned_genome<-binned_genome[-grep('JH',binned_genome$Chromosome),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$Chromosome,
                          IRanges(start = binned_genome$Start-2000,
                                  end = binned_genome$Start+2000))
binned_genome<-binned_genome[-grep('MT',binned_genome$Chromosome),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$Chromosome,
                         IRanges(start = binned_genome$Start,
                                 end = binned_genome$End))

# read metadata with bamPaths for all files
metadata <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/metadata_Eva.txt")

# count reads within bins for all bams
mouse_bins_counts <- matrix(data = NA, nrow=length(binned_chromosome_ggranges),ncol=1)
for (i in c(1:length(metadata$bamPath))){
bam.oi<-metadata$bamPath[[i]]
fragment_counts_mouse <- summarizeOverlaps(binned_chromosome_ggranges,bam.oi,inter.feature=F,ignore.strand=T) %>% assays() %>% .[[1]]
mouse_bins_counts<- cbind(mouse_bins_counts,fragment_counts_mouse)
colnames(mouse_bins_counts)[dim(mouse_bins_counts)[2]]<-bam.oi
}
mouse_bins_counts<-mouse_bins_counts[,2:dim(mouse_bins_counts)[2]]
colnames(mouse_bins_counts) <- metadata$SampleID
all_samples_names<-colnames(mouse_bins_counts) 

# normalise raw counts for MDS
libSize <- Matrix::colSums(mouse_bins_counts)
#write.table(libSize, "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/libSize_all_samples.txt", sep='\t',col.names=F,row.names=T,quote=F)
rpm <- sweep(mouse_bins_counts * 1e6, MARGIN=2, STATS=libSize, FUN="/")
normCounts <- log2(rpm+1)

histone<-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))]
condition <-unlist(strsplit(all_samples_names,"_"))[grep('H3',unlist(strsplit(all_samples_names,"_")))-1]

#plot MDS
pdf('/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/SFig1D.mouse_mds_binned_genome.pdf',width=5,height=4)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
plotMDS(normCounts, col=c("blue","orange","darkred","darkgreen")[as.numeric(as.factor(histone))],  cex = 1, pch = c(1,2,3,4,5,6,7,8)[as.numeric(as.factor(condition))])
legend(par("usr")[2], par("usr")[4], legend=levels(as.factor(condition)), pch=c(1,2,3,4,5,6,7,8), col="black", ncol=1, cex=0.6)
dev.off()
        
# Repeat for active marks only
pdf('/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig1B.mouse_mds_binned_genome.pdf',width=5,height=4)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
normCounts_active <- normCounts[,-grep('H3K27me3', colnames(normCounts))]
plotMDS(normCounts_active, col=c("blue","orange","darkgreen")[as.numeric(as.factor(histone[-grep('H3K27me3', colnames(normCounts))]))],  cex = 1, pch = c(1,2,3,4,5,6,7,8)[as.numeric(as.factor(condition[-grep('H3K27me3', colnames(normCounts))]))])
legend(par("usr")[2], par("usr")[4], legend=levels(as.factor(condition)), pch=c(1,2,3,4,5,6,7,8), col="black", ncol=1, cex=0.6)
dev.off()

# Fig 1C: Adhi
# SFig 1C: IGV screenshots

# Fig 1D, SFig1E-F

# Read all peaks used in this project
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks")
peak.list <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks",pattern="sorted.bed")
peaks<-list()
for (i in 1:length(peak.list)){
peaks[[i]] <- read.delim(peak.list[i])
colnames(peaks[[i]])<-c('peak.chr','peak.start','peak.end','peak.total.value','peak.max.value','peak.max.value.region')}
names(peaks) <- gsub("_0.01.stringent.bed.w.chr.sorted.bed","",peak.list[1:length(peak.list)])
names(peaks) <- gsub(".0.01.stringent.bed.w.chr.sorted.bed","",names(peaks))
names(peaks) <- gsub("_0.01.stringent.bed.w.chr.bed.sorted.bed","",names(peaks))

# For table S1 : Adhi

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
  peakN       <- data.frame(peakN = nrow(peakInfo), peakTotalBases = sum(peakInfo$width), histone = sampleInfo[2], biology = sampleInfo[1] , rep = sampleInfo[3]) %>% rbind(peakN, .)
}

# Annotate peaks to genomic regions
# Get TSS regions
promoter_mm <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=2000, downstream=2000)

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
subset_peak_anno.mouse <- peak_anno.mouse[grep('H3K18la',names(peak_anno.mouse))]
subset_peak_anno.mouse<-subset_peak_anno.mouse[c(3,4,7,10,13,14,17,20)]
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig1C.H3K18la.merged.genomic.feature.distribution.pdf", height=2, width=7)
plotAnnoBar(subset_peak_anno.mouse)
dev.off()

subset_peak_anno.mouse <- peak_anno.mouse[-grep('merged',names(peak_anno.mouse))]
PTM.order <- unlist(strsplit(names(subset_peak_anno.mouse),'_'))
PTM.order <- PTM.order[grep('H3',PTM.order)]
subset_peak_anno.mouse<-subset_peak_anno.mouse[order(PTM.order)]

pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/SFig1F-G.peaks.genomic.annotation.pdf", height=5.5, width=5)
plotAnnoBar(subset_peak_anno.mouse)
plotDistToTSS(subset_peak_anno.mouse,
              title = 'Distribution relative to TSS')
dev.off()

# Fig 1E: Adhi

# Fig 1F:

# Select peaks in promoters and annotate to genes
subset_peak_anno.mouse <- peak_anno.mouse[grepl('merged|PIM_H3K27ac_1|ESC.ser_H3K27ac',names(peak_anno.mouse))]
subset_peaks <- peaks[grepl('merged|PIM_H3K27ac_1|ESC.ser_H3K27ac',names(peak_anno.mouse))]
subset_peaks.promoters <-list()
subset_peak_anno.mouse.promoters <- list()
subset_peaks.promoters.ggranges <- list()

for (i in 1:length(subset_peak_anno.mouse)){
	subset_peaks.promoters[[i]]<-subset_peaks[[i]][which(subset_peak_anno.mouse[[i]]@detailGenomicAnnotation$Promoter),]
	subset_peaks.promoters.ggranges[[i]]<- makeGRangesFromDataFrame(subset_peaks.promoters[[i]],
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

genes.w.promoter.peaks <- lapply(subset_peak_anno.mouse.promoters, function(i) as.data.frame(i)$geneId)
genes.w.promoter.peaks <- genes.w.promoter.peaks[order(names(genes.w.promoter.peaks))]
ggvennlist_promoters <- list()

library(ggVennDiagram)

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


library(gridExtra)							
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig1F.andSFig3C.promoters.covered.by.peaks.venn.pdf",width=10,height=10)
do.call(grid.arrange,c(ggvennlist_promoters[1:6],nrow=2))
dev.off()

# Fig 1G: Adhi

# Fig 1H:

######## Correlation to gene expression: 

genes.w.promoter.peaks <- lapply(subset_peak_anno.mouse.promoters, function(i) as.data.frame(i)$SYMBOL)
genes.w.promoter.peaks <- genes.w.promoter.peaks[order(names(genes.w.promoter.peaks))]
names(genes.w.promoter.peaks) <- gsub(".w.chr.sorted.bed",'',names(genes.w.promoter.peaks))
library(gplots)

venn_GAS<- venn(genes.w.promoter.peaks[grep('GAS',names(genes.w.promoter.peaks))][c(1,2,4)], show.plot = T)
venn_PIM<- venn(genes.w.promoter.peaks[grep('PIM',names(genes.w.promoter.peaks))][c(1,2,4)], show.plot = T)
venn_ESC.2i<- venn(genes.w.promoter.peaks[grep('ESC.2i',names(genes.w.promoter.peaks))], show.plot = T)
venn_ESC.ser<- venn(genes.w.promoter.peaks[grep('ESC.ser',names(genes.w.promoter.peaks))][c(1,2,4)], show.plot = T)

names(attributes(venn_GAS)$intersection)<-gsub('_merged','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_GAS)$intersection)<-gsub('GAS_','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_merged','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_1','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('PIM_','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_1','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_public.bed','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_merged','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('ESC.ser_','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.2i)$intersection)<-gsub('_merged','',names(attributes(venn_ESC.2i)$intersection))
names(attributes(venn_ESC.2i)$intersection)<-gsub('ESC.2i_','',names(attributes(venn_ESC.2i)$intersection))

#### Correlate to gene expression
suppressPackageStartupMessages({
  library(biomaRt)
  library(DESeq2)
  library(dplyr)
  library(edgeR)
  library(ggfortify)
  library(ggplot2)
  library(Matrix)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(Rsubread)
  library(scales)
  library(tibble)
  library(tidyverse)
  library(WGCNA)
  library(scales)
  library(beeswarm)
  
})
suppressPackageStartupMessages({
  library(biomaRt)
  library(DESeq2)
  library(dplyr)
  library(edgeR)
  library(ggfortify)
  library(ggplot2)
  library(Matrix)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(Rsubread)
  library(scales)
  library(tibble)
  library(tidyverse)
  library(WGCNA)
})


# Add RNA seq data and annotate to gene names 
GAS <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_promoter/GAS_hPTM_RNA_Pr_normCounts.txt")
PIM <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_promoter/PIM_hPTM_RNA_Pr_normCounts.txt")
ESC.ser <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_RNA_promoter/ESC_hPTM_RNA_Pr_normCounts.txt")

# setwd("/home/egalle/public/EvaGalle/Data/RNA/Sequencing/Counts/")
# rawCounts_GAS <- readRDS("merged_rawCounts.SEQ00030.rds")  
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")
geneAnnotation <- getBM(attributes=c("ensembl_gene_id",'entrezgene_id', "external_gene_name","start_position", "end_position"), 
                        filters="ensembl_gene_id", 
                        values=PIM$EnsemblID, mart=ensembl)
geneAnnotation$width = geneAnnotation$end_position-geneAnnotation$start_position+1   
                      
# # Calculate RPKM values
# libSize <- Matrix::colSums(rawCounts_GAS)
# rpm <- sweep(rawCounts_GAS * 1e6, MARGIN=2, STATS=libSize, FUN="/")
rpm <- merge(PIM, geneAnnotation,
                    by.y="ensembl_gene_id", 
                    by.x="EnsemblID",all.y=T,all.x=T)
rpm$RNA<-(2^rpm$RNA)/rpm$width*1000
x <- list(unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la'])))
names(x) <- c('all','H3K27ac-H3K4me3','H3K4me3','H3K18la-H3K27ac','H3K18la-H3K4me3','H3K27ac','H3K18la')

saveRDS(x, "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/promoters/hPTM_promoter_occupancy_PIM.RDS")

rpm$Promoter <- NA
for (i in 1:length(names(x))){
	rpm$Promoter[rpm$GeneName%in%x[[i]]]<- names(x)[i]
}
rpm$Promoter[is.na(rpm$Promoter)] <-'no active mark'
rpm_PIM <- rpm


# Repeat for GAS
geneAnnotation <- getBM(attributes=c("ensembl_gene_id",'entrezgene_id', "external_gene_name","chromosome_name","start_position", "end_position"), 
                        filters="ensembl_gene_id", 
                        values=GAS$EnsemblID, mart=ensembl)
geneAnnotation$width = geneAnnotation$end_position-geneAnnotation$start_position   
                      
# # Calculate RPKM values
rpm <- merge(GAS, geneAnnotation,
                    by.y="ensembl_gene_id", 
                    by.x="EnsemblID",all.y=T,all.x=T)
rpm$RNA<-(2^rpm$RNA)/rpm$width*1000

x <- list(unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la:H3K27ac'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la:H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K27ac'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la'])))
names(x) <- c('all','H3K27ac-H3K4me3','H3K4me3','H3K18la-H3K27ac','H3K18la-H3K4me3','H3K27ac','H3K18la')
saveRDS(x, "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/promoters/hPTM_promoter_occupancy_GAS.RDS")

rpm$Promoter <- NA
for (i in 1:length(names(x))){
	rpm$Promoter[rpm$GeneName%in%x[[i]]]<- names(x)[i]
}
rpm$Promoter[is.na(rpm$Promoter)] <-'no active mark'
rpm_GAS <- rpm

# Repeat for ESC.ser
geneAnnotation <- getBM(attributes=c("ensembl_gene_id",'entrezgene_id', "external_gene_name","chromosome_name","start_position", "end_position"), 
                        filters="ensembl_gene_id", 
                        values=ESC.ser$EnsemblID, mart=ensembl)
geneAnnotation$width = geneAnnotation$end_position-geneAnnotation$start_position   
                      
# # Calculate RPKM values
rpm <- merge(ESC.ser, geneAnnotation,
                    by.y="ensembl_gene_id", 
                    by.x="EnsemblID",all.y=T,all.x=T)
rpm$RNA<-(2^rpm$RNA)/rpm$width*1000

x <- list(unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la'])))
names(x) <- c('all','H3K27ac-H3K4me3','H3K4me3','H3K18la-H3K27ac','H3K18la-H3K4me3','H3K27ac','H3K18la')
saveRDS(x, "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/promoters/hPTM_promoter_occupancy_mESC.ser.RDS")

rpm$Promoter <- NA
for (i in 1:length(names(x))){
	rpm$Promoter[rpm$GeneName%in%x[[i]]]<- names(x)[i]
}
rpm$Promoter[is.na(rpm$Promoter)] <-'no active mark'
rpm_ESC.ser<-rpm

colors_to_use <- cbind(c('red','blue','lightblue','darkblue','orange','orange3','darkgreen','red'), sort(unique(rpm_GAS$Promoter)))
colnames(colors_to_use) <- c('color','Promoter')
	 
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig1H.promoters.covered.by.active.peaks.RNA.boxplots.pdf",width=4,height=8)
par(mar=c(5,10,5,5),mfrow=c(3,1))
# Reorder groups to get them from high to low
new_order <- with(rpm_GAS,reorder(Promoter, RNA, median, na.rm=T))
no2 <- data.frame(rpm_GAS %>% group_by(Promoter) %>% summarise(median=median(RNA)))
no2 <- merge(no2, colors_to_use, by = "Promoter")
boxplot(log2(rpm_GAS$RNA+1)~new_order,outline=F,
	col=no2$color[order(no2$median)], ylab = '', xlab = 'Gene expression log2(RPKM)',notch=T,
	horizontal=T,las=1,main='GAS',ylim=c(0,7))
text(y=1:8, x=6, paste0('n=',table(new_order)),cex=0.7)
t.test(log2(rpm_GAS$RNA[rpm$Promoter=='all']+1),log2(rpm_GAS$RNA[rpm$Promoter=='H3K27ac-H3K4me3']+1))	


new_order <- with(rpm_PIM,reorder(Promoter, RNA, median, na.rm=T))
no2 <- data.frame(rpm_PIM %>% group_by(Promoter) %>% summarise(median=median(RNA)))
no2 <- merge(no2, colors_to_use, by = "Promoter")
boxplot(log2(rpm_PIM$RNA+1)~new_order,outline=F,
	col=no2$color[order(no2$median)], ylab = '', xlab = 'Gene expression log2(RPKM)',notch=T,
	horizontal=T,las=1,main='PIM',ylim=c(0,7))

text(y=1:8, x=6, paste0('n=',table(new_order)),cex=0.7)

new_order <- with(rpm_ESC.ser,reorder(Promoter, RNA, median, na.rm=T))
no2 <- data.frame(rpm_ESC.ser %>% group_by(Promoter) %>% summarise(median=median(RNA)))
no2 <- merge(no2, colors_to_use, by = "Promoter")
boxplot(log2(rpm_ESC.ser$RNA+1)~new_order,outline=F,
	col=no2$color[order(no2$median)], ylab = '', xlab = 'Gene expression log2(RPKM)',notch=T,
	horizontal=T,las=1,main='ESC-ser',ylim=c(0,7))

text(y=1:8, x=6, paste0('n=',table(new_order)),cex=0.7)

dev.off()


# Fig 1I: Adhi

# Supplemental Fig 2A,B,C : Adhi

# SFig3A: 

library(gplots)
names(genes.w.promoter.peaks) <- gsub('.w.chr.sorted.bed','', names(genes.w.promoter.peaks) )
venn_GAS<- venn(genes.w.promoter.peaks[grep('GAS',names(genes.w.promoter.peaks))], show.plot = T)
venn_PIM<- venn(genes.w.promoter.peaks[grep('PIM',names(genes.w.promoter.peaks))], show.plot = T)
venn_ESC.2i<- venn(genes.w.promoter.peaks[grep('ESC.2i',names(genes.w.promoter.peaks))], show.plot = T)
venn_ESC.ser<- venn(genes.w.promoter.peaks[grep('ESC.ser',names(genes.w.promoter.peaks))], show.plot = T)

names(attributes(venn_GAS)$intersection)<-gsub('_merged','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_GAS)$intersection)<-gsub('GAS_','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_merged','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_1','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('PIM_','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_merged','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_public.bed','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_1','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('ESC.ser_','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.2i)$intersection)<-gsub('_merged','',names(attributes(venn_ESC.2i)$intersection))
names(attributes(venn_ESC.2i)$intersection)<-gsub('ESC_','',names(attributes(venn_ESC.2i)$intersection))


x <- list(unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')
compMF_GAS_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="MF")

compBP_GAS_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="BP")

compCC_GAS_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="CC")


GO_GAS_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="all")         
                         

x <- list(unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')

GO_PIM_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="all")      
                               
compMF_PIM_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="MF")

compBP_PIM_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="BP")

compCC_PIM_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="CC")
                         
x <- list(unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')

GO_ESC.ser_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="all")      

                         
compMF_ESC.ser_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="MF")

compBP_ESC.ser_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="BP")

compCC_ESC.ser_sub <- compareCluster(geneCluster = x,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb='org.Mm.eg.db',
                         ont="CC")
                         

pdf("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/plots/promoters.covered.by.peaks.venn.GO.enrichment_ALL_sub_3_groups.pdf",width=9,height=10)

dotplot(GO_PIM_sub, showCategory=5, includeAll=T)
dotplot(GO_GAS_sub, showCategory=5, includeAll=T)
dotplot(GO_ESC.ser_sub, showCategory=5, includeAll=T)
dev.off()  

dp<-list()
dp[[1]] <- dotplot(compBP_GAS_sub, showCategory=5, title = "GAS BP GO Enrichment Analysis", includeAll=F)
dp[[2]] <- dotplot(compMF_GAS_sub, showCategory=5, title = "GAS MF GO Enrichment Analysis", includeAll=F)
dp[[3]] <- dotplot(compCC_GAS_sub, showCategory=5, title = "GAS CC GO Enrichment Analysis", includeAll=F)

dp[[4]] <- dotplot(compBP_PIM_sub, showCategory=5, title = "PIM BP GO Enrichment Analysis", includeAll=F)
dp[[5]] <- dotplot(compMF_PIM_sub,showCategory=5, title = "PIM MF GO Enrichment Analysis", includeAll=F)
dp[[6]] <- dotplot(compCC_PIM_sub, showCategory=5, title = "PIM CC GO Enrichment Analysis", includeAll=F)

dp[[7]] <- dotplot(compBP_ESC.ser_sub, showCategory=5, title = "ESC.ser BP GO Enrichment Analysis",includeAll=F)
dp[[8]] <- dotplot(compMF_ESC.ser_sub, showCategory=5, title = "ESC.ser MF GO Enrichment Analysis", includeAll=F)
dp[[9]] <- dotplot(compCC_ESC.ser_sub, showCategory=5, title = "ESC.ser CC GO Enrichment Analysis", includeAll=F)

dp[[11]] <- dotplot(compBP_GAS_sub, showCategory=5, title = "GAS BP GO Enrichment Analysis", includeAll=T)
dp[[12]] <- dotplot(compMF_GAS_sub, showCategory=5, title = "GAS MF GO Enrichment Analysis", includeAll=T)
dp[[13]] <- dotplot(compCC_GAS_sub, showCategory=5, title = "GAS CC GO Enrichment Analysis", includeAll=T)

dp[[14]] <- dotplot(compBP_PIM_sub, showCategory=5, title = "PIM BP GO Enrichment Analysis", includeAll=T)
dp[[15]] <- dotplot(compMF_PIM_sub,showCategory=5, title = "PIM MF GO Enrichment Analysis", includeAll=T)
dp[[16]] <- dotplot(compCC_PIM_sub, showCategory=5, title = "PIM CC GO Enrichment Analysis", includeAll=T)

dp[[17]] <- dotplot(compBP_ESC.ser_sub, showCategory=5, title = "ESC.ser BP GO Enrichment Analysis", includeAll=T)
dp[[18]] <- dotplot(compMF_ESC.ser_sub, showCategory=5, title = "ESC.ser MF GO Enrichment Analysis", includeAll=T)
dp[[19]] <- dotplot(compCC_ESC.ser_sub, showCategory=5, title = "ESC.ser CC GO Enrichment Analysis", includeAll=T)

dp[[21]] <- dotplot(compCC_GAS_sub, showCategory=15, title = "ESC.ser CC GO Enrichment Analysis", includeAll=T)
dp[[22]] <- dotplot(compCC_PIM_sub, showCategory=15, title = "ESC.ser CC GO Enrichment Analysis", includeAll=T)
dp[[23]] <- dotplot(compCC_ESC.ser_sub, showCategory=15, title = "ESC.ser CC GO Enrichment Analysis", includeAll=T)

library(gridExtra)

pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Sfig3A.GO.different.combinations.of.promoter.hPTM.pdf",width=18,height=13)
do.call("grid.arrange",c(dp[c(23,21,22)],ncol=3))
dev.off()

# Sfig3B: Cistrome output

# Example how to generate a Cistrome input file: 

x <- list(unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K4me3'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac'])),
		unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la'])))
names(x) <- c('H3K18la-H3K27ac-H3K4me3','H3K27ac-H3K4me3','H3K4me3','H3K18la-H3K27ac','H3K18la-H3K4me3','H3K27ac','H3K18la')

db<-c()
for (i in 1:length(subset_peak_anno.mouse.promoters)){
	db <- rbind(db,as.data.frame(subset_peak_anno.mouse.promoters[[i]]))	
}

genes.w.promoter.coord <- list()
for (i in 1:length(x)){
	genes.w.promoter.coord[[i]] <- 	db[db$SYMBOL %in% x[[i]],]
	genes.w.promoter.coord[[i]]<-genes.w.promoter.coord[[i]][,colnames(genes.w.promoter.coord[[i]]) %in% c('geneChr','geneStart',"geneEnd","SYMBOL")]
	genes.w.promoter.coord[[i]]$geneChr<-paste0('chr',genes.w.promoter.coord[[i]]$geneChr)
	genes.w.promoter.coord[[i]]<-genes.w.promoter.coord[[i]][!duplicated(genes.w.promoter.coord[[i]]$SYMBOL),]
	genes.w.promoter.coord[[i]]<-genes.w.promoter.coord[[i]][,1:3]
	genes.w.promoter.coord[[i]][,3]<-genes.w.promoter.coord[[i]][,2]+2000
	genes.w.promoter.coord[[i]][,2]<-genes.w.promoter.coord[[i]][,2]-2000
	colnames(genes.w.promoter.coord[[i]]) <- c('chr','promStart','promEnd')
	
	}
names(genes.w.promoter.coord) <- names(x) 

ptm.oi<-'H3K18la-H3K27ac-H3K4me3'
write.table(genes.w.promoter.coord[[ptm.oi]], paste0("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/promoters/promoters_w_",ptm.oi,"_PIM.bed"),quote=F,col.names=F,row.names=F,sep="\t")


# Sfig3C: see above







