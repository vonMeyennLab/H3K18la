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

########################################################################################
###################### Counting bams inside bins ######################################
########################################################################################

# load the binned genome and format 
binned_genome <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/binned_genome/binned_genome_hg38.bed.txt")
binned_genome <- data.frame(binned_genome$Chromosome,binned_genome$Start,binned_genome$End)
colnames(binned_genome)<-c('chr','start','end')
binned_genome$chr<-paste0('chr',binned_genome$chr)
binned_genome<-binned_genome[-grep('MT',binned_genome$chr),]
binned_chromosome_ggranges   <- GRanges(seqnames = binned_genome$chr, 
                         IRanges(start = binned_genome$start, 
                                 end = binned_genome$end))	
								 
# load the metadata for all human samples containing paths to bam files in column 'bamPath'
metadata <- read.delim("~/metadata.txt")
	
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
plotMDS(human_bins_norm.counts, col=c("blue","orange","darkred","darkgreen","grey")[as.numeric(as.factor(histone))], cex = 1,pch=20)


########################################################################################
###############################  From here on; Peaks ######################################
########################################################################################

# Load all human peak sets
peak.list <- list.files(path="~/peaks",pattern=".bed")
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

# Figure 3A/C/D/H/I: Adhi
plotAnnoBar(peak_anno.human[grep('merged',names(peak_anno.human))])
plotDistToTSS(peak_anno.human[grep('merged',names(peak_anno.human))],
              title = 'Distribution relative to TSS')

plotAnnoBar(peak_anno.human[-grep('merged',names(peak_anno.human))])
plotDistToTSS(peak_anno.human[-grep('merged',names(peak_anno.human))],
                title = 'Distribution relative to TSS')

# Make 1 master peak file containing all samples

peaks<-bind_rows(peaks, .id = "column_label")
peaks <- peaks[,c(2,3,4,5,6,7,1)]
write.table(data.frame(peaks),'masterpeak.human.bed',sep='\t',quote=F, col.names=F,row.names=F)




################################################################################################
################# Create files containing peak overlap with cCREs ##############################
################################################################################################

# in bash: define and sort inputs
refbed=/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_hs/GRCh38-ccREs.bed
sort -k1,1 -k2,2n $refbed > $refbed.sorted.bed

peakPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks_human/peaks
for i in *0.01.stringent.*
do
sort -k1,1 -k2,2n $peakPath/${i} > $resultsPath/peaks_human/${i}.sorted.bed
done

# in bash: use bedtools intersect function to calculate overlap between each peak set and the cCREs
resultsPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE/overlap.cCRE/human/
cd $peakPath
for i in *0.01.stringent.*
do
cd $resultsPath
bedtools intersect -wao -a $resultsPath/peaks_human/${i}.sorted.bed -b $refbed.sorted.bed -sorted > $i.overlap.human.cCREs.bed
cd $peakPath
done

# Analyse overlap in R
# read cCRE file
cCREs <- read.table("/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_hs/GRCh38-ccREs.bed.sorted.bed")
colnames(cCREs) <- c('chr.cCRE','start.cCRE','end.cCRE','id1.cCRE','id2.cCRE','type.cCRE')

# define overlap files (1 per peak set)
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE/overlap.cCRE/human/")
peaks.overlap.cCREs.files <- list.files(pattern='overlap.cCREs.bed.2.bed') # the '.2.bed' extension indicates for me the files that were corrected for ^M containing lines. (corrected wih sed -e 's/^M//' command)

# read overlap files and store in list
peaks.overlap.cCRE <- list()
for (i in 1:length(peak.list)){
  peaks.overlap.cCRE[[i]] <- read.table(paste0("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE/overlap.cCRE/human/",peaks.overlap.cCREs.files[i]))
  colnames(peaks.overlap.cCRE[[i]]) <- c('chr.peak',
                                         'start.peak',
                                         'end.peak',
                                         'total.signal.peak',
                                         'max.signal.peak',
                                         'max.signal.peak.region',
                                         'chr.cCRE',
                                         'start.cCRE',
                                         'end.cCRE',
                                         'id1.cCRE',
                                         'id2.cCRE',
                                         'type.cCRE',
                                         'bp.overlap')
   colnames(peaks[[i]]) <-colnames(peaks.overlap.cCRE[[i]])[1:6]
}
names(peaks.overlap.cCRE)<-gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)

# for normalization: calculate the amount of bp all peaks in 1 peak set spans
total.peak.bp<-list()
for (i in 1:length(peaks)){
	peaks[[i]]$width <- peaks[[i]]$end.peak-peaks[[i]]$start.peak
	total.peak.bp[[i]]<-sum(peaks[[i]]$width)}
total.peak.bp<-unlist(total.peak.bp)
names(total.peak.bp)<-names(peaks)

# list peaks that overlap with specific cCREs and calculate the overlap (in # of bp) of each peak set with specific cCREs
peaks.w.cCRE.hg38 <- list()
peaks.bp.overlap.w.cCRE.hg38 <-list()
library(dplyr)
for (i in 1:length(peaks.overlap.cCRE)){
  peaks.oi<-peaks.overlap.cCRE[[i]]
  
  peaks.bp.overlap.w.cCRE.hg38[[i]]<- data.frame(peaks.oi %>% group_by(type.cCRE) %>% summarize(sum=sum(bp.overlap)))
  
  peaks.w.cCRE.hg38[[i]]<-list(unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('PLS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('pELS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('dELS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('CTCF-only,CTCF-bound',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('H3K4me3',peaks.oi$type.cCRE)]),
	   unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[peaks.oi$bp.overlap==0]))
                   names(peaks.w.cCRE.hg38[[i]]) <- c('PLS', 'pELS','dELS','CTCF-only,CTCF-bound','DNase-H3K4me3','no cCRE')
}
names(peaks.w.cCRE.hg38) <-gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)
names(peaks.bp.overlap.w.cCRE.hg38) <-gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)

cols<-cbind(names(peaks.w.cCRE.hg38),NA)
cols[,2][grep('H3K18la',cols[,1])]<-'blue'
cols[,2][grep('H3K14la',cols[,1])]<-'lightblue'
cols[,2][grep('H4',cols[,1])]<-'darkblue'
cols[,2][grep('H3K27me3',cols[,1])]<-'darkred'
cols[,2][grep('H3K9me3',cols[,1])]<-'grey'
cols[,2][grep('H3K4me3',cols[,1])]<-'darkgreen'
cols[,2][grep('H3K36me3',cols[,1])]<-'green'
cols[,2][grep('ac',cols[,1])]<-'orange'

ggvennlist<-list()
j=1
for (i in c(1:length(peaks.w.cCRE.hg38))){
  ggvennlist[[i]]<-ggVennDiagram(peaks.w.cCRE.hg38[[i]],label_alpha = 0,set_size = 3,label_size = 0,label=NULL)+
    ggplot2::scale_fill_gradient(low="white",high = cols[j,2])+ theme(text = element_text(size = 7),plot.title = element_text(size=10))+ggtitle(names(peaks.w.cCRE.hg38)[i])  
  j=j+1
}
# Figure 3A/C/D/H/I: Adhi
library(gridExtra)
pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/Fig3E.cCRE.human.venn.pdf",width=15,height=15)
do.call("grid.arrange",ggvennlist)
do.call("grid.arrange",ggvennlist[c(grep('merged',names(peaks.w.cCRE.hg38)))])
dev.off()


# Normalize overlap bp for size of peaks 

peaks.bp.overlap.w.cCRE.hg382 <- bind_cols(peaks.bp.overlap.w.cCRE.hg38)
peaks.bp.overlap.w.cCRE.hg382 <-peaks.bp.overlap.w.cCRE.hg382[,c(1:c(dim(peaks.bp.overlap.w.cCRE.hg382)[2]/2))*2]
colnames(peaks.bp.overlap.w.cCRE.hg382) <- gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)
rownames(peaks.bp.overlap.w.cCRE.hg382)<-peaks.bp.overlap.w.cCRE.hg38[[1]]$type.cCRE
peaks.bp.overlap.w.cCRE.hg382 <-rbind(peaks.bp.overlap.w.cCRE.hg382,total.peak.bp)
peaks.bp.overlap.w.cCRE.hg383 <-t(t(as.matrix(peaks.bp.overlap.w.cCRE.hg382))/as.numeric(peaks.bp.overlap.w.cCRE.hg382[11,]))

# Normalize overlap bp for size of cCRE

cCRE.per.type<- data.frame(cCREs %>% group_by(type.cCRE) %>% summarize(sum=sum(end.cCRE-start.cCRE)))
peaks.bp.overlap.w.cCRE.hg384 <-as.matrix(peaks.bp.overlap.w.cCRE.hg383)/c(1,cCRE.per.type$sum,1)
peaks.bp.overlap.w.cCRE.hg385 <- peaks.bp.overlap.w.cCRE.hg384[2:10,]
peaks.bp.overlap.w.cCRE.hg386 <- peaks.bp.overlap.w.cCRE.hg385[,grep('merged',colnames(peaks.bp.overlap.w.cCRE.hg385))]
colnames(peaks.bp.overlap.w.cCRE.hg386) <- gsub('Muscle_','',colnames(peaks.bp.overlap.w.cCRE.hg386))
colnames(peaks.bp.overlap.w.cCRE.hg386) <- gsub('_merged.2.bed','',colnames(peaks.bp.overlap.w.cCRE.hg386))

# create enrichment heatmap

breaksList = seq(-2, 2, by = 0.1)
library(pheatmap)
library(RColorBrewer)
pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/Fig3F.cCRE.human.heatmap.pdf",width=8,height=8)
pheatmap(peaks.bp.overlap.w.cCRE.hg386,cellheight=13, scale='column',main='column-normalized',breaks = breaksList,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)))
dev.off()


## Tissue specific enhancer overlap

# This is in bash. Sort input files to make intersect faster
resultsPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers
QUERYbed=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.human.bed
sort -k1,1 -k2,2n $QUERYbed > $QUERYbed.sorted.bed
refbed1=/home/egalle/public/EvaGalle/Data/Public_Data/Williams_2021_human.muscle/human_muscle_enhancers.bed
sort -k1,1 -k2,2n $refbed1 > $refbed1.sorted.bed

# This is in bash. Intersect enhancers with peaks
cd $resultsPath
bedtools intersect -wao -a ${QUERYbed}.sorted.bed -b $refbed1.sorted.bed -sorted > ${QUERYbed}.overlap_w_muscle_enhancers.bed


# From here back in R:
R

setwd( "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers")
muscle.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/masterpeak.human.bed.overlap_w_muscle_enhancers.bed",header=F)
colnames(muscle.enhancers.covered.by.peaks) <- c('peak.chr','peak.start','peak.end',
'peak.total.value','peak.max.value','peak.max.value.region','peak.origin',"chr.enhancer",'start.enhancer','end.enhancer','width.enhancer','logFC (pre- vs post-training)',
'logCPM	p-value	FDR','p-value','FDR','overlap.bp')

enhancers.all<-read.delim("/home/egalle/public/EvaGalle/Data/Public_Data/Williams_2021_human.muscle/human_muscle_enhancers.bed")

muscle.enhancers.covered <- c()
muscle.enhancers.covered.coord <- list()
for (i in 1:length(unique(muscle.enhancers.covered.by.peaks$peak.origin)))
{
	peakOI<-unique(muscle.enhancers.covered.by.peaks$peak.origin)[i]
	muscle.enhancers.covered[i]<-length(unique(paste0(muscle.enhancers.covered.by.peaks$chr.enhancer,muscle.enhancers.covered.by.peaks$start.enhancer,muscle.enhancers.covered.by.peaks$end.enhancer)[muscle.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(paste0(enhancers.all[,1],enhancers.all[,2],enhancers.all[,3])))
	muscle.enhancers.covered.coord[[i]]<-unique(cbind(muscle.enhancers.covered.by.peaks$chr.enhancer,muscle.enhancers.covered.by.peaks$start.enhancer,muscle.enhancers.covered.by.peaks$end.enhancer)[muscle.enhancers.covered.by.peaks$peak.origin==peakOI,])
	colnames(muscle.enhancers.covered.coord[[i]]) <- c('chr','start','end')

}
names(muscle.enhancers.covered) <- gsub("Muscle_","",unique(muscle.enhancers.covered.by.peaks$peak.origin))
names(muscle.enhancers.covered)[names(muscle.enhancers.covered)=='.']<-'not covered'
names(muscle.enhancers.covered.coord)<-names(muscle.enhancers.covered)
cols<-cbind(names(muscle.enhancers.covered),NA)
cols[,2][grep('.',cols[,1])] <-'grey'
cols[,2][grep('H3K18la',cols[,1])]<-'blue'
cols[,2][grep('H3K14la',cols[,1])]<-'lightblue'
cols[,2][grep('H4',cols[,1])]<-'darkblue'
cols[,2][grep('H3K27me3',cols[,1])]<-'darkred'
cols[,2][grep('H3K4me3',cols[,1])]<-'darkgreen'
cols[,2][grep('H3K36me3',cols[,1])]<-'green'
cols[,2][grep('ac',cols[,1])]<-'orange'
cols[,2][grep('H3K9me3',cols[,1])]<-'grey'


pdf('/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/Fig3G.Fraction of human muscle enhancers covered by peaks.pdf')
par(mar=c(15,20,20,5),cex=0.7)
barplot(muscle.enhancers.covered[order(muscle.enhancers.covered)][grep('merged',names(muscle.enhancers.covered[order(muscle.enhancers.covered)]))],horiz=T,
        las=1,xlim=c(0,1),col=cols[,2][order(muscle.enhancers.covered)][grep('merged',names(muscle.enhancers.covered[order(muscle.enhancers.covered)]))],
        xlab='Fraction of muscle enhancers covered',)
dev.off()

# Figure 3A/C/D/H/I: Adhi

## Correlation promoters to gene expression:
#Adhi, can you add the script to make these files here?

Prom_vs_expr <- read.delim("/home/egalle/public/_Projects/MS_H3K18la//Human/hPTM_RNA_promoter/hsMuscle_hPTM_RNA_Pr_normCounts.txt")
pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/SupplFig4D.cor.hPTM.promoter.vs.RNA.human.pdf")
	ggplot(Prom_vs_expr,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('H3K18la')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() # + ylim(c(2.5,12))
	ggplot(Prom_vs_expr,aes(y=H3K27ac,x=RNA)) + geom_point(alpha = 0.3, color='orange') + ggtitle('H3K27ac')+ylab("H3K27ac promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(Prom_vs_expr,aes(y=H3K27me3,x=RNA)) + geom_point(alpha = 0.3, color='darkred') + ggtitle('H3K27me3')+ylab("H3K27me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(Prom_vs_expr,aes(y=H3K4me3,x=RNA)) + geom_point(alpha = 0.3, color='darkgreen') + ggtitle('H3K4me3')+ylab("H3K4me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(Prom_vs_expr,aes(y=H3K9me3,x=RNA)) + geom_point(alpha = 0.3, color='grey') + ggtitle('H3K9me3')+ylab("H3K9me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() 

	dev.off()
	
# Calculate spearman correlation

cor.test(Prom_vs_expr$H3K4me3, Prom_vs_expr$RNA, method='spearman')
cor.test(Prom_vs_expr$H3K27me3, Prom_vs_expr$RNA, method='spearman')
cor.test(Prom_vs_expr$H3K9me3, Prom_vs_expr$RNA, method='spearman')
cor.test(Prom_vs_expr$H3K27ac, Prom_vs_expr$RNA, method='spearman')
cor.test(Prom_vs_expr$H3K18la, Prom_vs_expr$RNA, method='spearman')

# Define highest genes expressed / genes with highest promoter hPTM values

highest.2000.genes <- Prom_vs_expr[order(Prom_vs_expr$RNA,decreasing=T),][1:2000,]
highest.2000.promoter.H3K18la <- Prom_vs_expr[order(Prom_vs_expr$H3K18la,decreasing=T),][1:2000,]
highest.2000.promoter.H3K27ac <- Prom_vs_expr[order(Prom_vs_expr$H3K27ac,decreasing=T),][1:2000,]
highest.2000.promoter.H3K4me3 <- Prom_vs_expr[order(Prom_vs_expr$H3K4me3,decreasing=T),][1:2000,]
highest.2000.promoter.H3K27me3 <- Prom_vs_expr[order(Prom_vs_expr$H3K27me3,decreasing=T),][1:2000,]
highest.2000.promoter.H3K9me3 <- Prom_vs_expr[order(Prom_vs_expr$H3K9me3,decreasing=T),][1:2000,]

# Perform GO of highest genes expressed / genes with highest promoter hPTM values

GO.RNA <- enrichGO(highest.2000.genes$EntrezID,pvalueCutoff  = 0.05,
	                         pAdjustMethod = "BH",
	                         OrgDb='org.Hs.eg.db',
	                         ont="ALL")
	 
GO.promoter.H3K18la <- enrichGO(highest.2000.promoter.H3K18la$EntrezID,pvalueCutoff  = 0.05,
	 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.promoter.H3K27ac <- enrichGO(highest.2000.promoter.H3K27ac$EntrezID,pvalueCutoff  = 0.05,
	 	 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.promoter.H3K4me3 <- enrichGO(highest.2000.promoter.H3K4me3$EntrezID,pvalueCutoff  = 0.05,
		 	pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.promoter.H3K9me3 <- enrichGO(highest.2000.promoter.H3K9me3$EntrezID,pvalueCutoff  = 0.05,
		 	pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.promoter.H3K27me3 <- enrichGO(highest.2000.promoter.H3K27me3$EntrezID,pvalueCutoff  = 0.05,
		 	 	 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")

library(clusterProfiler)
pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/SupplFig4E.F.GO human promoter plot.15.terms.pdf",width=10,height=5)	 
dotplot(GO.RNA,showCategory=15,title='top 2000 expressed genes') 
dotplot(GO.promoter.H3K18la,showCategory=15,title='top H3K18la promoters') 
dotplot(GO.promoter.H3K27ac,showCategory=15,title='top H3K27ac promoters') 
dotplot(GO.promoter.H3K4me3,showCategory=15,title='top H3K4me3 promoters') 
dotplot(GO.promoter.H3K27me3,showCategory=15,title='top H3K27me3 promoters') 
dotplot(GO.promoter.H3K9me3,showCategory=15,title='top H3K9me3 promoters') 
dev.off()


######## CGI promoters

CGI_Prom_vs_expr <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Human/hPTM_RNA_CGIpromoter/hsMuscle_hPTM_RNA_CGIpr_normCounts.txt")

pdf("/home/egalle/public/_Projects/MS_H3K18la/GitHub/try_output/SupplFig4G.cor.hPTM.CGI.promoter.vs.RNA.human.pdf")
	ggplot(CGI_Prom_vs_expr,aes(y=H3K18la,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('H3K18la')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() # + ylim(c(2.5,12))
	ggplot(Prom_vs_expr,aes(y=H3K27ac,x=RNA)) + geom_point(alpha = 0.3, color='orange') + ggtitle('H3K27ac')+ylab("H3K27ac promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(CGI_Prom_vs_expr,aes(y=H3K27me3,x=RNA)) + geom_point(alpha = 0.3, color='darkred') + ggtitle('H3K27me3')+ylab("H3K27me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(CGI_Prom_vs_expr,aes(y=H3K4me3,x=RNA)) + geom_point(alpha = 0.3, color='darkgreen') + ggtitle('H3K4me3')+ylab("H3K4me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() #+ ylim(c(2.5,12))
	ggplot(CGI_Prom_vs_expr,aes(y=H3K9me3,x=RNA)) + geom_point(alpha = 0.3, color='grey') + ggtitle('H3K9me3')+ylab("H3K9me3 promoter peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() 

	dev.off()
	
# Calculate spearman correlation

cor.test(CGI_Prom_vs_expr$H3K4me3, Prom_vs_expr$RNA, method='spearman')
cor.test(CGI_Prom_vs_expr$H3K27me3, Prom_vs_expr$RNA, method='spearman')
cor.test(CGI_Prom_vs_expr$H3K9me3, Prom_vs_expr$RNA, method='spearman')
cor.test(CGI_Prom_vs_expr$H3K27ac, Prom_vs_expr$RNA, method='spearman')
cor.test(CGI_Prom_vs_expr$H3K18la, Prom_vs_expr$RNA, method='spearman')

######## suppl. fig 4H: Adhi?
# Figure 3A/C/D/H/I: Adhi


##### Link peak-enhancers-overlap files to closest non-overlapping gene promoter in bash

QUERYbed=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.human.bed
refbed1=/home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Human/Gene_promoter_3000_hg38_Seqmonk.bed
refbed2=/home/egalle/public/EvaGalle/Data/Public_Data/Williams_2021_human.muscle/human_muscle_enhancers.bed

sort -k1,1 -k2,2n $refbed1 > $refbed1.sorted.bed
sort -k1,1 -k2,2n $refbed2 > $refbed2.sorted.bed
n=1
bedtools closest -io -k $n -a $refbed2.sorted.bed -b $refbed1.sorted.bed -D ref > /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Human/muscle.enhancers.${n}.closest.non.overlapping.promoters.bed
n=5
bedtools closest -io -k $n -a $refbed2.sorted.bed -b $refbed1.sorted.bed -D ref > /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Human/muscle.enhancers.${n}.closest.non.overlapping.promoters.bed

peakPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks_human/peaks/
cd $peakPath

for i in Muscle*
do
resultPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/human
n=1
bedtools intersect -wao -a /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Human/muscle.enhancers.${n}.closest.non.overlapping.promoters.bed -b ${i} > ${resultPath}/${i}.overlap_muscle.enhancers.${n}.closest.non.overlapping.promoters.bed.bed
n=5
bedtools intersect -wao -a /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Human/muscle.enhancers.${n}.closest.non.overlapping.promoters.bed -b ${i} > ${resultPath}/${i}.overlap_muscle.enhancers.${n}.closest.non.overlapping.promoters.bed.bed
done	

# switch to R

R
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/human/")
xi <-list()
xnames <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/human",pattern='Muscle')

for (i in 1:length(xnames)){
	xi[[i]] <- read.delim(xnames[i],header=F)
	colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','width.enhancer','logFC (pre- vs post-training)',
	'logCPM	p-value	FDR','p-value','FDR','linked.promoter.chr','linked.promoter.start','linked.promoter.end','strand.promoter','gene.id','ensembl.id','gene.description','gene.strand',
	'distance.enhancer.to.promoter','peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
	xi[[i]]$peak.id <- paste0(xi[[i]]$peak.chr,':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
}
names(xi) <-xnames

rawPeakCounts.H3K18la <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_H3K18la_rawPeakCounts.rds")
rawPeakCounts.H3K4me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_H3K4me3_rawPeakCounts.rds")
rawPeakCounts.H3K27me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_H3K27me3_rawPeakCounts.rds")
rawPeakCounts.H3K27ac <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_H3K27ac_rawPeakCounts.rds")
rawPeakCounts.H3K9me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_H3K9me3_rawPeakCounts.rds")
normPeakCounts.H3K18la <- edgeR::cpm((rawPeakCounts.H3K18la[,1]+rawPeakCounts.H3K18la[,2])/2, log=TRUE)
normPeakCounts.H3K4me3 <- edgeR::cpm((rawPeakCounts.H3K4me3[,1]+rawPeakCounts.H3K4me3[,2])/2, log=TRUE)
normPeakCounts.H3K27me3 <- edgeR::cpm((rawPeakCounts.H3K27me3[,1]+rawPeakCounts.H3K27me3[,2])/2, log=TRUE)
normPeakCounts.H3K9me3 <- edgeR::cpm((rawPeakCounts.H3K9me3[,1]+rawPeakCounts.H3K9me3[,2])/2, log=TRUE)
normPeakCounts.H3K27ac <- edgeR::cpm((rawPeakCounts.H3K27ac[,1]+rawPeakCounts.H3K27ac[,2])/2, log=TRUE)


normPeakCounts <- list(normPeakCounts.H3K18la,
	normPeakCounts.H3K27ac,
	normPeakCounts.H3K27me3,
	normPeakCounts.H3K4me3,
	normPeakCounts.H3K9me3)
names(normPeakCounts) <- c('normPeakCounts.H3K18la','normPeakCounts.H3K27ac','normPeakCounts.H3K27me3','normPeakCounts.H3K4me3','normPeakCounts.H3K9me3')

for (i in c(6,12,18,24,30)){
	xi[[i]] <- merge(xi[[i]], normPeakCounts[[i/6]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
	colnames(xi[[i]])[length(colnames(xi[[i]]))] <- 'peak.log2CPM'
}

for (i in c(6,12,18,24,30)){
	xi[[i-1]] <- merge(xi[[i-1]], normPeakCounts[[i/6]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
	colnames(xi[[i-1]])[length(colnames(xi[[i-1]]))] <- 'peak.log2CPM'
}


RNA_and_PTM.Promoter <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Human/hsMuscle_hist_RNA_comb_all.txt")

x<-list()
for (i in 1:length(xnames)){
	x[[i]] <- merge(xi[[i]],RNA_and_PTM.Promoter, by.x='ensembl.id',by.y='ensembl_gene_id', all.x=T,all.y=T)
}
x.copy <- x
lapply(x,dim)

library(ggplot2)
#SupplFig5.A : Adhi
pdf("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/plots/SupplFig5.B.cor.hPTM.enhancer.vs.RNA.human.pdf")
	ggplot(x[[6]],aes(y=peak.log2CPM,x=RNA)) + geom_point(alpha = 0.3, color='blue') + ggtitle('H3K18la')+ylab("H3K18la enhancer peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() + ylim(c(2.5,12))
	ggplot(x[[12]],aes(y=peak.log2CPM,x=RNA)) + geom_point(alpha = 0.3, color='orange') + ggtitle('H3K27ac')+ylab("H3K27ac enhancer peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() + ylim(c(2.5,12))
	ggplot(x[[18]],aes(y=peak.log2CPM,x=RNA)) + geom_point(alpha = 0.3, color='darkred') + ggtitle('H3K27me3')+ylab("H3K27me3 enhancer peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() + ylim(c(2.5,12))
	ggplot(x[[24]],aes(y=peak.log2CPM,x=RNA)) + geom_point(alpha = 0.3, color='darkgreen') + ggtitle('H3K4me3')+ylab("H3K4me3 enhancer peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() + ylim(c(2.5,12))
	ggplot(x[[30]],aes(y=peak.log2CPM,x=RNA)) + geom_point(alpha = 0.3, color='grey') + ggtitle('H3K9me3')+ylab("H3K9me3 enhancer peak (log2 CPM)") + xlab("Gene expression (log2 CPM)") + theme_bw() + ylim(c(2.5,12))

	dev.off()


Rho.gene.expression.enhancer.hPTM <- rbind(c(cor.test(x[[6]]$peak.log2CPM,x[[6]]$RNA,method='spearman')[[4]],
cor.test(x[[12]]$peak.log2CPM,x[[12]]$RNA,method='spearman')[[4]],
cor.test(x[[18]]$peak.log2CPM,x[[18]]$RNA,method='spearman')[[4]],
cor.test(x[[24]]$peak.log2CPM,x[[24]]$RNA,method='spearman')[[4]],
cor.test(x[[30]]$peak.log2CPM,x[[30]]$RNA,method='spearman')[[4]]),
c(cor.test(x[[6-1]]$peak.log2CPM,x[[6-1]]$RNA,method='spearman')[[4]],
cor.test(x[[12-1]]$peak.log2CPM,x[[12-1]]$RNA,method='spearman')[[4]],
cor.test(x[[18-1]]$peak.log2CPM,x[[18-1]]$RNA,method='spearman')[[4]],
cor.test(x[[24-1]]$peak.log2CPM,x[[24-1]]$RNA,method='spearman')[[4]],
cor.test(x[[30-1]]$peak.log2CPM,x[[30-1]]$RNA,method='spearman')[[4]]))

rownames(Rho.gene.expression.enhancer.hPTM) <- c('1 linked gene','5 linked genes')
colnames(Rho.gene.expression.enhancer.hPTM) <- c('H3K18la','H3K27ac','H3K27me3','H3K4me3','H3K9me3')


highest.2000.enhancer.H3K18la <-x[[6]][order(x[[6]]$peak.log2CPM,decreasing=T),][1:2000,]
highest.2000.enhancer.H3K27ac <-x[[12]][order(x[[12]]$peak.log2CPM,decreasing=T),][1:2000,]
highest.2000.enhancer.H3K4me3 <-x[[18]][order(x[[18]]$peak.log2CPM,decreasing=T),][1:2000,]
highest.2000.enhancer.H3K27me3 <-x[[24]][order(x[[24]]$peak.log2CPM,decreasing=T),][1:2000,]
highest.2000.enhancer.H3K9me3 <-x[[30]][order(x[[30]]$peak.log2CPM,decreasing=T),][1:2000,]

GO.enhancer.H3K18la <- enrichGO(highest.2000.enhancer.H3K18la$entrezgene_id,pvalueCutoff  = 0.05,
	 	 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.enhancer.H3K27ac <- enrichGO(highest.2000.enhancer.H3K27ac$entrezgene_id,pvalueCutoff  = 0.05,
		 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.enhancer.H3K4me3 <- enrichGO(highest.2000.enhancer.H3K4me3$entrezgene_id,pvalueCutoff  = 0.05,
		 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.enhancer.H3K27me3 <- enrichGO(highest.2000.enhancer.H3K27me3$entrezgene_id,pvalueCutoff  = 0.05,
		 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")
GO.enhancer.H3K9me3 <- enrichGO(highest.2000.enhancer.H3K9me3$entrezgene_id,pvalueCutoff  = 0.05,
		 pAdjustMethod = "BH",OrgDb='org.Hs.eg.db', ont="ALL")

pdf("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/plots/Fig3.J.GO human plot.25.terms.pdf",width=10,height=7)	 
dotplot(GO.enhancer.H3K18la,showCategory=25,title='top H3K18la enhancers') 
dotplot(GO.enhancer.H3K27ac,showCategory=25,title='top H3K27ac enhancers') 
dotplot(GO.enhancer.H3K4me3,showCategory=25,title='top H3K4me3 enhancers') 
dotplot(GO.enhancer.H3K27me3,showCategory=25,title='top H3K27me3 enhancers') 
dotplot(GO.enhancer.H3K9me3,showCategory=25,title='top H3K9me3 enhancers') 
dev.off()
