
# First: compute overlap between peaks and cCRE in bash:
# define input
resultsPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE
refbed1=/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_mm10/mm10-ccREs.bed.sorted.bed
peakPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks/
peakPath1=/home/egalle/public/_Projects/MS_H3K18la/Mouse/peaks/
peakPath2=/home/egalle/public/_Projects/MS_H3K18la/Mouse/mergedPeaks/

cd $peakPath
for i in *bam*
do
# add 'chr' to chromosome column
    awk '{print "chr"$0}' $peakPath2/$i > $peakPath/$i.w.chr.bed
# sort bed files to speed up the computation
    sort -k1,1 -k2,2n $peakPath/$i.w.chr.bed > $peakPath/$i.w.chr.sorted.bed
# use bedtools intersect function to calculate overlap
    cd $resultsPath
    bedtools intersect -wao -a $peakPath/$i.w.chr.sorted.bed -b /home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_mm10/mm10-ccREs.bed.sorted.bed -sorted > $i.overlap.cCREs.bed
done

# analyse in R
R

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
library('nVennR')

# Read cCRE
cCREs <- read.table("/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_mm10/mm10-ccREs.bed.sorted.bed")
colnames(cCREs) <- c('chr.cCRE','start.cCRE','end.cCRE','id1.cCRE','id2.cCRE','type.cCRE')

# Read peak-cCRE-overlaps and peaks
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks")
peaks.files <- list.files(pattern='bed.w.chr.sorted.bed')
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE/overlap.cCRE/")
peaks.overlap.cCREs.files <- list.files(pattern='overlap.cCREs.bed')
#metadata <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/metadata_Eva.txt")
peaks.overlap.cCRE <- list()
peaks <- list()
for (i in 1:length(peaks.files)){
  peaks.overlap.cCRE[[i]] <- read.delim(paste0("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/cCRE/overlap.cCRE/",peaks.overlap.cCREs.files[i]),header=F)
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
peaks[[i]] <- read.table(paste0("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks/",peaks.files[i]))
colnames(peaks[[i]]) <-colnames(peaks.overlap.cCRE[[i]])[1:6]
}

# Rename peaks and overlap-files to remove superfluous name parts
names(peaks)<-gsub("_0.01.stringent.bed","",peaks.files)
names(peaks)<-gsub(".w.chr.sorted.bed","",names(peaks))
names(peaks)<-gsub(".0.01.stringent.bed","",names(peaks))

names(peaks.overlap.cCRE)<-gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)
names(peaks.overlap.cCRE)<-gsub("_0.01.stringent.bed.w.chr.sorted.bed.overlap.cCREs.bed","",names(peaks.overlap.cCRE))
names(peaks.overlap.cCRE)<-gsub(".0.01.stringent.bed.overlap.cCREs.bed","",names(peaks.overlap.cCRE))


# Calculate total peak bp for normalisation
total.peak.bp<-list()
for (i in 1:length(peaks)){
    peaks[[i]]$width <- peaks[[i]]$end.peak-peaks[[i]]$start.peak
    total.peak.bp[[i]]<-sum(peaks[[i]]$width )}

total.peak.bp<-unlist(total.peak.bp)
names(total.peak.bp)<-names(peaks)
write.table(total.peak.bp,"/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/total.peak.bp.txt",sep="\t")

# Compute lists of peaks that overlap with each cCRE 
peaks.w.cCRE.mm10 <- list()
peaks.bp.overlap.w.cCRE.mm10 <-list()
library(dplyr)
for (i in 1:length(peaks.overlap.cCRE)){
  peaks.oi<-peaks.overlap.cCRE[[i]]
  peaks.oi$bp.overlap <- as.numeric(peaks.oi$bp.overlap)
  peaks.bp.overlap.w.cCRE.mm10[[i]]<- data.frame(peaks.oi %>% group_by(type.cCRE) %>% summarize(sum=sum(bp.overlap)))
  
  peaks.w.cCRE.mm10[[i]]<-list(unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('PLS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('pELS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('dELS',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('CTCF-only,CTCF-bound',peaks.oi$type.cCRE)]),
                               unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[grep('H3K4me3',peaks.oi$type.cCRE)]),
                                unique(paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak)[peaks.oi$bp.overlap==0]))
                              # unique(paste0(peaks[[i]]$chr.peak,'.',peaks[[i]]$start.peak,'.',peaks[[i]]$end.peak))[!c(paste0(peaks[[i]]$chr.peak,'.',peaks[[i]]$start.peak,'.',peaks[[i]]$end.peak) %in% paste0(peaks.oi$chr.peak,'.',peaks.oi$start.peak,'.',peaks.oi$end.peak))])
  names(peaks.w.cCRE.mm10[[i]]) <- c('PLS', 'pELS','dELS','CTCF','DNase-H3K4me3','no cCRE')
}
names(peaks.w.cCRE.mm10) <-gsub("_0.01.stringent.bed.overlap.cCREs","",peaks.overlap.cCREs.files)
names(peaks.bp.overlap.w.cCRE.mm10) <-gsub("_0.01.stringent.bed.overlap.cCREs","",peaks.overlap.cCREs.files)

# normalise overlap bp for each pair of cCRE-peaks for total bp in specific cCRE and total bp in specific peak set
# create data frame with overlap bp for each pair of cCRE (rows) - peaks (columns)
peaks.bp.overlap.w.cCRE.mm102 <- bind_cols(peaks.bp.overlap.w.cCRE.mm10)
# remove redundant columns with cCRE annotations and add the latter as rownames
peaks.bp.overlap.w.cCRE.mm102 <-peaks.bp.overlap.w.cCRE.mm102[,c(1:c(dim(peaks.bp.overlap.w.cCRE.mm102)[2]/2))*2]
rownames(peaks.bp.overlap.w.cCRE.mm102)<-peaks.bp.overlap.w.cCRE.mm10[[1]]$type.cCRE
# Simplify colnames
colnames(peaks.bp.overlap.w.cCRE.mm102) <- gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)

# Load total peak bp sums
total.peak.bp<-read.table("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/total.peak.bp.txt",sep="\t")
total.peak.bp<-data.frame(rownames(total.peak.bp),total.peak.bp$x)
colnames(total.peak.bp) <- c('peak','bp')
total.peak.bp$peak <- gsub(".w.chr.sorted.bed","",total.peak.bp$peak)
total.peak.bp <- total.peak.bp[order(total.peak.bp$peak),]
rownames(total.peak.bp)<-total.peak.bp$peak

# calculate total number of bp per peak set that is overlapping with (any) cCRE
total.bp.in.cCRE <- colSums(peaks.bp.overlap.w.cCRE.mm102)

# add total number of bp per peak set that is overlapping with (any) cCRE and total peak bp sums to working dataframe
peaks.bp.overlap.w.cCRE.mm102 <-rbind(peaks.bp.overlap.w.cCRE.mm102,total.bp.in.cCRE,total.peak.bp$bp)
rownames(peaks.bp.overlap.w.cCRE.mm102)[11:12] <- c("total.bp.in.cCRE","total.peak.bp")

# normalise overlap bp for each pair of cCRE-peaks for total bp in specific peak set
peaks.bp.overlap.w.cCRE.mm103 <-t(t(as.matrix(peaks.bp.overlap.w.cCRE.mm102))/as.numeric(peaks.bp.overlap.w.cCRE.mm102[12,]))

# Calculate the fraction of bp in peaks that are not lying in cCRE and add to working dataframe
peaks.bp.overlap.w.cCRE.mm103 <- rbind(peaks.bp.overlap.w.cCRE.mm103, peaks.bp.overlap.w.cCRE.mm103[12,]-peaks.bp.overlap.w.cCRE.mm103[11,])
rownames(peaks.bp.overlap.w.cCRE.mm103)[13] <- "peak.bp.not.in.cCRE"

# calculate total number of bp per cCRE category
cCRE.per.type<- data.frame(cCREs %>% group_by(type.cCRE) %>% summarize(sum=sum(end.cCRE-start.cCRE)))

# normalise overlap bp for each pair of cCRE-peaks for total bp per cCRE category
cCRE.per.type <- cCRE.per.type[order(cCRE.per.type$type.cCRE),]
peaks.bp.overlap.w.cCRE.mm104 <-as.matrix(peaks.bp.overlap.w.cCRE.mm103)/c(1,as.numeric(cCRE.per.type$sum),1,sum(as.numeric(cCRE.per.type$sum)),1)

# make subset dataframe with only non-duplicate-information
peaks.bp.overlap.w.cCRE.mm105 <- peaks.bp.overlap.w.cCRE.mm104[c(2:10),grep('merged|MB|BMDM',colnames(peaks.bp.overlap.w.cCRE.mm104))]

# make subset dataframe where major cCRE categories are taken together
peaks.bp.overlap.w.cCRE.mm106 <- rbind(peaks.bp.overlap.w.cCRE.mm105[1,],
    peaks.bp.overlap.w.cCRE.mm105[2,]+peaks.bp.overlap.w.cCRE.mm105[3,],
    peaks.bp.overlap.w.cCRE.mm105[4,]+peaks.bp.overlap.w.cCRE.mm105[5,],
    peaks.bp.overlap.w.cCRE.mm105[6,]+peaks.bp.overlap.w.cCRE.mm105[7,],
    peaks.bp.overlap.w.cCRE.mm105[8,]+ peaks.bp.overlap.w.cCRE.mm105[9,])
rownames(peaks.bp.overlap.w.cCRE.mm106) <- c('CTCF-only','dELS','DNase-H3K4me3','pELS','PLS')

# plot enrichment
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-2, 2, by = 0.1)
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig2B.cCRE.mouse.heatmap.pdf",width=8,height=8)
pheatmap(peaks.bp.overlap.w.cCRE.mm106,cellheight=10,cellwidth=10, scale='column',breaks = breaksList,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)))
dev.off()






# Fig 2C:

# Create a masterpeak file to compute overlap with tissue-specific enhancers:

masterPeak_df <- bind_rows(peaks, .id = "column_label")
masterPeak_df <- masterPeak_df[,c(2,3,4,5,6,7,8,1)]
write.table(masterPeak_df,"/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.bed",quote=F,col.names=F,sep="\t",row.names=F)

# Calculate overlap with tissue specific enhancers in bash

resultsPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers
QUERYbed=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.bed
sort -k1,1 -k2,2n $QUERYbed > $QUERYbed.sorted.bed

# calculate overlap for all merged peaks assembled in 1 data frame:

# define paths to refbeds containing sorted enhancer coordinates for the different tissues of interest
# for ESC's, we use ENCODE's dataset of cCREs which we'll subset further down to select ELS specifically.
refbed2=/home/egalle/public/EvaGalle/Data/Public_Data/Rovito_2021_NAR_mouse.muscle/enhancers_gastrocnemius.txt
refbed3=/home/egalle/public/EvaGalle/Data/Public_Data/Denisenko_2017_Epigeneticsandchromatin/BMDM.activation.enhancers.bed
refbed4=/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_mm10/mm10-ccREs.E14/ENCFF872IES_ENCFF857GJE_ENCFF163HEV.7group.bed
refbed5=/home/egalle/public/EvaGalle/Data/Public_Data/Blum_2012_GenesDev_myoblasts.myotubes/MB_enhancers.mm10.bed.sorted.txt
refbed6=/home/egalle/public/EvaGalle/Data/Public_Data/Blum_2012_GenesDev_myoblasts.myotubes/MT_enhancers.mm10.bed.sorted.txt
refbed7=/home/egalle/public/EvaGalle/Data/Public_Data/Rovito_2021_NAR_mouse.muscle/enhancers_adiposetissue.txt

# calculate overlap with bedtools intersect function
cd $resultsPath
bedtools intersect -wao -a $refbed2 -b $QUERYbed.sorted.bed > gastroc.enhancers.covered.by.peaks.bed
bedtools intersect -wao -a $refbed3 -b $QUERYbed.sorted.bed > BMDM.Denisenko.enhancers.covered.by.peaks.bed
bedtools intersect -wao -a $refbed4 -b $QUERYbed.sorted.bed > E14.cCREs.covered.by.peaks.bed
bedtools intersect -wao -a $refbed7 -b $QUERYbed.sorted.bed > ADI.enhancers.covered.by.peaks.bed
bedtools intersect -wao -a $refbed5 -b $QUERYbed.sorted.bed > MB.enhancers.covered.by.peaks.bed
bedtools intersect -wao -a $refbed6 -b $QUERYbed.sorted.bed > MT.enhancers.covered.by.peaks.bed

# Process in R: calculate fraction of tissue-specific enhancers that is covered by each peak set
# GAS
R
setwd( "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers")
gastroc.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/gastroc.enhancers.covered.by.peaks.bed",header=F)
colnames(gastroc.enhancers.covered.by.peaks) <- c('chr.enhancer','start.enhancer','end.enhancer','enhancer.nr','peak.chr','peak.start','peak.end',
'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
gastroc.enhancers.covered <- c()
for (i in 1:length(unique(gastroc.enhancers.covered.by.peaks$peak.origin))){
    peakOI<-unique(gastroc.enhancers.covered.by.peaks$peak.origin)[i]
    gastroc.enhancers.covered[i]<-length(unique(gastroc.enhancers.covered.by.peaks$enhancer.nr[gastroc.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(gastroc.enhancers.covered.by.peaks$enhancer.nr))
    
}
names(gastroc.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(gastroc.enhancers.covered.by.peaks$peak.origin))
names(gastroc.enhancers.covered)[names(gastroc.enhancers.covered)=='.']<-'Not covered by any peak'

# ADIP
setwd( "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers")
ADI.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/ADI.enhancers.covered.by.peaks.bed",header=F)
colnames(ADI.enhancers.covered.by.peaks) <- c('chr.enhancer','start.enhancer','end.enhancer','enhancer.nr','peak.chr','peak.start','peak.end',
'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
ADI.enhancers.covered <- c()
for (i in 1:length(unique(ADI.enhancers.covered.by.peaks$peak.origin))){
    peakOI<-unique(ADI.enhancers.covered.by.peaks$peak.origin)[i]
    ADI.enhancers.covered[i]<-length(unique(ADI.enhancers.covered.by.peaks$enhancer.nr[ADI.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(ADI.enhancers.covered.by.peaks$enhancer.nr))
   
}
names(ADI.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(ADI.enhancers.covered.by.peaks$peak.origin))
names(ADI.enhancers.covered)[names(ADI.enhancers.covered)=='.']<-'Not covered by any peak'

# MB/MT

setwd( "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers")
MB.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/MB.enhancers.covered.by.peaks.bed",header=F)
colnames(MB.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.id','peak.chr','peak.start','peak.end',
'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
MT.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/MT.enhancers.covered.by.peaks.bed",header=F)
colnames(MT.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.id','peak.chr','peak.start','peak.end',
'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
        
MT.enhancers.covered <- c()
MB.enhancers.covered <- c()
for (i in 1:length(unique(MB.enhancers.covered.by.peaks$peak.origin))){
    peakOI<-unique(MB.enhancers.covered.by.peaks$peak.origin)[i]
    MB.enhancers.covered[i]<-length(unique(MB.enhancers.covered.by.peaks$enhancer.id[MB.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(MB.enhancers.covered.by.peaks$enhancer.id))
    peakOI<-unique(MT.enhancers.covered.by.peaks$peak.origin)[i]
    MT.enhancers.covered[i]<-length(unique(MT.enhancers.covered.by.peaks$enhancer.id[MT.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(MT.enhancers.covered.by.peaks$enhancer.id))
    
}
names(MB.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(MB.enhancers.covered.by.peaks$peak.origin))
names(MT.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(MT.enhancers.covered.by.peaks$peak.origin))
names(MB.enhancers.covered)[names(MB.enhancers.covered)=='.']<-'Not covered by any peak'
names(MT.enhancers.covered)[names(MT.enhancers.covered)=='.']<-'Not covered by any peak'

# BMDM

BMDM.Denisenko.enhancers.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/BMDM.Denisenko.enhancers.covered.by.peaks.bed",header=F)
colnames(BMDM.Denisenko.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.target.chr','enhancer.target.start','enhancer.target.end','enhancer.strand','target.gene','chr.peak','start.peak','end.peak','peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')

BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr <- paste0(BMDM.Denisenko.enhancers.covered.by.peaks$chr.enhancer,BMDM.Denisenko.enhancers.covered.by.peaks$start.enhancer,BMDM.Denisenko.enhancers.covered.by.peaks$end.enhancer)
BMDM.Denisenko.enhancers.covered <- c()
for (i in 1:length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin))){
    peakOI<-unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin)[i]
    BMDM.Denisenko.enhancers.covered[i]<-length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr[BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr))
    
}
names(BMDM.Denisenko.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin))
BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin[BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin=='.'] <- 'Not covered by any peak'
names(BMDM.Denisenko.enhancers.covered)[names(BMDM.Denisenko.enhancers.covered)=='.'] <- 'Not covered by any peak'

# ESC

E14.cCREs.covered.by.peaks <- read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/tissue_enhancers/E14.cCREs.covered.by.peaks.bed",header=F)
colnames(E14.cCREs.covered.by.peaks) <- c("chr.cCRE",'start.cCRE','end.cCRE','cCRE.id',rep('NA',5),'type.cCRE','NA',
'chr.peak','start.peak','end.peak','peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
E14.cCREs.covered.by.peaks<-E14.cCREs.covered.by.peaks[,-grep('NA',colnames(E14.cCREs.covered.by.peaks))]


E14.enhancers.covered.by.peaks<-E14.cCREs.covered.by.peaks[grep('ELS',E14.cCREs.covered.by.peaks$type.cCRE),]
E14.enhancers.covered <- c()
for (i in 1:length(unique(E14.enhancers.covered.by.peaks$peak.origin))){
    peakOI<-unique(E14.enhancers.covered.by.peaks$peak.origin)[i]
    E14.enhancers.covered[i]<-length(unique(E14.enhancers.covered.by.peaks$cCRE.id[E14.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(E14.enhancers.covered.by.peaks$cCRE.id))
}
names(E14.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(E14.enhancers.covered.by.peaks$peak.origin))
E14.enhancers.covered.by.peaks$peak.origin[E14.enhancers.covered.by.peaks$peak.origin=='.'] <- 'Not covered by any peak'

enhancers.covered <- cbind(BMDM.Denisenko.enhancers.covered[order(names(BMDM.Denisenko.enhancers.covered))], 
E14.enhancers.covered[order(names(E14.enhancers.covered))],
gastroc.enhancers.covered[order(names(gastroc.enhancers.covered))],
ADI.enhancers.covered[order(names(ADI.enhancers.covered))],
MB.enhancers.covered[order(names(MB.enhancers.covered))],
MT.enhancers.covered[order(names(MT.enhancers.covered))])
colnames(enhancers.covered) <- c('BMDM','E14','GAS','ADI','MB','MT')

enhancers.covered.selection <- enhancers.covered[grep('merged|ESC.ser_H3K27ac|BMDM|MB|PIM_H3K27ac',rownames(enhancers.covered)),]
cols<-cbind(rownames(enhancers.covered.selection),NA)
cols[,2][grep('H3K18la',cols[,1])]<-'blue'
cols[,2][grep('H3K27me3',cols[,1])]<-'darkred'
cols[,2][grep('H3K4me3',cols[,1])]<-'darkgreen'
cols[,2][grep('ac',cols[,1])]<-'orange'

library(tidyr)
temp <- data.frame(x=rownames(enhancers.covered.selection))
tissue <- temp %>% separate(x, c('Tissue','Histone','Rep'), sep ='_')

pdf('/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig2C.Fraction of tissue enhancers covered by peaks.pdf',height=6,width=25)
par(mar=c(5,5,2,1),mfrow=c(1,6),cex=1)
for (i in 1:dim(enhancers.covered.selection)[2]){
	barplot(enhancers.covered.selection[,i][order(enhancers.covered.selection[,i])],horiz=T,
        las=1,xlim=c(0,1),col=cols[,2][order(enhancers.covered.selection[,i])],
        names = tissue[,1][order(enhancers.covered.selection[,i])],
        xlab=paste0('Fraction of ',colnames(enhancers.covered.selection)[i],' enhancers covered'))
}
dev.off()


# Fig 2D

# list which peaks are overlapping with the resp. cCRE per tissue
x.dELS <-list()
x.pELS <-list()
x.PLS <-list()
for (i in c(1:length(peaks.overlap.cCRE))){
	peaks.overlap.cCRE[[i]]$enhancer.id <- paste0(peaks.overlap.cCRE[[i]]$chr.cCRE,":",
													peaks.overlap.cCRE[[i]]$start.cCRE,"-",
													peaks.overlap.cCRE[[i]]$end.cCRE)
  x.dELS[[i]] <- peaks.overlap.cCRE[[i]]$enhancer.id[peaks.overlap.cCRE[[i]]$type.cCRE=='dELS']
  x.pELS[[i]] <- peaks.overlap.cCRE[[i]]$enhancer.id[peaks.overlap.cCRE[[i]]$type.cCRE=='pELS']
  x.PLS[[i]] <- peaks.overlap.cCRE[[i]]$enhancer.id[peaks.overlap.cCRE[[i]]$type.cCRE=='PLS']
  }
  names(x.dELS) <- names(peaks.overlap.cCRE)
  names(x.pELS) <- names(peaks.overlap.cCRE)
  names(x.PLS) <- names(peaks.overlap.cCRE)
  
  
# Create venn plots with the overlaps of interest
ggvennlist.cCRE<-list()
x<-x[order(names(x))]

  ggvennlist.cCRE[[1]]<-ggVennDiagram(x.dELS[names(x.dELS) %in% c('PIM_H3K18la_merged','PIM_H3K4me3_merged',"PIM_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('dELS PIM')+ scale_color_manual(values=c('blue','orange','darkgreen'))

  ggvennlist.cCRE[[2]]<-ggVennDiagram(x.pELS[names(x.pELS) %in% c('PIM_H3K18la_merged','PIM_H3K4me3_merged',"PIM_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('pELS PIM')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								
  ggvennlist.cCRE[[3]]<-ggVennDiagram(x.PLS[names(x.PLS) %in% c('PIM_H3K18la_merged','PIM_H3K4me3_merged',"PIM_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('PLS PIM')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								
  ggvennlist.cCRE[[4]]<-ggVennDiagram(x.dELS[names(x.dELS) %in% c('GAS_H3K18la_merged','GAS_H3K4me3_merged',"GAS_H3K27ac_merged")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('dELS GAS')+ scale_color_manual(values=c('blue','orange','darkgreen'))

  ggvennlist.cCRE[[5]]<-ggVennDiagram(x.pELS[names(x.pELS) %in% c('GAS_H3K18la_merged','GAS_H3K4me3_merged',"GAS_H3K27ac_merged")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('pELS GAS')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								
  ggvennlist.cCRE[[6]]<-ggVennDiagram(x.PLS[names(x.PLS) %in% c('GAS_H3K18la_merged','GAS_H3K4me3_merged',"GAS_H3K27ac_merged")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('PLS GAS')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								
  ggvennlist.cCRE[[7]]<-ggVennDiagram(x.dELS[names(x.dELS) %in% c('ESC.ser_H3K18la_merged','ESC.ser_H3K4me3_merged',"ESC.ser_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('dELS ESC.ser')+ scale_color_manual(values=c('blue','orange','darkgreen'))

  ggvennlist.cCRE[[8]]<-ggVennDiagram(x.pELS[names(x.pELS) %in% c('ESC.ser_H3K18la_merged','ESC.ser_H3K4me3_merged',"ESC.ser_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('pELS ESC.ser')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								
  ggvennlist.cCRE[[9]]<-ggVennDiagram(x.PLS[names(x.PLS) %in% c('ESC.ser_H3K18la_merged','ESC.ser_H3K4me3_merged',"ESC.ser_H3K27ac_1")],
  								label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  								ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  								theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  								ggtitle('PLS ESC.ser')+ scale_color_manual(values=c('blue','orange','darkgreen'))
  								

pdf('/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig2D.cCRE.covered.by.merged.TISSUE.peaks.venn.pdf',width=10,height=10)
do.call(grid.arrange,c(ggvennlist.cCRE[1:9],ncol=3))
dev.off()

# Fig 2E : Adhi

# Fig 2F+ SupplFig4C:

# Perform GO analysis on genes closest to the top 2000 H3K18lactylated dELS
# To do this, we must first link each enhancer to its closest gene
#       (to be precise: we linked each enhancer to its closest-non-overlapping gene promoter)

# list of sorted dELS
ME7=/home/egalle/public/EvaGalle/Data/Public_Data/ENCODE_cCRE_mm10/mm10-ccREs.dELS.bed.sorted.bed

# compute closest-non-overlapping (n) gene promoter(s) to each enhancer
n=1
mm10_promoters=/home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Mouse/Mus_musculus.GRCm38.102_ensDb.giga.promoters.protein.coding.sorted.bed
bedtools closest -io -k $n -a ${ME7} -b ${mm10_promoters} -D ref > /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Mouse/dELS.${n}.closest.non.overlapping.promoters.bed

# overlap these enhancers-linked-to-genes with peak files of the matching tissue
peakPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks
resultPath=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression
cd $peakPath

for i in *GAS*merged* *ESC.ser*merged* *PIM*merged* ESC.ser_H3K27ac* PIM_H3K27ac*
do
    bedtools intersect -wao -a /home/egalle/public/EvaGalle/Data/Public_Data/Genomic_Annotation_Data/Mouse/dELS.${n}.closest.non.overlapping.promoters.bed -b ${i} > ${resultPath}/dELS.${n}.closest.non.overlapping.promoters.overlap.$i
done


# Analyse in R
# Read bed files containing enhancers-linked-to-genes overlapped with peak files of the matching tissue
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/")
xi <-list()
xnames <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/",pattern='GAS')[1:4]

for (i in 1:length(xnames)){
    xi[[i]] <- read.delim(xnames[i],header=F)

     colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
     'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
     'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
    xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
    
}
names(xi) <-xnames

# read rawPeakCounts and normalise
# Adhi can you add how you made these?
rawPeakCounts.H3K18la <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/GAS_H3K18la_rawPeakCounts.rds")
rawPeakCounts.H3K4me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/GAS_H3K4me3_rawPeakCounts.rds")
rawPeakCounts.H3K27me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/GAS_H3K27me3_rawPeakCounts.rds")
normPeakCounts.H3K18la <- edgeR::cpm((rawPeakCounts.H3K18la[,1]+rawPeakCounts.H3K18la[,2])/2, log=TRUE)
normPeakCounts.H3K4me3 <- edgeR::cpm((rawPeakCounts.H3K4me3[,1]+rawPeakCounts.H3K4me3[,2])/2, log=TRUE)
normPeakCounts.H3K27me3 <- edgeR::cpm((rawPeakCounts.H3K27me3[,1]+rawPeakCounts.H3K27me3[,2])/2, log=TRUE)

normPeakCounts <- list(normPeakCounts.H3K18la,normPeakCounts.H3K18la,normPeakCounts.H3K27me3,normPeakCounts.H3K4me3)

for (i in 1:length(xnames)){
    xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
    colnames(xi[[i]])[colnames(xi[[i]])=='V1'] <- 'peak.log2CPM'
}

# Add gene expression info
RPKM_Counts <-read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/RNA/GAS_RPKM_TPM.txt")

x<-list()
for (i in 1:length(xnames)){
    x[[i]] <- merge(xi[[i]],RPKM_Counts, by.x='linked.ensemblid', by.y='Gene_ID',all.x=T,all.y=T)
    colnames(x[[i]])[length(colnames(x[[i]]))] <- 'gene.expression.log2CPM'
}
GAS <- x

# find GO terms of genes closest to highest 2000 enhancers per hPTM
highest.2000.enhancer.H3K18la <-GAS[[1]][order(GAS[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

library(clusterProfiler)
library(org.Mm.eg.db)
x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.GAS <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")
         

# Repeat for ESC
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/")
xi <-list()
xnames <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/",pattern='ESC')[1:4]

for (i in 1:length(xnames)){
    xi[[i]] <- read.delim(xnames[i],header=F)
    colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
     'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
     'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
    xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
    
}
names(xi) <-xnames



# read rawPeakCounts and normalise
# Adhi can you add how you made these?
rawPeakCounts.H3K18la <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/ESC_H3K18la_rawPeakCounts.rds")
rawPeakCounts.H3K4me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/ESC_H3K4me3_rawPeakCounts.rds")
rawPeakCounts.H3K27me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/ESC_H3K27me3_rawPeakCounts.rds")
normPeakCounts.H3K18la <- edgeR::cpm((rawPeakCounts.H3K18la[,1]+rawPeakCounts.H3K18la[,2])/2, log=TRUE)
normPeakCounts.H3K4me3 <- edgeR::cpm((rawPeakCounts.H3K4me3[,1]+rawPeakCounts.H3K4me3[,2])/2, log=TRUE)
normPeakCounts.H3K27me3 <- edgeR::cpm((rawPeakCounts.H3K27me3[,1]+rawPeakCounts.H3K27me3[,2])/2, log=TRUE)

normPeakCounts <- list(normPeakCounts.H3K18la,normPeakCounts.H3K18la,normPeakCounts.H3K27me3,normPeakCounts.H3K4me3)

for (i in 1:length(xnames)){
    xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
    colnames(xi[[i]])[colnames(xi[[i]])=='V1'] <- 'peak.log2CPM'
}

# Add gene expression info
rawCounts <- readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/RNA/S046_rawCounts.rds")
normCounts <- data.frame(edgeR::cpm(rowMeans(rawCounts[,grep('SL',colnames(rawCounts))]), log=TRUE))

RPKM_Counts <-read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/RNA/ESC.ser_ESC.2i_RPKM_TPM.txt")

PIM_RPKM_TPM.txt 

x<-list()
for (i in 1:length(xnames)){
   # x[[i]] <- merge(xi[[i]],normCounts, by.x='linked.ensemblid', by.y='row.names',all.x=T,all.y=T)
    x[[i]] <- merge(xi[[i]],RPKM_Counts, by.x='linked.ensemblid', by.y='Gene_ID',all.x=T,all.y=T)
    colnames(x[[i]])[length(colnames(x[[i]]))] <- 'gene.expression.log2CPM'
}
ESC <- x

# find GO terms of genes closest to highest 2000 H3K18lactylated enhancers
highest.2000.enhancer.H3K18la <-ESC[[1]][order(ESC[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

library(org.Mm.eg.db)
x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.ESC <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")

# Repeat for PIM
# Read bed files containing enhancer linked to closest 1 or 5 genes, overlapped with relevant peaks
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/")
xi <-list()
xnames <- list.files(path="/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/expression/",pattern='PIM')[1:4]
for (i in 1:length(xnames)){
    xi[[i]] <- read.delim(xnames[i],header=F)
    colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
     'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
     'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
    xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
}
names(xi) <-xnames


# read rawPeakCounts and normalise
rawPeakCounts.H3K18la <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/PIM_H3K18la_rawPeakCounts.rds")
rawPeakCounts.H3K4me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/PIM_H3K4me3_rawPeakCounts.rds")
rawPeakCounts.H3K27me3 <-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/PIM_H3K27me3_rawPeakCounts.rds")
normPeakCounts.H3K18la <- edgeR::cpm((rawPeakCounts.H3K18la[,1]+rawPeakCounts.H3K18la[,2])/2, log=TRUE)
normPeakCounts.H3K4me3 <- edgeR::cpm((rawPeakCounts.H3K4me3[,1]+rawPeakCounts.H3K4me3[,2])/2, log=TRUE)
normPeakCounts.H3K27me3 <- edgeR::cpm((rawPeakCounts.H3K27me3[,1]+rawPeakCounts.H3K27me3[,2])/2, log=TRUE)

normPeakCounts <- list(rowMeans(normPeakCounts.H3K18la),rowMeans(normPeakCounts.H3K27me3),rowMeans(normPeakCounts.H3K4me3))
names(normPeakCounts) <- c('normPeakCounts.H3K18la','normPeakCounts.H3K27me3','normPeakCounts.H3K4me3')
for (i in 1:3){
    xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
    colnames(xi[[i]])[length(colnames(xi[[i]]))] <- 'peak.log2CPM'
}

# Add gene expression info
RPKM_Counts <-read.delim("/home/egalle/public/_Projects/MS_H3K18la/Mouse/RNA/PIM_RPKM_TPM.txt")

x<-list()
for (i in 1:length(xnames)){
    x[[i]] <- merge(xi[[i]],RPKM_Counts, by.x='linked.ensemblid', by.y='Gene_ID',all.x=T,all.y=T)
    colnames(x[[i]])[length(colnames(x[[i]]))] <- 'gene.expression.log2CPM'
}


PIM<-x

# find GO terms of genes closest to highest 2000 enhancers per hPTM
highest.2000.enhancer.H3K18la <-PIM[[1]][order(PIM[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.PIM <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")
                   
                   
for (i in 1:length(GAS)){
	GAS[[i]]$gene.expression.RPKM <- rowMeans(cbind(GAS[[i]]$GAS_Con_1_RPKM,GAS[[i]]$GAS_Con_2_RPKM,GAS[[i]]$GAS_Con_3_RPKM))
		ESC[[i]]$gene.expression.RPKM <- rowMeans(cbind(ESC[[i]]$E14_SL_1_RPKM,ESC[[i]]$E14_SL_2_RPKM,ESC[[i]]$E14_SL_3_RPKM))
			PIM[[i]]$gene.expression.RPKM <- rowMeans(cbind(PIM[[i]]$pfkfb3.3d..HLI..macrophage.fastq_RPKM,PIM[[i]]$pfkfb3.3d..HLI..macrophage1.fastq_RPKM,PIM[[i]]$pfkfb3.3d..HLI..macrophage2.fastq_RPKM))
}                   

library(ggplot2)
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/SFig4C.cor.H3K18la.enhancer.vs.RNA.new.pdf")

    ggplot(ESC[[1]],aes(y=peak.log2CPM,x=log2(gene.expression.RPKM))) + geom_point(alpha = 0.3, color='blue') + ggtitle('ESC')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 RPKM)") + theme_bw() +  xlim(c(-10,15)) + ylim(c(2.5,12))
     ggplot(GAS[[1]],aes(y=peak.log2CPM,x=log2(gene.expression.RPKM))) + geom_point(alpha = 0.3, color='blue') + ggtitle('GAS')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 RPKM)")+ theme_bw() +  xlim(c(-10,15))#+ ylim(c(2.5,12))
     ggplot(PIM[[1]],aes(y=peak.log2CPM,x=log2(gene.expression.RPKM))) + geom_point(alpha = 0.3, color='blue') + ggtitle('PIM')+ylab("H3K18la promoter peak (log2 CPM)") + xlab("Gene expression (log2 RPKM)")+ theme_bw() + xlim(c(-10,15)) #+ ylim(c(2.5,12))

dev.off()

# Calculate Rho values for these correlations: 
Rho.gene.expression.enhancer.hPTM <- c(cor.test(ESC[[1]]$peak.log2CPM,ESC[[1]]$gene.expression.RPKM,method='spearman')[[4]],
        cor.test(GAS[[1]]$peak.log2CPM,GAS[[1]]$gene.expression.RPKM,method='spearman')[[4]],
        cor.test(PIM[[1]]$peak.log2CPM,PIM[[1]]$gene.expression.RPKM,method='spearman')[[4]])
names(Rho.gene.expression.enhancer.hPTM) <- c('ESC','GAS','PIM')


# Plot GO enrichments:
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig2F.GO enhancers top PTMs.new.pdf",width=7,height=5)
dotplot(GO.enhancer.H3K18la.ESC,showCategory=10,title='top ESC H3K18la enhancers')
dotplot(GO.enhancer.H3K18la.GAS,showCategory=10,title='top GAS H3K18la enhancers')
dotplot(GO.enhancer.H3K18la.PIM,showCategory=10,title='top PIM H3K18la enhancers')
dev.off()
       


# Suppl. Fig 4A

# in bash: create file containing all overlaps between peaks and public peak sets:

# Define master peak with all peaks (of unique data; only merged data sets)
QUERYbed=/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.bed.sorted.bed

# Define all public peak sets
refbed1=BMDM_CTCF.sorted.bed
refbed2=BMDM_H3K27ac.sorted.bed
refbed3=BMDM_H3K4me1.sorted.bed
refbed4=BMDM_H3K4me3.sorted.bed
refbed5=BMDM_PolIIa.sorted.bed
refbed6=ESC_H3K27ac.sorted.bed
refbed7=ESC_H3K27ac_Yan.sorted.bed
refbed8=ESC_H3K27me3_Perino.sorted.bed
refbed9=ESC_H3K36me3.sorted.bed
refbed10=ESC_H3K4me1.sorted.bed
refbed11=ESC_H3K4me1_Yan.sorted.bed
refbed12=ESC_H3K4me3_Perino.sorted.bed
refbed13=ESC_H3K4me3.sorted.bed
refbed14=ESC_H3K9me3.sorted.bed
refbed15=GAS_ATAC.sorted.bed
refbed16=GAS_CTCF.sorted.bed
refbed17=GAS_H3K27ac.sorted.bed
refbed18=GAS_H3K4me1.sorted.bed
refbed19=GAS_H3K4me3.sorted.bed
refbed20=MB_H3K18ac.bed.sorted.bed
refbed21=MB_H3K27me3.bed.sorted.bed
refbed22=MB_H3K36me3.bed.sorted.bed
refbed23=MB_H3K4me1.bed.sorted.bed
refbed24=MB_H3K4me2.bed.sorted.bed
refbed25=MB_H3K4me3.bed.sorted.bed
refbed26=MB_H3K9ac.bed.sorted.bed
refbed27=MB_H4K12ac.bed.sorted.bed
refbed28=MB_PolII.bed.sorted.bed
refbed29=MT_H3K18ac.bed.sorted.bed
refbed30=MT_H3K27me3.bed.sorted.bed
refbed31=MT_H3K36me3.bed.sorted.bed
refbed32=MT_H3K4me1.bed.sorted.bed
refbed33=MT_H3K4me2.bed.sorted.bed
refbed34=MT_H3K4me3.bed.sorted.bed
refbed35=MT_H3K9ac.bed.sorted.bed
refbed36=MT_H4K12ac.bed.sorted.bed
refbed37=MT_PolII.bed.sorted.bed
refbed38=BMDM_H3K18ac_M0.sorted.bed
refbed39=BMDM_H3K18ac_M1.sorted.bed
refbed40=BMDM_H3K18la_M0.sorted.bed
refbed41=BMDM_H3K18la_M1.sorted.bed
refbed42=ESC.2i_H3K4me1.sorted.bed
refbed43=ESC.2i_H3K4me3.sorted.bed
refbed44=ESC.2i_H3K27ac.sorted.bed
refbed45=ESC.2i_H3K27me3.sorted.bed

#!/bin/bash
cd /home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/used_public_peaks/

FILES=*bed

for f in $FILES
do
    cut -f1-3 < "$f" > "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/used_public_peaks/sorted/trimmed/$f"
done

# Go to folder where all public peaks are found (alternatively, add this to their name)
cd /home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/used_public_peaks/sorted/trimmed

# Create intersect file


bedtools intersect -wao -a ${QUERYbed} -b $refbed1 $refbed11 $refbed7 $refbed2 $refbed3 $refbed4 $refbed5 $refbed6 $refbed8 $refbed9 $refbed10  $refbed12 $refbed13 $refbed14 $refbed15 $refbed16 $refbed17 $refbed18 $refbed19 $refbed20 $refbed21 $refbed22 $refbed23 $refbed24 $refbed25 $refbed26 $refbed27 $refbed28 $refbed29 $refbed30 $refbed31 $refbed32 $refbed33 $refbed34 $refbed35 $refbed36 $refbed37 $refbed38 $refbed39 $refbed40 $refbed41 $refbed42 $refbed43 $refbed44 $refbed45 -sorted -filenames > ${QUERYbed}.overlap_w_all_publicPeaks.bed

# Analyse in R
# Read intersect file
R
library(dplyr)
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/overlap")
# overlap_publicpeaks<-read.delim("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/masterPeak/masterpeak.bed.sorted.bed.overlap_w_all_publicPeaks.bed" ,header=F)

# colnames(overlap_publicpeaks) <- c('chr.peak',
                                           # 'start.peak',
                                           # 'end.peak',
                                           # 'total.peak.value',
                                           # 'max.peak.value',
                                           # 'max.peak.value.region',
                                           # 'length.peak',
                                           # 'type.peak',
                                           # 'type.publicpeak',
                                           # 'chr.publicpeak',
                                           # 'start.publicpeak',
                                           # 'end.publicpeak',
                                            # 'bp.overlap')
# saveRDS(overlap_publicpeaks, "/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/overlap.RDS")

overlap_publicpeaks<-readRDS("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/overlap.RDS")
                                       
# Add the width of the the input files 
overlap_publicpeaks[,14] <-  overlap_publicpeaks[,13]
overlap_publicpeaks[,13] <-  overlap_publicpeaks[,12]-overlap_publicpeaks[,11]+1
colnames(overlap_publicpeaks) [13:14] <- c('length.publicpeak','bp.overlap')
# Calculate the overlap between any public peak set and any of our CUT&Tag peak sets
total.overlap<-data.frame(overlap_publicpeaks %>% group_by(type.peak,type.publicpeak) %>% summarise(sum=sum(as.numeric(paste0(bp.overlap)))))

# remove rows that contain overlap with no ref peak
total.overlap.f<-total.overlap[!total.overlap[,2]=='.',]

# read total peak bp in  merged peaks and adjust for merging with the total overlap file
total.peak.bp<-read.table("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/total.peak.bp.txt",sep="\t")
total.peak.bp <- cbind(total.peak.bp,rownames(total.peak.bp))
colnames(total.peak.bp) <- c("total.bp",'type.peak')

#merge for normalisation later on
total.overlap.f<-merge(total.overlap.f,total.peak.bp,by='type.peak')
total.overlap.f$type.publicpeak<-gsub(".sorted.bed","",total.overlap.f$type.publicpeak )
total.overlap.f$type.publicpeak<-gsub(".bed","",total.overlap.f$type.publicpeak )
# Calculate total public peak bp and save for speed
setwd("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/used_public_peaks/sorted/trimmed/")
# public.peaks.files <- list.files(pattern='.bed')
# public.peaks <- list()
# for (i in 1:length(public.peaks.files)){
  # public.peaks[[i]] <- read.table(paste0("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/public_peaks/used_public_peaks/sorted/trimmed/",
                                  # public.peaks.files[i]))
# }
# names(public.peaks)<-gsub(".sorted.bed","",public.peaks.files)
# names(public.peaks)<-gsub(".bed","",names(public.peaks))
# total.public.peak.bp<-list()
# for (i in 1:length(public.peaks)){
    # public.peaks[[i]]$width <- public.peaks[[i]][,3]-public.peaks[[i]][,2]
    # total.public.peak.bp[[i]]<-sum(public.peaks[[i]]$width )
# }
# total.public.peak.bp<-unlist(total.public.peak.bp)
# names(total.public.peak.bp)<-names(public.peaks)
# write.table(total.public.peak.bp,"/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/total.public.peak.bp.txt",sep="\t")
total.public.peak.bp<-read.table("/home/egalle/public/_Projects/MS_H3K18la/Eva_analysis/total.public.peak.bp.txt",sep="\t")
total.public.peak.bp <- cbind(total.public.peak.bp,rownames(total.public.peak.bp))
colnames(total.public.peak.bp) <- c("total.public.peak.bp",'type.publicpeak')
unique(total.overlap.f$type.publicpeak) %in% unique(total.public.peak.bp$type.publicpeak)
total.overlap.f<-merge(total.overlap.f,total.public.peak.bp,by='type.publicpeak')

# normalise total overlap for total size of our C&T peaks and for public peak sizes to obtain enrichment
total.overlap.f$normalized.sum <- total.overlap.f$sum/total.overlap.f$total.bp/total.overlap.f$total.public.peak.bp

# save all overlap information
total.overlap.temp<-total.overlap.f[,c(1,2,6)]

# reshape column normalized.sum (containing enrichment) to matrix with CUT&Tag peaks and public peaks in rows and columns resp.
library(reshape2)
total.overlap.ff <- dcast(total.overlap.f,type.peak~type.publicpeak, value.var='normalized.sum')
rownames(total.overlap.ff)<-total.overlap.ff[,1]
total.overlap.ff<-total.overlap.ff[,2:dim(total.overlap.ff)[2]]

library(pheatmap)
library(ggplot2)

breaksList = seq(-2, 2, by = 0.1)

total.overlap.subset<-total.overlap.ff[,-grep('CTCF|ATAC',colnames(total.overlap.ff))]
total.overlap.subset<-total.overlap.subset[,-grep('Yan|H3K4me3_Perino',colnames(total.overlap.subset))]
total.overlap.subset<-total.overlap.subset[,-grep('PolII|H3K9ac',colnames(total.overlap.subset))]
total.overlap.subset<-total.overlap.subset[,-grep('H3K4me2|H4K12ac',colnames(total.overlap.subset))]
total.overlap.subset<-total.overlap.subset[,-grep('H3K36|H3K9',colnames(total.overlap.subset))]

 total.overlap.subset$BMDM_H3K18la <- (total.overlap.subset$BMDM_H3K18la_M0+total.overlap.subset$BMDM_H3K18la_M1)/2
 total.overlap.subset$BMDM_H3K18ac <- (total.overlap.subset$BMDM_H3K18ac_M0+total.overlap.subset$BMDM_H3K18ac_M1)/2
 total.overlap.subset<-total.overlap.subset[,-grep('_M1|_M0',colnames(total.overlap.subset))]

total.overlap.fff<-total.overlap.subset
rownames(total.overlap.subset) <- gsub("_merged","",rownames(total.overlap.subset))
total.overlap.ffff<-total.overlap.subset[rownames(total.overlap.subset) %in% c("ADI_H3K18la","ADI_H3K27ac","BMDM_H3K18la_1","BMDM_H3K27ac_1","BMDM_H3K27me3_1","ESC.2i_H3K18la","ESC.2i_H3K27me3","ESC.2i_H3K4me3","ESC.ser_H3K18la","ESC.ser_H3K27ac_1","ESC.ser_H3K27me3","ESC.ser_H3K4me3","GAS_H3K18la","GAS_H3K27ac","GAS_H3K27me3","GAS_H3K4me3","MB_H3K18la_1","MB_H3K27me3_1","MT_H3K18la","MT_H3K27me3","PIM_H3K18la","PIM_H3K27ac_1","PIM_H3K27me3","PIM_H3K4me3"),]

library(RColorBrewer)
# plot enrichment
pdf("/home/egalle/public/_Projects/MS_H3K18la/Revision/GitHub/try_output/SFig4A.heatmap.overlap.peaks.pdf",height=15)
pheatmap(as.matrix(total.overlap.ffff),
 cluster_rows=T,cluster_cols=T,legend=T,scale='row',cellheight=10,main='overlap our C&T peaks vs public ChIPseq peaks',breaks = breaksList,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)))
dev.off()











