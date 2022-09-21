
# This R script has been used to generate the following figures:
### Figures 2D, 4H

### Compute overlaps between CnT peaks and cCREs in BASH

# Define input
# resultsPath=../CnT_cCRE_overlaps
# refbed=../mm10-ccREs.bed.sorted.bed
# peakPath=../peaks/

# Calculate overlaps
# cd $peakPath
# for i in *bam*
#   do
  # awk '{print "chr"$0}' $peakPath2/$i > $peakPath/$i.w.chr.bed  # add 'chr' to chromosome column
  # sort -k1,1 -k2,2n $peakPath/$i.w.chr.bed > $peakPath/$i.w.chr.sorted.bed  # sort bed files to speed up the computation
  # cd $resultsPath 
  # bedtools intersect -wao -a $peakPath/$i.w.chr.sorted.bed -b refbed -sorted > $i.overlap.cCREs.bed
# done

######################################################################################################################################

library(GenomicRanges)
library(ggVennDiagram)
library(gplots)
library(gridExtra)							
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Read cCRE genomic coordinates
cCREs <- read.table("../mm10-ccREs.bed.sorted.bed")
colnames(cCREs) <- c('chr.cCRE','start.cCRE','end.cCRE','id1.cCRE','id2.cCRE','type.cCRE')

# Load peak and cCRE overlapping peak info
peaks.files <- list.files(pattern='bed.w.chr.sorted.bed')
peaks.overlap.cCREs.files <- list.files(pattern='overlap.cCREs.bed')

peaks.overlap.cCRE <- list()
peaks <- list()

for (i in 1:length(peaks.files)){
  peaks.overlap.cCRE[[i]] <- read.delim(paste0("../CnT_cCRE_overlaps/",peaks.overlap.cCREs.files[i]),header=F)
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
  peaks[[i]] <- read.table(paste0("../peaks/",peaks.files[i]))
  colnames(peaks[[i]]) <-colnames(peaks.overlap.cCRE[[i]])[1:6]
}

# Rename peaks and overlap-files to remove superfluous name parts
names(peaks)<-gsub("_0.01.stringent.bed","",peaks.files)
names(peaks)<-gsub(".w.chr.sorted.bed","",names(peaks))
names(peaks)<-gsub(".0.01.stringent.bed","",names(peaks))

names(peaks.overlap.cCRE)<-gsub("_0.01.stringent.bed.overlap.cCREs.bed","",peaks.overlap.cCREs.files)
names(peaks.overlap.cCRE)<-gsub("_0.01.stringent.bed.w.chr.sorted.bed.overlap.cCREs.bed","",names(peaks.overlap.cCRE))
names(peaks.overlap.cCRE)<-gsub(".0.01.stringent.bed.overlap.cCREs.bed","",names(peaks.overlap.cCRE))

# list which peaks are overlapping with the pELS and dELS per tissue
x.dELS <-list()
x.pELS <-list()

for (i in c(1:length(peaks.overlap.cCRE))){
  peaks.overlap.cCRE[[i]]$enhancer.id <- paste0(peaks.overlap.cCRE[[i]]$chr.cCRE,":",
                                                peaks.overlap.cCRE[[i]]$start.cCRE,"-",
                                                peaks.overlap.cCRE[[i]]$end.cCRE)
  x.dELS[[i]] <- peaks.overlap.cCRE[[i]]$enhancer.id[peaks.overlap.cCRE[[i]]$type.cCRE=='dELS']
  x.pELS[[i]] <- peaks.overlap.cCRE[[i]]$enhancer.id[peaks.overlap.cCRE[[i]]$type.cCRE=='pELS']
}

names(x.dELS) <- names(peaks.overlap.cCRE)
names(x.pELS) <- names(peaks.overlap.cCRE)

# Plot Venn diagram showing the pELS/dELS overlaps
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

ggvennlist.cCRE[[3]]<-ggVennDiagram(x.dELS[names(x.dELS) %in% c('GAS_H3K18la_merged','GAS_H3K4me3_merged',"GAS_H3K27ac_merged")],
                                    label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('dELS GAS')+ scale_color_manual(values=c('blue','orange','darkgreen'))

ggvennlist.cCRE[[4]]<-ggVennDiagram(x.pELS[names(x.pELS) %in% c('GAS_H3K18la_merged','GAS_H3K4me3_merged',"GAS_H3K27ac_merged")],
                                    label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('pELS GAS')+ scale_color_manual(values=c('blue','orange','darkgreen'))

ggvennlist.cCRE[[5]]<-ggVennDiagram(x.dELS[names(x.dELS) %in% c('ESC.ser_H3K18la_merged','ESC.ser_H3K4me3_merged',"ESC.ser_H3K27ac_1")],
                                    label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('dELS ESC.ser')+ scale_color_manual(values=c('blue','orange','darkgreen'))

ggvennlist.cCRE[[6]]<-ggVennDiagram(x.pELS[names(x.pELS) %in% c('ESC.ser_H3K18la_merged','ESC.ser_H3K4me3_merged',"ESC.ser_H3K27ac_1")],
                                    label_alpha = 0,set_size = 3,label_size = 3,set_color = c(1:3))+
  ggplot2::scale_fill_gradient(low="white",high = 'pink')+ 
  theme(text = element_text(size = 7),plot.title = element_text(size=10))+
  ggtitle('pELS ESC.ser')+ scale_color_manual(values=c('blue','orange','darkgreen'))

do.call(grid.arrange,c(ggvennlist.cCRE[1:6],ncol=3))
