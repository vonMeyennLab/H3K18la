# This R script has been used to generate the following figures:
### Figures 2D, 4H

library(GenomicRanges)
library(ggVennDiagram)
library(gplots)
library(gridExtra)		
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Read cCRE genomic coordinates
cCREs <- read.table("../mm10-ccREs.bed.sorted.bed")
colnames(cCREs) <- c('chr.cCRE','start.cCRE','end.cCRE','id1.cCRE','id2.cCRE','type.cCRE')

# Calculate overlaps between CnT peaks and cCREs (done in BASH) using bedtools
## bedtools intersect -wao -a CnT_peak.bed -b cCRE.bed -sorted > <tissue>.overlap.cCREs.bed

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

######################################################################################################################################

### Figures 2C, 4J

# Create a master peak file for calculating overlaps with tissue-specific enhancers
masterPeak_df <- bind_rows(peaks, .id = "column_label")
masterPeak_df <- masterPeak_df[,c(2,3,4,5,6,7,8,1)]

# Calculate overlaps between CnT peaks and tissue-specific enhancers (done in BASH) using bedtools
# For ESC, ENCODE cCRE dataset has been used which was further filtered for ELS (pELS + dELS)
## bedtools intersect -wao -a $tissueEnhancers.bed -b masterPeak.bed > <tissue>.enhancers.covered.by.peaks.bed

# Calculate fraction of gastroc enhancers that is covered by each CnT peak set
gastroc.enhancers.covered.by.peaks <- read.delim("../gastroc.enhancers.covered.by.peaks.bed",header=F)
colnames(gastroc.enhancers.covered.by.peaks) <- c('chr.enhancer','start.enhancer','end.enhancer','enhancer.nr','peak.chr','peak.start','peak.end',
                                                  'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
gastroc.enhancers.covered <- c()
for (i in 1:length(unique(gastroc.enhancers.covered.by.peaks$peak.origin)))
{
  peakOI<-unique(gastroc.enhancers.covered.by.peaks$peak.origin)[i]
  gastroc.enhancers.covered[i]<-length(unique(gastroc.enhancers.covered.by.peaks$enhancer.nr[gastroc.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(gastroc.enhancers.covered.by.peaks$enhancer.nr))
  
}

names(gastroc.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(gastroc.enhancers.covered.by.peaks$peak.origin))
names(gastroc.enhancers.covered)[names(gastroc.enhancers.covered)=='.']<-'Not covered by any peak'

# Calculate fraction of adipose tissue enhancers that is covered by each CnT peak set
ADI.enhancers.covered.by.peaks <- read.delim("../ADI.enhancers.covered.by.peaks.bed",header=F)
colnames(ADI.enhancers.covered.by.peaks) <- c('chr.enhancer','start.enhancer','end.enhancer','enhancer.nr','peak.chr','peak.start','peak.end',
                                              'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
ADI.enhancers.covered <- c()
for (i in 1:length(unique(ADI.enhancers.covered.by.peaks$peak.origin)))
{
  peakOI<-unique(ADI.enhancers.covered.by.peaks$peak.origin)[i]
  ADI.enhancers.covered[i]<-length(unique(ADI.enhancers.covered.by.peaks$enhancer.nr[ADI.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(ADI.enhancers.covered.by.peaks$enhancer.nr))
  
}

names(ADI.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(ADI.enhancers.covered.by.peaks$peak.origin))
names(ADI.enhancers.covered)[names(ADI.enhancers.covered)=='.']<-'Not covered by any peak'

# Calculate fraction of MB/MT enhancers that is covered by each CnT peak set
MB.enhancers.covered.by.peaks <- read.delim("../MB.enhancers.covered.by.peaks.bed",header=F)
colnames(MB.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.id','peak.chr','peak.start','peak.end',
                                             'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')
MT.enhancers.covered.by.peaks <- read.delim("../MT.enhancers.covered.by.peaks.bed",header=F)
colnames(MT.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.id','peak.chr','peak.start','peak.end',
                                             'peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')

MT.enhancers.covered <- c()
MB.enhancers.covered <- c()

for (i in 1:length(unique(MB.enhancers.covered.by.peaks$peak.origin)))
{
  peakOI<-unique(MB.enhancers.covered.by.peaks$peak.origin)[i]
  MB.enhancers.covered[i]<-length(unique(MB.enhancers.covered.by.peaks$enhancer.id[MB.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(MB.enhancers.covered.by.peaks$enhancer.id))
  peakOI<-unique(MT.enhancers.covered.by.peaks$peak.origin)[i]
  MT.enhancers.covered[i]<-length(unique(MT.enhancers.covered.by.peaks$enhancer.id[MT.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(MT.enhancers.covered.by.peaks$enhancer.id))
  
}

names(MB.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(MB.enhancers.covered.by.peaks$peak.origin))
names(MT.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(MT.enhancers.covered.by.peaks$peak.origin))
names(MB.enhancers.covered)[names(MB.enhancers.covered)=='.']<-'Not covered by any peak'
names(MT.enhancers.covered)[names(MT.enhancers.covered)=='.']<-'Not covered by any peak'

# Calculate fraction of BMDM enhancers that is covered by each CnT peak set
BMDM.Denisenko.enhancers.covered.by.peaks <- read.delim("../BMDM.Denisenko.enhancers.covered.by.peaks.bed",header=F)
colnames(BMDM.Denisenko.enhancers.covered.by.peaks) <- c("chr.enhancer",'start.enhancer','end.enhancer','enhancer.target.chr','enhancer.target.start','enhancer.target.end','enhancer.strand','target.gene','chr.peak','start.peak','end.peak','peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')

BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr <- paste0(BMDM.Denisenko.enhancers.covered.by.peaks$chr.enhancer,BMDM.Denisenko.enhancers.covered.by.peaks$start.enhancer,BMDM.Denisenko.enhancers.covered.by.peaks$end.enhancer)
BMDM.Denisenko.enhancers.covered <- c()

for (i in 1:length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin)))
{
  peakOI<-unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin)[i]
  BMDM.Denisenko.enhancers.covered[i]<-length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr[BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(BMDM.Denisenko.enhancers.covered.by.peaks$enhancer.nr))
  
}
names(BMDM.Denisenko.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin))
BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin[BMDM.Denisenko.enhancers.covered.by.peaks$peak.origin=='.'] <- 'Not covered by any peak'
names(BMDM.Denisenko.enhancers.covered)[names(BMDM.Denisenko.enhancers.covered)=='.'] <- 'Not covered by any peak'

# Calculate fraction of ESC enhancers that is covered by each CnT peak set
E14.cCREs.covered.by.peaks <- read.delim("../E14.cCREs.covered.by.peaks.bed",header=F)
colnames(E14.cCREs.covered.by.peaks) <- c("chr.cCRE",'start.cCRE','end.cCRE','cCRE.id',rep('NA',5),'type.cCRE','NA',
                                          'chr.peak','start.peak','end.peak','peak.total.value','peak.max.value','peak.max.value.region','peak.width','peak.origin','overlap.bp')

E14.cCREs.covered.by.peaks<-E14.cCREs.covered.by.peaks[,-grep('NA',colnames(E14.cCREs.covered.by.peaks))]
E14.enhancers.covered.by.peaks<-E14.cCREs.covered.by.peaks[grep('ELS',E14.cCREs.covered.by.peaks$type.cCRE),]
E14.enhancers.covered <- c()

for (i in 1:length(unique(E14.enhancers.covered.by.peaks$peak.origin)))
{
  peakOI<-unique(E14.enhancers.covered.by.peaks$peak.origin)[i]
  E14.enhancers.covered[i]<-length(unique(E14.enhancers.covered.by.peaks$cCRE.id[E14.enhancers.covered.by.peaks$peak.origin==peakOI]))/length(unique(E14.enhancers.covered.by.peaks$cCRE.id))
}
names(E14.enhancers.covered) <- gsub(".w.chr.sorted.bed","",unique(E14.enhancers.covered.by.peaks$peak.origin))
E14.enhancers.covered.by.peaks$peak.origin[E14.enhancers.covered.by.peaks$peak.origin=='.'] <- 'Not covered by any peak'

# Combine the tissue specific enhancer overlaps
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

temp <- data.frame(x=rownames(enhancers.covered.selection))
tissue <- temp %>% separate(x, c('Tissue','Histone','Rep'), sep ='_')

# For each tissue-specific enhancer, plot the fraction covered by each CnT peak set
par(mar=c(5,5,2,1),mfrow=c(1,6),cex=1)
for (i in 1:dim(enhancers.covered.selection)[2])
{
  barplot(enhancers.covered.selection[,i][order(enhancers.covered.selection[,i])],horiz=T,
          las=1,xlim=c(0,1),col=cols[,2][order(enhancers.covered.selection[,i])],
          names = tissue[,1][order(enhancers.covered.selection[,i])],
          xlab=paste0('Fraction of ',colnames(enhancers.covered.selection)[i],' enhancers covered'))
}
