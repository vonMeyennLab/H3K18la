# This R script has been used to generate the following figures:
### Figures 2F, 4K
### Supp. Figure S7B

library(clusterProfiler)
library(GenomicRanges)
library(gridExtra)
library(org.Mm.eg.db)
library(tidyr)

# Identify closest non-overlapping (n) gene promoter(s) to each dELS (done in BASH) using bedtools
## bedtools closest -io -k 1 -a ${dELS.bed} -b ${mm10_promoters} -D ref > dELS.<tissue>.closest.non.overlapping.promoters.bed

# Calculate overlaps between CnT peaks and dELS linked to closest genes (done in BASH) using bedtools for each tissue
## bedtools intersect -wao -a dELS.<tissue>.closest.non.overlapping.promoters.bed -b CnT_peak.bed > dELS.${n}.closest.non.overlapping.promoters.overlap

# Read bed files containing enhancers-linked-to-genes overlapped with peak files for GAS samples
setwd("$WORKDIR")
xi <-list()
xnames <- list.files(path="$WORKDIR",pattern='GAS')[1:4]

for (i in 1:length(xnames))
{
  xi[[i]] <- read.delim(xnames[i],header=F)
  
  colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
                         'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
                         'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
  xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
  
}
names(xi) <-xnames

# Load normalized gene counts from RNAseq and promoter based CnT data for each histone mark
# Refer to R script 5B_CnT_RNA_correlation_inELS.R to check how this count table was generated
normPeakCounts.H3K18la <-readRDS("../GAS_H3K18la_normPeakCounts.rds")
normPeakCounts.H3K4me3 <-readRDS("../GAS_H3K4me3_normPeakCounts.rds")
normPeakCounts.H3K27me3 <-readRDS("../GAS_H3K27me3_normPeakCounts.rds")

normPeakCounts <- list(normPeakCounts.H3K18la,normPeakCounts.H3K18la,normPeakCounts.H3K27me3,normPeakCounts.H3K4me3)

for (i in 1:length(xnames)){
  xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
  colnames(xi[[i]])[colnames(xi[[i]])=='V1'] <- 'peak.log2RPKM'
}

GAS <- xi

# Perform GO enrichment analysis of 2000 genes closest to dELS with highest normalized RPKM for each histone mark
highest.2000.enhancer.H3K18la <-GAS[[1]][order(GAS[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.GAS <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")


# Read bed files containing enhancers-linked-to-genes overlapped with peak files for ESC samples
setwd("$WORKDIR")
xi <-list()
xnames <- list.files(path="$WORKDIR",pattern='ESC')[1:4]

for (i in 1:length(xnames))
{
  xi[[i]] <- read.delim(xnames[i],header=F)
  colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
                         'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
                         'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
  xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
  
}
names(xi) <-xnames

# Load normalized gene counts from RNAseq and promoter based CnT data for each histone mark
# Refer to R script 5B_CnT_RNA_correlation_inELS.R to check how this count table was generated
normPeakCounts.H3K18la <-readRDS("../ESC_H3K18la_normPeakCounts.rds")
normPeakCounts.H3K4me3 <-readRDS("../ESC_H3K4me3_normPeakCounts.rds")
normPeakCounts.H3K27me3 <-readRDS("../ESC_H3K27me3_normPeakCounts.rds")

normPeakCounts <- list(normPeakCounts.H3K18la,normPeakCounts.H3K18la,normPeakCounts.H3K27me3,normPeakCounts.H3K4me3)

for (i in 1:length(xnames)){
  xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
  colnames(xi[[i]])[colnames(xi[[i]])=='V1'] <- 'peak.log2CPM'
}

ESC <- xi

# Perform GO enrichment analysis of 2000 genes closest to dELS with highest normalized RPKM for each histone mark
highest.2000.enhancer.H3K18la <-ESC[[1]][order(ESC[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.ESC <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")

# Read bed files containing enhancers-linked-to-genes overlapped with peak files for PIM samples
setwd("$WORKDIR")
xi <-list()
xnames <- list.files(path="$WORKDIR",pattern='PIM')[1:4]

for (i in 1:length(xnames))
{
  xi[[i]] <- read.delim(xnames[i],header=F)
  colnames(xi[[i]]) <- c('enhancer.chr','enhancer.start','enhancer.end','enhancer.id1','enhancer.id2','enhancer.type','linked.promoter.chr','linked.promoter.start','linked.promoter.end',
                         'linked.gene.length','linked.gene.strand','linked.ensemblid','linked.gene.id','linked.gene.type','distance.enhancer.to.promoter',
                         'peak.chr','start.peak','end.peak','total.peak.value','max.peak.value','max.peak.region','overlap.peak.enhancer.bp')
  xi[[i]]$peak.id <- paste0(gsub('chr','',xi[[i]]$peak.chr),':',xi[[i]]$start.peak, '-',xi[[i]]$end.peak)
}
names(xi) <-xnames

# Load normalized gene counts from RNAseq and promoter based CnT data for each histone mark
# Refer to R script 5B_CnT_RNA_correlation_inELS.R to check how this count table was generated
normPeakCounts.H3K18la <-readRDS("../PIM_H3K18la_normPeakCounts.rds")
normPeakCounts.H3K4me3 <-readRDS("../PIM_H3K4me3_normPeakCounts.rds")
normPeakCounts.H3K27me3 <-readRDS("../PIM_H3K27me3_normPeakCounts.rds")

normPeakCounts <- list(rowMeans(normPeakCounts.H3K18la),rowMeans(normPeakCounts.H3K27me3),rowMeans(normPeakCounts.H3K4me3))
names(normPeakCounts) <- c('normPeakCounts.H3K18la','normPeakCounts.H3K27me3','normPeakCounts.H3K4me3')

for (i in 1:3){
  xi[[i]] <- merge(xi[[i]], normPeakCounts[[i]], by.x='peak.id',by.y='row.names', all.x=T,all.y=T)
  colnames(xi[[i]])[length(colnames(xi[[i]]))] <- 'peak.log2CPM'
}

PIM <- xi

# Perform GO enrichment analysis of 2000 genes closest to dELS with highest normalized RPKM for each histone mark
highest.2000.enhancer.H3K18la <-PIM[[1]][order(PIM[[1]]$peak.log2CPM,decreasing=T),][1:2000,]

x<-mapIds(org.Mm.eg.db, highest.2000.enhancer.H3K18la$linked.ensemblid, 'ENTREZID', 'ENSEMBL')
GO.enhancer.H3K18la.PIM <- enrichGO(unique(as.numeric(unlist(x))),pvalueCutoff  = 0.05,
                                    pAdjustMethod = "BH",OrgDb='org.Mm.eg.db', ont="ALL")
# Plot GO enrichments
dotplot(GO.enhancer.H3K18la.ESC,showCategory=10,title='top ESC H3K18la enhancers')
dotplot(GO.enhancer.H3K18la.GAS,showCategory=10,title='top GAS H3K18la enhancers')
dotplot(GO.enhancer.H3K18la.PIM,showCategory=10,title='top PIM H3K18la enhancers')
