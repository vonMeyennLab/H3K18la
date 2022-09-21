# This R script has been used to generate the following figures:
### Figures 1F, 4D
### Supp. Figures 3C, 6F

library(biomaRt)
library(ChIPseeker)
library(ChIPpeakAnno)
library(clusterProfiler)
library(GenomicRanges)
library(ggVennDiagram)
library(gplots)
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

######################################################################################################################################

### Supp. Figure 3A
                                 
names(genes.w.promoter.peaks) <- gsub('.w.chr.sorted.bed','', names(genes.w.promoter.peaks))

# Enumerate sets of all possible groups
venn_GAS<- venn(genes.w.promoter.peaks[grep('GAS',names(genes.w.promoter.peaks))], show.plot = T)
venn_PIM<- venn(genes.w.promoter.peaks[grep('PIM',names(genes.w.promoter.peaks))], show.plot = T)
venn_ESC.ser<- venn(genes.w.promoter.peaks[grep('ESC.ser',names(genes.w.promoter.peaks))], show.plot = T)

# Update sample names
names(attributes(venn_GAS)$intersection)<-gsub('_merged','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_GAS)$intersection)<-gsub('GAS_','',names(attributes(venn_GAS)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_merged','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('_1','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_PIM)$intersection)<-gsub('PIM_','',names(attributes(venn_PIM)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_merged','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_public.bed','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('_1','',names(attributes(venn_ESC.ser)$intersection))
names(attributes(venn_ESC.ser)$intersection)<-gsub('ESC.ser_','',names(attributes(venn_ESC.ser)$intersection))

# GO CC enrichment analysis for GAS samples
x <- list(unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_GAS)$intersection[names(attributes(venn_GAS)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')

compCC_GAS_sub <- compareCluster(geneCluster = x,
                                 fun           = "enrichGO",
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH",
                                 OrgDb='org.Mm.eg.db',
                                 ont="CC")

# GO CC enrichment analysis for PIM samples
x <- list(unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_PIM)$intersection[names(attributes(venn_PIM)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')


compCC_PIM_sub <- compareCluster(geneCluster = x,
                                 fun           = "enrichGO",
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH",
                                 OrgDb='org.Mm.eg.db',
                                 ont="CC")

# GO CC enrichment analysis for ESC.ser samples
x <- list(unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K4me3'])))
names(x) <- c('group1','group2','group3')

compCC_ESC.ser_sub <- compareCluster(geneCluster = x,
                                     fun           = "enrichGO",
                                     pvalueCutoff  = 0.05,
                                     pAdjustMethod = "BH",
                                     OrgDb='org.Mm.eg.db',
                                     ont="CC")

# Combine the dot plots showing the GO CC enrichment analyses
dp<-list()
dp[[1]] <- dotplot(compCC_GAS_sub, showCategory=15, title = "GAS CC GO Enrichment Analysis", includeAll=T)
dp[[2]] <- dotplot(compCC_PIM_sub, showCategory=15, title = "PIM CC GO Enrichment Analysis", includeAll=T)
dp[[3]] <- dotplot(compCC_ESC.ser_sub, showCategory=15, title = "ESC.ser CC GO Enrichment Analysis", includeAll=T)

do.call("grid.arrange",c(dp[c(3,1,2)],ncol=3))
                                 
######################################################################################################################################

### Figures 1H, 4F
                                 
# Load raw gene counts from RNAseq and promoter based CnT data
ESC.ser <- read.delim("../ESC_hPTM_RNA_prom_rawCounts.txt")

# Annotate gene IDs with gene names
geneAnnotation <- getBM(attributes=c("ensembl_gene_id",'entrezgene_id', "external_gene_name","chromosome_name","start_position", "end_position"), 
                        filters="ensembl_gene_id", 
                        values=ESC.ser$EnsemblID, mart=ensembl)
geneAnnotation$width = geneAnnotation$end_position-geneAnnotation$start_position   

# Calculate RPKM values
rpm <- merge(ESC.ser, geneAnnotation,
             by.y="ensembl_gene_id", 
             by.x="EnsemblID",all.y=T,all.x=T)
rpm$RNA<-(2^rpm$RNA)/rpm$width*1000

# Generate all possible groups of CnT histone marks
x <- list(unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac:H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K27ac'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la:H3K4me3'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K27ac'])),
          unique(unlist(attributes(venn_ESC.ser)$intersection[names(attributes(venn_ESC.ser)$intersection)=='H3K18la'])))
names(x) <- c('all','H3K27ac-H3K4me3','H3K4me3','H3K18la-H3K27ac','H3K18la-H3K4me3','H3K27ac','H3K18la')

rpm$Promoter <- NA
for (i in 1:length(names(x))){
  rpm$Promoter[rpm$GeneName%in%x[[i]]]<- names(x)[i]
}
rpm$Promoter[is.na(rpm$Promoter)] <-'no active mark'
rpm_ESC.ser<-rpm

# Define group colors
colors_to_use <- cbind(c('red','blue','lightblue','darkblue','orange','orange3','darkgreen','red'), sort(unique(rpm_ESC.ser$Promoter)))
colnames(colors_to_use) <- c('color','Promoter')

# Box plot showing the RPKM distribution of genes for each CnT histone group
new_order <- with(rpm_ESC.ser,reorder(Promoter, RNA, median, na.rm=T))
no2 <- data.frame(rpm_ESC.ser %>% group_by(Promoter) %>% summarise(median=median(RNA)))
no2 <- merge(no2, colors_to_use, by = "Promoter")
boxplot(log2(rpm_ESC.ser$RNA+1)~new_order,outline=F,
        col=no2$color[order(no2$median)], ylab = '', xlab = 'Gene expression log2(RPKM)',notch=T,
        horizontal=T,las=1,main='ESC-ser',ylim=c(0,7))

text(y=1:8, x=6, paste0('n=',table(new_order)),cex=0.7)


                              
