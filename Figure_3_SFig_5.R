# Fig 3A + SFig5A

library(ggpubr)
library(vioplot)

# load raw peak counts (union of MB and MT /ESC 2i and ser peaks) and normalise to RPM

MB_MT_H3K18la_rawPeakCounts<-readRDS("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/MB_MT_H3K18la_rawPeakCounts.rds")
libSize <- Matrix::colSums(MB_MT_H3K18la_rawPeakCounts)
MB_MT_H3K18la_rpmCounts <- sweep(MB_MT_H3K18la_rawPeakCounts * 1e6, MARGIN=2, STATS=libSize, FUN="/")

ESC_H3K18la_rawPeakCounts<-readRDS("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Mouse/hPTM_rawPeakCounts/ESC_H3K18la_rawPeakCounts.rds")
libSize <- Matrix::colSums(ESC_H3K18la_rawPeakCounts)
ESC_H3K18la_rpmCounts <- sweep(ESC_H3K18la_rawPeakCounts * 1e6, MARGIN=2, STATS=libSize, FUN="/")

# Calculate FC from MB to MT and from 2i to ser

MB_MT_H3K18la_rpmFC <- (1+MB_MT_H3K18la_rpmCounts)/(MB_MT_H3K18la_rpmCounts[,1]+1)
MB_MT_H3K18la_rpmFC.mean<-rowMeans(MB_MT_H3K18la_rpmFC[,2:3])

ESC_H3K18la_rpmFC <- ESC_H3K18la_rpmCounts/(1+rowMeans(ESC_H3K18la_rpmCounts[,3:4]))
ESC_H3K18la_rpmFC.mean<-rowMeans(ESC_H3K18la_rpmFC[,1:2])

# Get coordinates of union peaks
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

# Annotate peaks to genomic elements

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

# annotate peaks to genomic elements with ChIPseeker
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
promoter_mm <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=2000, downstream=2000)

peaks_anno_MB_MT <- annotatePeak(MB_MT_peaks.ggranges,
                                tssRegion=c(-2000, 2000),
                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                annoDb='org.Mm.eg.db')
peaks_anno_ESC <- annotatePeak(ESC_peaks.ggranges,
                                 tssRegion=c(-2000, 2000),
                                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                 annoDb='org.Mm.eg.db')
                                 
                                 
color <- peaks_anno_MB_MT@anno@elementMetadata$annotation
color[grep('Intron',color)] <- 'Intron'
color[grep('Exon',color)] <- 'Exon'
color_MB_MT <- color

color <- peaks_anno_ESC@anno@elementMetadata$annotation
color[grep('Intron',color)] <- 'Intron'
color[grep('Exon',color)] <- 'Exon'
color_ESC <- color
                                 
                                 
pdf("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig3A.peak_rpm_per_genomic_feature.pdf",width=10,height = 7)
par(mfrow=c(2,2),mar=c(6,10,5,5))
boxplot(abs(log2(rowMeans(MB_MT_H3K18la_rpmCounts[,2:3])/(1+MB_MT_H3K18la_rpmCounts[,1])))~color_MB_MT,
        horizontal=T,las=1, main = 'MB to MT',ylab='',xlab='abs(log2 FC)',notch=T,outpch=16,outcex=0.4)
abline(v=median(abs(log2(rowMeans(MB_MT_H3K18la_rpmCounts[,2:3])/(1+MB_MT_H3K18la_rpmCounts[,1])))[color_MB_MT=='Promoter (<=1kb)']),lty=2,col='red')
dev.off()

pdf("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Revision/GitHub/try_output/SFig5A.peak_rpm_per_genomic_feature.pdf",width=10,height = 7)
par(mfrow=c(2,2),mar=c(6,10,5,5))
boxplot(abs(log2(rowMeans(ESC_H3K18la_rpmCounts[,1:2])/(1+rowMeans(ESC_H3K18la_rpmCounts[,3:4]))))~color_ESC,
        horizontal=T,las=1, main = '2i to ser',ylab='',xlab='abs(log2 FC)',notch=T,outpch=16,outcex=0.4)
abline(v=median(abs(log2(rowMeans(ESC_H3K18la_rpmCounts[,3:4])/(1+ESC_H3K18la_rpmCounts[,1:2])))[color_ESC=='Promoter (<=1kb)']),lty=2,col='red')
dev.off()

# Fig 3B,D,E,F and SFig5B,D,E,F: Adhi
# Fig 3C and SFig5B : IGV screenshots

# Fig 3G

#### Now look at dynamic events in peaks localised in tissue-specific enhancers:
# First, intersect the enhancer bed files with the peak files:
# peak=/home/egalle/public/_Projects/MS_H3K18la/Mouse/MT_vs_MB/MB_MT_union_peaks_rpm.bed
# 
# sort -k1,1 -k2,2n $peak > $peak.sorted.bed
# 
# refbed5=/home/egalle/public/EvaGalle/Data/Public_Data/Blum_2012_GenesDev_myoblasts.myotubes/MB_enhancers.mm10.bed
# refbed6=/home/egalle/public/EvaGalle/Data/Public_Data/Blum_2012_GenesDev_myoblasts.myotubes/MT_enhancers.mm10.bed
# 
# bedtools intersect -wao -a $peak.sorted.bed -b $refbed5.sorted.bed $refbed6.sorted.bed -sorted -filenames > $peak.overlap.MB_MT.enhancers.bed
# 
setwd("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Mouse/MT_vs_MB/")
MB_MT.peaks.overlap.w.tissue.enhancers <- read.table(list.files(pattern="overlap.MB_MT.enhancers.bed"),header=F)
colnames(MB_MT.peaks.overlap.w.tissue.enhancers) <- c('chr.peak','start.peak','end.peak','MB_1','MT_1','MT_2','enhancer.origin','chr.enhancer','start.enhancer','end.enhancer','overlap.bp')

pdf("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig3G.MT vs MB H3K18la in myogenic enhancers.pdf",width=5,height = 7)
par(mai=c(1.5,0.5,0.5,0.5),mfrow=c(2,3),xpd=T)
boxplot(log2(rowMeans(MB_MT.peaks.overlap.w.tissue.enhancers[,5:6])/MB_MT.peaks.overlap.w.tissue.enhancers[,4])~paste0(c(MB_MT.peaks.overlap.w.tissue.enhancers$overlap.bp==0),MB_MT.peaks.overlap.w.tissue.enhancers$enhancer.origin),
        notch=T,outline=T,names = c('MB enh','MT enh','no enh'), xlab = '', ylab = 'log2FC',main='MT vs MB H3K18la in myogenic enhancers',outpch=16, col = c('lightblue', 'blue4','seagreen2'))
dev.off()


# Fig 3H

#Define genes with H3K18la peak in promoter in MB/MT samples
peaks.MB.H3K18la <- read.table("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks/MB_H3K18la_1_0.01.stringent.bed.w.chr.sorted.bed")
colnames(peaks.MB.H3K18la) <- c('chr','start','end')
peaks.MB.H3K18la_ggranges <- makeGRangesFromDataFrame(peaks.MB.H3K18la,
                                                      keep.extra.columns=TRUE,
                                                      ignore.strand=FALSE,
                                                      seqinfo=NULL,
                                                      seqnames.field=c("chromosome", "chrom",
                                                                       "chr", "chromosome_name"),
                                                      start.field="start",
                                                      end.field=c("end", "stop"),
                                                      strand.field="strand",
                                                      starts.in.df.are.0based=FALSE)
peaks.MB.H3K18la_anno <- annotatePeak(peaks.MB.H3K18la_ggranges,
                                      tssRegion=c(-2000, 2000),
                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      annoDb='org.Mm.eg.db')  
genes.w.MB.H3K18la_promoter_peak <- data.frame(peaks.MB.H3K18la_anno@anno$ENSEMBL,peaks.MB.H3K18la_anno@anno$SYMBOL)[peaks.MB.H3K18la_anno@detailGenomicAnnotation$Promoter==T,]
genes.w.MB.H3K18la_promoter_peak <- unique(genes.w.MB.H3K18la_promoter_peak)

peaks.MT.H3K18la <- read.table("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Eva_analysis/used_Peaks/peaks/MT_H3K18la_merged_0.01.stringent.bed.w.chr.sorted.bed")
colnames(peaks.MT.H3K18la) <- c('chr','start','end')
peaks.MT.H3K18la_ggranges <- makeGRangesFromDataFrame(peaks.MT.H3K18la,
                                                      keep.extra.columns=TRUE,
                                                      ignore.strand=FALSE,
                                                      seqinfo=NULL,
                                                      seqnames.field=c("chromosome", "chrom",
                                                                       "chr", "chromosome_name"),
                                                      start.field="start",
                                                      end.field=c("end", "stop"),
                                                      strand.field="strand",
                                                      starts.in.df.are.0based=FALSE)
peaks.MT.H3K18la_anno <- annotatePeak(peaks.MT.H3K18la_ggranges,
                                      tssRegion=c(-2000, 2000),
                                      TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      annoDb='org.Mm.eg.db')  
genes.w.MT.H3K18la_promoter_peak <- data.frame(peaks.MT.H3K18la_anno@anno$ENSEMBL,peaks.MT.H3K18la_anno@anno$SYMBOL)[peaks.MT.H3K18la_anno@detailGenomicAnnotation$Promoter==T,]
genes.w.MT.H3K18la_promoter_peak <- unique(genes.w.MT.H3K18la_promoter_peak)

rawCounts <- read.table("/Volumes/green_groups_nme_public/EvaGalle/Data/RNA/Sequencing/Counts/counts.SEQ00040_and_51.txt", header=TRUE, sep="\t")

# Current pipeline 
# Remove genes with no count across all samples
rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]
rawCounts <- rawCounts[,grep('MB',colnames(rawCounts))]

samples <- colnames(rawCounts)
n<-length(samples)
supplementation <- unlist(strsplit(samples,"_"))[c(1:n)*6-1]
condition <- unlist(strsplit(samples,"_"))[c(1:n)*6-2]
replicate <- unlist(strsplit(samples,"_"))[c(1:n)*6]
seqnr <- unlist(strsplit(samples,"_"))[c(1:n)*6-5]

# Add sample information
metadata <- data.frame(condition, supplementation, replicate, seqnr,row.names = samples)

### RPM
# Remove genes with < 1 read count
 keep       <- rowSums(rawCounts) >= 10
 rawCounts  <- rawCounts[keep,]

# Calculate RPM values
libSize <- Matrix::colSums(rawCounts)

rpm <- sweep(rawCounts * 1e6, MARGIN=2, STATS=libSize, FUN="/")

normCounts <- log2(rpm + 1)

# Define group
group <- paste0(condition,"_",supplementation)
#group <- paste0(supplementation)
#group<- factor(group)

# Create DGEList object
y <- DGEList(counts=rawCounts, group=group)


## Filter low expressed genes
# Remove genes with low expression in > 50% samples
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

```

## Normalize counts
```{r}

# Apply TMM normalization
y <- calcNormFactors(y)

```

## Differential gene expression analysis
```{r}

# Construct design matrix
 design <- model.matrix(~ 0 + group+seqnr)
# design <- model.matrix(~ 0 + group+condition+seqnr)

# Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)                   # y$common.dispersion <- 0.1 (samples with no replicates)

# Estimate quasi-likelihood (QL) dispersions
fitGlm <- glmQLFit(y, design, robust = TRUE)

# Define contrast of interest
contr <- makeContrasts(grouphg_lac - grouphg_c, levels=design)

# Apply QL F-test for differential expression
qlf <- glmQLFTest(fitGlm, contrast=contr)

# Extract HFD vs CHD DE results
res <- topTags(qlf, n = Inf, sort.by = "PValue")$table               

# Gene mapping
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")

geneAnnotation <- getBM(attributes=c("ensembl_gene_id","entrezgene_id", "external_gene_name"), 
                    filters="ensembl_gene_id", 
                    values=rownames(res), mart=ensembl)

# Combine DE results with gene names
comb <- merge(geneAnnotation, data.frame(res), 
              by.x="ensembl_gene_id", 
              by.y="row.names")

# Rename columns
colnames(comb)[1:2] <- c("geneID", "geneName")

# Annotate genes based on log2FC / adjusted p-value thresholds
comb$regulation <- ""
comb[which(comb$logFC > 0.5 & comb$PValue <= 0.01),"regulation"] <- "up"
comb[which(comb$logFC < (-0.5) & comb$PValue <= 0.01),"regulation"] <- "down"
comb[which(comb$PValue > 0.01),"regulation"] <- "non-significant"
comb$regulation[comb$regulation==""]<- "non-significant"
# Select significant DE genes based on FDR corrected p values
sig_deg <- subset(comb, PValue < 0.05 & abs(logFC) > 0.5)

# add H3K18la promoter peak info to the comb dataframe
comb$H3K18la.promoter.peak_MT <- "No"
comb[which(comb$geneID %in% genes.w.MT.H3K18la_promoter_peak$peaks.MT.H3K18la_anno.anno.ENSEMBL),"H3K18la.promoter.peak_MT"] <- "Yes"
comb$H3K18la.promoter.peak_MB <- "No"
comb[which(comb$geneID %in% genes.w.MB.H3K18la_promoter_peak$peaks.MB.H3K18la_anno.anno.ENSEMBL),"H3K18la.promoter.peak_MB"] <- "Yes"

comb$H3K18la.promoter.peak_MT <- "No"
comb[which(comb$geneID %in% genes.w.MT.H3K18la_promoter_peak$peaks.MT.H3K18la_anno.anno.ENSEMBL),"H3K18la.promoter.peak_MT"] <- "Yes"
comb$H3K18la.promoter.peak_MB <- "No"
comb[which(comb$geneID %in% genes.w.MB.H3K18la_promoter_peak$peaks.MB.H3K18la_anno.anno.ENSEMBL),"H3K18la.promoter.peak_MB"] <- "Yes"
DEGs<-list()
DEGs[[i]] <- comb
library(ggplot2)
ggplot(DEGs[[1]], aes(x = logFC, y = -log10(PValue), color = regulation, shape = H3K18la.promoter.peak_MT)) +
  geom_point() + theme_minimal() +
  scale_color_manual(values = c("#00AFBB","#DCDCDC","#FC4E07")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold"), 
        legend.title = element_text(size=10, face="italic")) +
  labs(x="log fold change", y="-log(nom p-value)") +
  geom_vline(xintercept=c(-1, 1), lty=2, colour = "#323232") +
  geom_hline(yintercept=-log10(0.05), lty=2, colour = "#323232" ) 

#number of genes with hyper- and hypolactylated promoters
table(DEGs[[1]]$regulation)

# Read RNA 
RNA_Pr_MB_MT<-read.table("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Eva_analysis/RNA_PromH3K18la_MB_MT.txt")

DEGs[[1]] <- merge(DEGs[[1]],RNA_Pr_MB_MT, by.x = 'geneID',by.y='ensembl_gene_id')


pdf("/Volumes/green_groups_nme_public/_Projects/MS_H3K18la/Revision/GitHub/try_output/Fig3H.RNA.logFC(MB+Lact).per.promoter.peak.pdf",width=2,height =4)
boxplot(DEGs[[i]]$logFC~DEGs[[i]]$H3K18la.promoter.peak_MB+DEGs[[i]]$H3K18la.promoter.peak_MT, outline=T,col=c('grey','lightblue','lightblue','blue'),notch=T, pch=20
        ,xlab = 'H3K18la promoter peak in MB.MT',ylab = 'RNA logFC',outcol=alpha(c('grey','lightblue','lightblue','blue'),0.1))
abline(h=0,lty=2)
dev.off()



                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
