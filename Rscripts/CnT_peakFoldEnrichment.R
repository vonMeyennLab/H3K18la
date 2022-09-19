
# This R script has been used to generate the following figures:
### Figures 1E, 4C

library(GenomicRanges)
library(pheatmap)

# Load sample metadata
# Sample metadata table has the following columns:
# ID  tissue  histone replicate bamFilePath peakFilePath
metadata <- read.table("../metadata.txt", header = TRUE, sep = "\t")

# List all the samples
samples <- metadata$ID

# Load mm10 genomic annotation tracks
genome_annotation <- readRDS("../genomeAnnotation_mm10.rds")

# Select only the genomic annotation tracks of interest
genome_annotation_subset <-
  genome_annotation %>% .[c("CGI Promoters", "Non-CGI Promoters", "3' UTRs", "5' UTRs",
                            "Exons", "Introns", "Intergenic regions", 
                            "PLS", "pELS", "dELS", "CTCF", "DNase_H3K4me3" )] %>% GRangesList()

# Function calculating genome size based on genic and intergenic regions
genome_size <- 
  sum(
    GenomicRanges::reduce(genome_annotation$`Gene bodies`, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum,
    GenomicRanges::reduce(genome_annotation$`Intergenic regions`, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
  )

calc_fold_enrichment <- 
  function(gr1, gr2, genome_size){
    
    A <- 
      GenomicRanges::intersect(
        GenomicRanges::reduce(gr2, ignore.strand = TRUE),
        GenomicRanges::reduce(gr1, ignore.strand = TRUE)
      ) %>% GenomicRanges::width() %>% sum
    
    B <- 
      GenomicRanges::reduce(gr1, ignore.strand = TRUE) %>%
      GenomicRanges::width() %>% sum
    
    C <- 
      GenomicRanges::reduce(gr2, ignore.strand = TRUE) %>% 
      GenomicRanges::width() %>% sum
    
    log2((A/B)/(C/genome_size))
  }

# Calculate peak fold enrichment of all samples
sampleFE <- list()

for (sample in samples)
{
  hist <- subset(peakSummary, ID == sample)
  hist$chr <- gsub("chr", "", hist$V1)
  histGR <- GRanges(hist$V1, IRanges(start = hist$V2, end = hist$V3), strand = "*")
  
  histFE <- lapply(genome_annotation_subset, function(x){calc_fold_enrichment(histGR, x, genome_size = genome_size)})
  sampleFE[[sample]] <- histFE
}

sampleFE <- lapply(sampleFE, function(x){x <- do.call("rbind", x)})

# Combine peak fold enrichment values of all samples in table format
comb_FE <- do.call("cbind", sampleFE)
colnames(comb_FE) <- names(sampleFE)

# Group genomic regions into genomic features and cCREs
genomicFeatures <- c("CGI Promoters","Non-CGI Promoters", "Exons", "Introns", "3' UTRs", "Intergenic regions")
cCREs <- c("CTCF", "DNase_H3K4me3", "PLS","pELS", "dELS")

# Heat map visualization of peak fold enrichment of all samples
breaksList <- seq(-2, 2, by = 0.1)

pheatmap::pheatmap(t(comb_FE[genomicFeatures,samples]),
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "row",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(length(breaksList)),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 40, 
                   cellwidth         = 40,
                   border_color      = 'black',
                   breaks            = breaksList)

pheatmap::pheatmap(t(comb_FE[cCREs,samples]),
                   silent            = F,
                   cutree_rows       = 2,
                   cutree_cols       = 2,
                   scale             = "row",
                   color             = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(length(breaksList)),
                   fontsize_row      = 10, 
                   fontsize_col      = 10,
                   display_numbers   = FALSE,
                   fontsize_number   = 10, 
                   cluster_cols      = FALSE,
                   cluster_rows      = FALSE,
                   cellheight        = 40, 
                   cellwidth         = 40,
                   border_color      = 'black',
                   breaks            = breaksList)
