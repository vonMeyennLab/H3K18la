# H3K18 lactylation marks tissue-specific active enhancers 

This repository contains all scripts to reproduce the results of the manuscript [H3K18 lactylation marks tissue-specific active enhancers](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02775-y). 

### Abstract

#### Background

Histone lactylation has been recently described as a novel histone post-translational modification linking cellular metabolism to epigenetic regulation.

#### Results

Given the expected relevance of this modification and current limited knowledge of its function, we generate genome-wide datasets of H3K18la distribution in various in vitro and in vivo samples, including mouse embryonic stem cells, macrophages, adipocytes, and mouse and human skeletal muscle. We compare them to profiles of well-established histone modifications and gene expression patterns. Supervised and unsupervised bioinformatics analysis shows that global H3K18la distribution resembles H3K27ac, although we also find notable differences. H3K18la marks active CpG island-containing promoters of highly expressed genes across most tissues assessed, including many housekeeping genes, and positively correlates with H3K27ac and H3K4me3 as well as with gene expression. In addition, H3K18la is enriched at active enhancers that lie in proximity to genes that are functionally important for the respective tissue.

#### Conclusions

Overall, our data suggests that H3K18la is not only a marker for active promoters, but also a mark of tissue specific active enhancers.

### Directory structure

- `/Rscripts/`: R scripts for reproducing all figures of the manuscript
- `/ChromHMM/*.sh`: Commands for Cut&Tag based ChromHMM analysis
- `/ChromHMM/*.bed`: Chromatin state regions for each ChromHMM analysis
    - `mmESC.ser_GRCm38_7_segments.bed`: List of chromatin states defined based on H3K18la, H3K27ac, H3K4me3 and H3K27me3 for mouse ESC serum/lif samples
    - `mmGAS_GRCm38_7_segments.bed`: List of chromatin states defined based on H3K18la, H3K27ac, H3K4me3 and H3K27me3 for mouse GAS samples
    - `mmPIM_GRCm38_7_segments.bed`: List of chromatin states defined based on H3K18la, H3K27ac, H3K4me3 and H3K27me3 for mouse PIM samples
    - `mmH3K18la_GRCm38_10_segments.bed`: List of chromatin states defined based on all mouse H3K18la samples
    - `hsMuscle_GRCh38_7_segments.bed`: List of chromatin states defined based on H3K18la, H3K27ac, H3K4me3, H3K27me3 and H3K9me3 for human muscle samples
