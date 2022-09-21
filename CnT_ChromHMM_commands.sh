##################################################
Stacked design
##################################################

# BinarizeBam
"java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
BinarizeBam \
-paired /cluster/work/nme/software/apps/ChromHMM/1.22/CHROMSIZES/mm10_nochr.txt \
/scratch/../ChromHMM_design_stacked_ESC_samples.txt \
/scratch/../binary_files" \
/

# LearnModel
"seq 1 30 | parallel java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
LearnModel \
-s 123 -nobrowser -noimage -nobed -noenrich -noautoopen \
/scratch/../binary_files/ \
/scratch/../model/ \
{} \
GRCm38" \
/

# CompareModels
"java -mx20000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
CompareModels \
/scratch/../model/emissions_30.txt \
/scratch/../model/ \
/scratch/../compare_models/stacked_all_samples_1_to_30" \
/

# MakeSegmentation
"ls /scratch/../model/model_7.txt | parallel java -mx5000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
MakeSegmentation \
{} \
/scratch/../binary_files \
/scratch/../segmentation" \
/

# OverlapEnrichment
"java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
segmentation/GRCm38_7_segments.bed \
/../ChromHMM/annotatr_mm/ \
annotatr_enrichment/ESC_genomicFeature_enrichment" \
/

# OverlapEnrichment
"java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
segmentation/GRCm38_7_segments.bed \
/../ChromHMM/ENCODE_mmESC/ \
cCRE_mESC_enrichment/ESC_cCRE_enrichment" \
/

# OverlapEnrichment
"java -mx4000M -jar -Djava.io.tmpdir=\$TMPDIR /cluster/work/nme/software/apps/ChromHMM/1.22/ChromHMM.jar \
OverlapEnrichment \
segmentation/GRCm38_7_segments.bed \
/../ChromHMM_publicData/ \
chromHMM_pData_enrichment/ESC_ChIPseq_enrichment" \
