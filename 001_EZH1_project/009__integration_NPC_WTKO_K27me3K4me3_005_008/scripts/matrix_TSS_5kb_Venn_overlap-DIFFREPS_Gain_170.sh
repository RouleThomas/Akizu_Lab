#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_Venn_overlap_DIFFREPS_Gain_170.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_Venn_overlap-DIFFREPS_Gain_170.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap-DIFFREPS_Gain_170.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap-DIFFREPS_Gain_170_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotProfile -m output/deeptools/matrix_TSS_5kb_Venn_overlap-DIFFREPS_Gain_170.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap-DIFFREPS_Gain_170_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS"



