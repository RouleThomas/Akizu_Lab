#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_Venn_overlap_SICER2_Gain_1249.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_Venn_overlap-SICER2_Gain_1249.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap-SICER2_Gain_1249.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap-SICER2_Gain_1249_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotProfile -m output/deeptools/matrix_TSS_5kb_Venn_overlap-SICER2_Gain_1249.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap-SICER2_Gain_1249_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS"



