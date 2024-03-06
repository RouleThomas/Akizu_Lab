#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_macs2_H3K27me3_WTKO_pool_qval2.30103_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotProfile -m output/deeptools/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_5kb_H3K27me3_median_THOR_genePeaks_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



