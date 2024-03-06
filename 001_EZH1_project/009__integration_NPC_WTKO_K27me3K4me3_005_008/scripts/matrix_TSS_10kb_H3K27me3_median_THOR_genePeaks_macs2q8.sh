#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_macs2_H3K27me3_WTKO_pool_qval8_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_macs2q8_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



