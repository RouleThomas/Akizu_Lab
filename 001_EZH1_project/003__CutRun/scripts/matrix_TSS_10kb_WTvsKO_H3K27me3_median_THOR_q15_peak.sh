#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00

computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_WTvsKO_unique_Keepdup/THOR_qval15.bed \
    -S output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_q15_peak_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


