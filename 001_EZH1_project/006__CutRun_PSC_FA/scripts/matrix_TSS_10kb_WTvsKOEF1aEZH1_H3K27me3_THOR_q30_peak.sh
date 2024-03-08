#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/THOR_qval30.bed \
    -S output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/PSCWTvsKOEF1aEZH1H3K27me3-s1-rep0.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/PSCWTvsKOEF1aEZH1H3K27me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak.gz \
    -p 5



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak_heatmap.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak_heatmap1.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak_heatmap2.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4