#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50.bed \
    -S output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s1-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -p 5



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak_heatmap1.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10


plotProfile -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
