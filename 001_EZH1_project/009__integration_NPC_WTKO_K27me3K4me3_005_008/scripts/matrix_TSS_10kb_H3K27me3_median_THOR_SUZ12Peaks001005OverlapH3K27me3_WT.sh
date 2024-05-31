#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT_heatmap.pdf \
    --samplesLabel "WT" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT_heatmap2.pdf \
    --samplesLabel "WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT_heatmap3.pdf \
    --samplesLabel "WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10


plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WT_profile.pdf \
    --samplesLabel "WT" \
    --perGroup \
    --colors black \
    --refPointLabel "center" \
    -T "H3K27me3 read density" \
    -z "" \
    --yMax 100 \
    --yMin 30 \
    --plotWidth 7