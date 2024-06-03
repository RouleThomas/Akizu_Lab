#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../005__CutRun_NPC_PSC/output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks_overlap_001009_NPC_WT_H3K27me3_broad2.3.broadPeak \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10


plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_SUZ12Peaks001005OverlapH3K27me3_WTKOwithinWTpeaks_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "center" \
    -T "H3K27me3 read density" \
    -z ""


