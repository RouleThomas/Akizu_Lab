#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_macs2_H3K27me3_WTKO_pool_qval2.30103_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_heatmap_colorSmall2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_heatmap_colorSmall3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_genePeaks_heatmap_colorSmall4.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10
