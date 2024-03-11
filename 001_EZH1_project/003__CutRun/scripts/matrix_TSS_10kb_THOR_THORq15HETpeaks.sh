#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_WTvsHET_unique_Keepdup/THOR_qval15.bed \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15

plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15

plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap_kmean2.png \
    --samplesLabel "WT" "HET" "KO" \
    --kmeans 2 \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15


plotProfile -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_profile.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --perGroup \
    --colors black blue red \
    --refPointLabel "0" \
    -T "H3K27me3 read density" \
    -z ""


# 20240310 Naiara plot Slack


plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap1.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap2.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15HETpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15HETpeaks_heatmap3.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue \
    --colorNumber 3