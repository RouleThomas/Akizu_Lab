#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15



plotProfile -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_profile.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --perGroup \
    --colors black blue red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


# 202309 Naiara tasks
plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2_color.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2