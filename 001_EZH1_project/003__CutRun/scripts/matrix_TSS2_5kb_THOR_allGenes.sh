#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS2_5kp_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS2_5kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_5kp_THOR_allGenes_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList white,grey,black \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15



# 20240119 Naiara tasks


plotHeatmap -m output/deeptools/matrix_TSS2_5kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_5kp_THOR_allGenes_heatmap2_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS2_5kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_5kp_THOR_allGenes_heatmap2_colorSmall_kmeans20.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --kmeans 20 \
    --heatmapHeight 20 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS2_5kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_5kp_THOR_allGenes_heatmap2_colorSmall_kmeans8.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --kmeans 8 \
    --heatmapHeight 20 \
    --heatmapWidth 2


