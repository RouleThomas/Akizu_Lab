#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorList bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes_heatmap2.png \
    --samplesLabel "WT" "KO" \
    --colorList bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKO_H3K27me3_median_THOR_allGenes_heatmap3.png \
    --samplesLabel "WT" "KO" \
    --colorList bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


