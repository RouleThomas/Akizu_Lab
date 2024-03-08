#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/PSCWTvsKOEF1aEZH1H3K27me3-s1-rep0.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/PSCWTvsKOEF1aEZH1H3K27me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes_heatmap2.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes_heatmap3.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4



plotProfile -m output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes_profile.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



