#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-FergusonUniqueNorm99_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


