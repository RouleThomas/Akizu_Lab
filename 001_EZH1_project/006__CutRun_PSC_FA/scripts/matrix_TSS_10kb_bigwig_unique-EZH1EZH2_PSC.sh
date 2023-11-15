#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI-EZH1EZH2__EZH1only.gtf meta/ENCFF159KBI-EZH1EZH2__EZH1andEZH2.gtf meta/ENCFF159KBI-EZH1EZH2__EZH2only.gtf \
    -S output/bigwig/PSC_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_unique-EZH1EZH2_PSC.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_unique-EZH1EZH2_PSC.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique-EZH1EZH2_PSC_heatmap.pdf \
    --samplesLabel "EZH1" "EZH2" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'plot, heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2


   # --zMin 0 0 0 --zMax 30 30 30 

