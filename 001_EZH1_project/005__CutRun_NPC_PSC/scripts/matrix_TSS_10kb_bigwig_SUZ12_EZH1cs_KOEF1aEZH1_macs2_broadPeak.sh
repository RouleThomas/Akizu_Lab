#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak \
    -S output/bigwig/PSC_KOEF1aEZH1_SUZ12.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH1cs.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans2.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans4.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 4

plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans8.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 8

plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans20.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 20 \
    --outFileSortedRegions output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans20.bed



plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.pdf \
    --samplesLabel "SUZ12" "EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 25 \
    --outFileSortedRegions output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.bed
