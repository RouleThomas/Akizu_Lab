#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval3/PSC_KOEF1aEZH1_EZH1_pool_peaks.broadPeak \
    -S output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive
plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 20

plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 0.5


plotProfile -m output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black blue \
    --perGroup


plotProfile -m output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOEF1aEZH1-initialBigwig_EZH2-peakKOEF1aEZH1_EZH1macs2q3_plotProfile2.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black blue \
    --perGroup \
    --plotHeight 10 \
    --plotWidth 7

