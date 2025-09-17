#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_OEKO_vs_ESC_WT-H3K27me3.bed output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_OEKO_vs_ESC_WT-H3K27me3.bed \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_KO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 6



# interactive
plotHeatmap -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1-noXchr_thresh2_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 1 1 1



