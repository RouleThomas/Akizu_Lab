#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bpnoXchrthresh1-ESC_KO_vs_ESC_WT-EZH2.bed output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bpnoXchrthresh1-ESC_KO_vs_ESC_WT-EZH2.bed \
    -S output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1_heatmap.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colors black red blue \
    --perGroup \
    --plotWidth 6



# interactive
plotHeatmap -m output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1.gz \
    -out output/deeptools/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2-noXchr_thresh1_heatmap1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 1 1 1 



