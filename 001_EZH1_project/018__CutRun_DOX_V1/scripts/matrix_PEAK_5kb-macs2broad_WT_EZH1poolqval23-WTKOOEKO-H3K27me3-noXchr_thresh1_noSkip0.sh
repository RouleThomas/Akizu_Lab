#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WT_EZH1_noXchr_pool_peaks.broadPeak \
    -S output/bigwig_Ferguson/ESC_WT_H3K27me3_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_KO_H3K27me3_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_OEKO_H3K27me3_noXchr_unique_norm99_initialBigwig_thresh1_median.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0_plotProfile1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8






# interactive


plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WT_EZH1poolqval23-WTKOOEKO-H3K27me3-noXchr_thresh1_noSkip0_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2




