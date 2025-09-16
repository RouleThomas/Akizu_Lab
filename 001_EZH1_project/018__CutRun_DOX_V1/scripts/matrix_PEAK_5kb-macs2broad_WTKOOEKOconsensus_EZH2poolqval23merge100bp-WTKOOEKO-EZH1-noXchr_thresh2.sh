#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_noXchr_pool_peaks.sorted.merge100bp.bed \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_KO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8



# interactive


plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3



plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-EZH1-noXchr_thresh2_heatmap2.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 0.5 0.5 0.5



