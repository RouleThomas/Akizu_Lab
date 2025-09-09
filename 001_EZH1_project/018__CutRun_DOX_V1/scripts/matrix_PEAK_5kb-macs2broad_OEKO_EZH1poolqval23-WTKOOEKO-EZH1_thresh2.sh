#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_EZH1_pool_peaks.broadPeak \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_unique_norm99_initialBigwig_median_thresh2.bw output/bigwig_Ferguson/ESC_KO_EZH1_unique_norm99_initialBigwig_median_thresh2.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_unique_norm99_initialBigwig_median_thresh2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4


plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2_plotProfile2.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8


# interactive


plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-EZH1_thresh2_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 0.5 0.5 0.5




