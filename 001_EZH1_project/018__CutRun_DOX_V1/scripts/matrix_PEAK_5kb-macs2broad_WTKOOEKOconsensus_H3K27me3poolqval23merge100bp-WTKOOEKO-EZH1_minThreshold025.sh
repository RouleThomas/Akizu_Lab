#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_KO_EZH1_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_unique_norm99_initialBigwig_median.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -p 6 \
    --minThreshold 0.25





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_plotProfile2.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8 \
    --yMax 0.5



# interactive


plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 0.5 0.5 0.5




plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_heatmap3.pdf \
    --kmeans 3 --zMin 0 --zMax 0.5 \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --heatmapHeight 10 \
    --heatmapWidth 2




plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_heatmap3.pdf \
    --kmeans 3 \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr




plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_heatmap4.pdf \
    --kmeans 3 \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --heatmapWidth 2





plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025_plotProfile3.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8 \
    --yMin 1.5


