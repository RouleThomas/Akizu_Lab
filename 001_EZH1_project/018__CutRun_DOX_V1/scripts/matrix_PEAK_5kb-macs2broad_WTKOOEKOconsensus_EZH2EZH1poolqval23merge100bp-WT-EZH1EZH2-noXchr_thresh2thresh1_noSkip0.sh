#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2EZH1_noXchr_pool_peaks.sorted.merge100bp.bed \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_heatmap.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colors black grey \
    --perGroup \
    --plotWidth 8



plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3



plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_heatmap2.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 0.5 0.5


# clustering




plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_heatmap2kmeans3.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 0.5 0.5 \
    --kmeans 3




plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_plotProfile1kmeans3.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colors black grey \
    --perGroup \
    --plotWidth 8 \
    --kmeans 3





plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2EZH1poolqval23merge100bp-WT-EZH1EZH2-noXchr_thresh2thresh1_noSkip0_plotProfile1kmeans8.pdf \
    --samplesLabel "WT_EZH1" "WT_EZH2" \
    --colors black grey \
    --perGroup \
    --plotWidth 8 \
    --kmeans 8


