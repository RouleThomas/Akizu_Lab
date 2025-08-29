#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/ENCFF159KBI.gtf \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_KO_EZH1_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4



# interactive
plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 1 1 1



