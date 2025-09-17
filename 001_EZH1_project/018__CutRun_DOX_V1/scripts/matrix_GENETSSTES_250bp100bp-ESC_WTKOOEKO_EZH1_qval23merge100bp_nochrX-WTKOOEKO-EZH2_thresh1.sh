#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/gencode_ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_heatmap.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_plotProfile2.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_plotProfile3.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4



# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_heatmap1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3




plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp_nochrX-WTKOOEKO-EZH2_thresh1_heatmap2.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 2 2 2
