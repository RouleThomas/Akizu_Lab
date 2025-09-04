#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf \
    -S output/bigwig_Ferguson/ESC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_KO_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_heatmap.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2




plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_plotProfile2.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8


# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_heatmap1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "OEKO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 1 1 1




