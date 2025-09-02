#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/ENCFF159KBI_ESC001019_upregulated_q05fc058_ESC_KO_vs_ESC_WT.gtf meta/ENCFF159KBI_ESC001019_downregulated_q05fc058_ESC_KO_vs_ESC_WT.gtf \
    -S output/bigwig_Ferguson/ESC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_KO_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_plotProfile1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colors black red blue black red blue black red blue \
    --perGroup \
    --plotWidth 8

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colors white white white black red blue white white white \
    --perGroup \
    --plotWidth 4

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_plotProfile3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colors black red blue white white white white white white \
    --perGroup \
    --plotWidth 4



# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3




plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 2 2 2 
