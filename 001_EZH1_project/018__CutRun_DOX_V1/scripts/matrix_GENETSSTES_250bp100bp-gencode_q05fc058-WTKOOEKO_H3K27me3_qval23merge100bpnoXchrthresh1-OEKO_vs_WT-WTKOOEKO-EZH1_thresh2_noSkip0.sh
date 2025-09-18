#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/gencode_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf meta/gencode_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf \
    -S output/bigwig_Ferguson/ESC_WT_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_KO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw output/bigwig_Ferguson/ESC_OEKO_EZH1_noXchr_unique_norm99_initialBigwig_thresh2_median.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_heatmap.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_plotProfile1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 8

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_plotProfile2.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_plotProfile3.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colors black red blue \
    --perGroup \
    --plotWidth 4



# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_heatmap1.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3




plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-OEKO_vs_WT-WTKOOEKO-EZH1_thresh2_noSkip0_heatmap2.pdf \
    --samplesLabel "WT_EZH1" "KO_EZH1" "OEKO_EZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 1 1 1 
