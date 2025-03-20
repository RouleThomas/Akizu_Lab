#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Gain-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed output/diffreps/merged_intervals-padj001_nb_pval0001_log2FC1_initialBigwig-5kb2kb1kb500bp250bp-WTvsKO-Lost-macs2qval2.3_PSC_WTKO_EZH2_pool_peaks.bed \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2





plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 4 4



plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 2 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak_heatmap3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 125 125 15 15


plotProfile -m output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-WTKO_H3K27me3__padj001_nb_pval0001_log2FC1_initialBigwig-macs2qval2.3_WTKO_EZH2-H3K27me3_EZH2-peak_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colors white white black red \
    --perGroup \
    --plotWidth 4 \
    --yMax 6


