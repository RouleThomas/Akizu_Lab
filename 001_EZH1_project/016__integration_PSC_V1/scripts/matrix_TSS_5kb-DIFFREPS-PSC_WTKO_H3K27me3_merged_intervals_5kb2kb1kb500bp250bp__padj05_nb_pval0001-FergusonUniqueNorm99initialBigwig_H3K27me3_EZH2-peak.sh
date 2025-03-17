#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Gain.txt output/diffreps/merged_intervals-padj05_nb_pval0001-5kb2kb1kb500bp250bp-WTvsKO_initialBigwig-Lost.txt \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 4 4



plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_nb_pval0001-FergusonUniqueNorm99initialBigwig_H3K27me3_EZH2-peak_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 2 2


