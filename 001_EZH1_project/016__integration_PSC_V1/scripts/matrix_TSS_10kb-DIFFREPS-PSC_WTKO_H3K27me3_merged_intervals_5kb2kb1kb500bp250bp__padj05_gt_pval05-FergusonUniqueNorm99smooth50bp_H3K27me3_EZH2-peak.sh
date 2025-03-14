#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Gain.txt output/diffreps/merged_intervals-padj05_gt_pval05-5kb2kb1kb500bp250bp-WTvsKO-Lost.txt \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth50bp.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 4 4



plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12 2 2



plotProfile -m output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak_plotProfile1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colors black red white white \
    --perGroup



plotProfile -m output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-PSC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj05_gt_pval05-FergusonUniqueNorm99smooth50bp_H3K27me3_EZH2-peak_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "WT_EZH2" "KO_EZH2" \
    --colors white white black red \
    --perGroup \
    --yMax 2

