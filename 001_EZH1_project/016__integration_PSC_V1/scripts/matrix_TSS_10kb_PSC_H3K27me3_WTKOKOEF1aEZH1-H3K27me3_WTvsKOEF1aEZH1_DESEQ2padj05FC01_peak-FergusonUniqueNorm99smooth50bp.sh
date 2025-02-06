#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3_positivepadj05FC01.bed output/edgeR/DESEQ2-WTKOKOEF1aEZH1_H3K27me3_pool_peaks-PSC_KOEF1aEZH1_vs_PSC_WT-H3K27me3_negativepadj05FC01.bed \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth50bp.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive


plotProfile -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-H3K27me3_WTvsKOEF1aEZH1_DESEQ2padj05FC01_peak-FergusonUniqueNorm99smooth50bp_plotProfile.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colors black red blue \
    --perGroup

