#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_positive.bed output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_negative.bed \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth50bp.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth50bp.bw \
    output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth50bp.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" "KO_H3K27me3"  "WT_EZH2" "KOEF1aEZH1_EZH2" "KO_EZH2"  "WT_SUZ12" "KOEF1aEZH1_SUZ12" "KO_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 3 3 3 1 1 1 1 1 1


# interactive


plotProfile -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp_plotProfile1.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" "KO_H3K27me3"  "WT_EZH2" "KOEF1aEZH1_EZH2" "KO_EZH2"  "WT_SUZ12" "KOEF1aEZH1_SUZ12" "KO_SUZ12" \
    --colors white white white white white white black blue red \
    --perGroup \
    --yMax 0.4


plotProfile -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3EZH2SUZ12_WTKOEF1aEZH1KO-H3K27me3_WTvsKOEF1aEZH1_THORq30_peak-FergusonUniqueNorm99smooth50bp_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" "KO_H3K27me3"  "WT_EZH2" "KOEF1aEZH1_EZH2" "KO_EZH2"  "WT_SUZ12" "KOEF1aEZH1_SUZ12" "KO_SUZ12" \
    --colors white white white black blue red white white white \
    --perGroup \
    --yMax 0.4