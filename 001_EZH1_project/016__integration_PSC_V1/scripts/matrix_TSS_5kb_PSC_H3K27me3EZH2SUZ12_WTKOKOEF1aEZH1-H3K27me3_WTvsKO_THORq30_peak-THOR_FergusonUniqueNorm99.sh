#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_positive.bed output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/THOR_qval30_negative.bed \
    -S output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.bw \
    output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEZH2FergusonUniqueNorm99noInput-s1_median.bw output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEZH2FergusonUniqueNorm99noInput-s2_median.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1EZH2FergusonUniqueNorm99noInput-s2_median.bw \
    output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/PSCWTvsKOSUZ12FergusonUniqueNorm99noInput-s1_median.bw output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput/PSCWTvsKOSUZ12FergusonUniqueNorm99noInput-s2_median.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1SUZ12FergusonUniqueNorm99noInput-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-THOR_FergusonUniqueNorm99.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-THOR_FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_H3K27me3EZH2SUZ12_WTKOKOEF1aEZH1-H3K27me3_WTvsKO_THORq30_peak-THOR_FergusonUniqueNorm99_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "KOEF1aEZH1_H3K27me3" "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


