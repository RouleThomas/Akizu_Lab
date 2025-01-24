#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


