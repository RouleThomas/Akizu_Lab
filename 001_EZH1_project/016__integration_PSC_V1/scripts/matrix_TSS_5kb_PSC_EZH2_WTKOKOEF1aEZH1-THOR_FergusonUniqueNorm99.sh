#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEZH2FergusonUniqueNorm99noInput-s1_median.bw output/THOR/THOR_PSC_WTvsKO_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEZH2FergusonUniqueNorm99noInput-s2_median.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1EZH2FergusonUniqueNorm99noInput-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive

plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 0.5


plotProfile -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-THOR_FergusonUniqueNorm99_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black red blue \
    --perGroup



