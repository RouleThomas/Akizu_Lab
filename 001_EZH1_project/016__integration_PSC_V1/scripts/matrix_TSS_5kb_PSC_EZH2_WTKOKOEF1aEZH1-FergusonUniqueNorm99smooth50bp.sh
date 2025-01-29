#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth50bp.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



# interactive

plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 0.5


plotProfile -m output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_EZH2_WTKOKOEF1aEZH1-FergusonUniqueNorm99smooth50bp_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KO_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black red blue \
    --perGroup





