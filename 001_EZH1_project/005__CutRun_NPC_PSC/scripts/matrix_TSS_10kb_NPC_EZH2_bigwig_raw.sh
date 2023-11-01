#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig/NPC_WT_EZH2.dupmark.sorted.bw output/bigwig/NPC_KO_EZH2.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_NPC_EZH2_bigwig_raw.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_NPC_EZH2_bigwig_raw.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_EZH2_bigwig_raw_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2


