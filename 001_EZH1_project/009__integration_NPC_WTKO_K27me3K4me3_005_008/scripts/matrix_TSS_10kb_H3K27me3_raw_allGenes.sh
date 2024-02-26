#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_005.unique.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_008.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_raw_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_raw_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_raw_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "WT_005" "WT_008" "KO_005" "KO_008" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_raw_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_raw_allGenes_profile.pdf \
    --samplesLabel "WT_005" "WT_008" "KO_005" "KO_008" \
    --perGroup \
    --colors black darkgrey red darkred \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



