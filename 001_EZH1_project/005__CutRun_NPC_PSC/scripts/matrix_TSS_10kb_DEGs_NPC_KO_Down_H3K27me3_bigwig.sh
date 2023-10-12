#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Down.gtf \
    -S output/bigwig/NPC_WT_H3K27me3.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig.gz \
    -out output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig.gz \
    -out output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


