#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Down.gtf \
    -S output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-s1-rep0.bw \
    output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR.gz \
    -out output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR.gz \
    -out output/deeptools/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K4me3 read density" \
    -z ""


