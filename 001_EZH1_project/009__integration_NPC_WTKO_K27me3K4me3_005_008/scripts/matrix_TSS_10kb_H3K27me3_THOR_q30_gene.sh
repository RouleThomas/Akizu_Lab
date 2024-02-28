#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00




computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_H3K27me3_q30_pos_Promoter_5.gtf meta/ENCFF159KBI_H3K27me3_q30_neg_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep1.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep1.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_q30_gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_q30_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_q30_gene_heatmap.pdf \
    --samplesLabel "WT_005" "WT_008" "KO_005" "KO_008" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


