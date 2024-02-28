#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_H3K4me3_q20_pos_Promoter_5.gtf meta/ENCFF159KBI_H3K4me3_q20_neg_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s1-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s1-rep1.bw output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s2-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s2-rep1.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_H3K4me3_THOR_q20_gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_H3K4me3_THOR_q20_gene.gz \
    -out output/deeptools/matrix_TSS_5kb_H3K4me3_THOR_q20_gene_heatmap.pdf \
    --samplesLabel "WT_005" "WT_008" "KO_005" "KO_008" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


