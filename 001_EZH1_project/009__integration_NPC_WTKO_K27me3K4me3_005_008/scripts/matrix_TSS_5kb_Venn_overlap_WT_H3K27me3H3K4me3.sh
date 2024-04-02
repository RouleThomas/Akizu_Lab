#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3only.gtf meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3.gtf meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3only.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_heatmap_colorSmall.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_profile.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --perGroup \
    --colors orange purple \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""



