#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_H3K27me3_q30_pos_Promoter_5.gtf meta/ENCFF159KBI_H3K27me3_q30_neg_Promoter_5.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene_heatmap1.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2 




plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4 



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10 



plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_median_THOR_gainLost_gene_profile.pdf \
    --samplesLabel "WT" "KO" \
    --colors black darkgrey \
    -T "Read density"



