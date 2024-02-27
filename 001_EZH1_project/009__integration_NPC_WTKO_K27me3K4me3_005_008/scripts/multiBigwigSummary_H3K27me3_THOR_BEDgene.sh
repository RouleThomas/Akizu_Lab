#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary BED-file -b output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep1.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep1.bw -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_BEDgene.npz --BED /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_gene.bed


# Plot
## PCA
plotPCA -in output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_BEDgene.npz \
    --transpose \
    --ntop 0 \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --colors black darkgrey red darkred \
    --markers o o o o \
    -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_BEDgene_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_BEDgene.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_BEDgene_heatmap.pdf





