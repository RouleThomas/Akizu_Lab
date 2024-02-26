#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1-rep1.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2-rep1.bw -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR.npz


# Plot
## PCA
plotPCA -in output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR.npz \
    --transpose \
    --ntop 0 \
    --labels WT_005 WT_008 KO_005 KO_008 \
    -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/THOR/THOR_NPC_WTvsKO_H3K27me3/multiBigwigSummary_H3K27me3_THOR_heatmap.pdf





