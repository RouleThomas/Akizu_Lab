#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s1-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s1-rep1.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s1-rep2.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s2-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s2-rep1.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_TMM/PSCWTvsKOH3K27me3TMM-s2-rep2.bw -o output/deeptools/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM.npz


# Plot
## PCA
plotPCA -in output/deeptools/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 \
    --colors black black black red red red \
    -o output/deeptools/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/deeptools/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/deeptools/multiBigwigSummary_H3K27me3_WTvsKO_THOR_TMM_heatmap.pdf





