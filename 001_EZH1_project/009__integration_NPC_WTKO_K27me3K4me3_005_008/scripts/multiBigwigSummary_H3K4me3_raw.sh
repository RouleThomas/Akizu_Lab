#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bw output/bigwig/NPC_KO_H3K4me3_005.unique.dupmark.sorted.bw output/bigwig/NPC_KO_H3K4me3_008.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_H3K4me3_raw.npz


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_H3K4me3_raw.npz \
    --transpose \
    --ntop 0 \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --colors black darkgrey red darkred \
    -o output/bigwig/multiBigwigSummary_H3K4me3_raw_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_H3K4me3_raw.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_H3K4me3_raw_heatmap.pdf





