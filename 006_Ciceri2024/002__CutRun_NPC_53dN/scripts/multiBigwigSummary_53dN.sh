#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/53dN_WT_H3K27ac_R1.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K27ac_R2.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K9me3_R1.unique.dupmark.sorted.bw output/bigwig/53dN_WT_H3K9me3_R2.unique.dupmark.sorted.bw output/bigwig/53dN_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/53dN_WT_IGG_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_53dN.npz


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_53dN.npz \
    --transpose \
    --ntop 0 \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --colors "#00bfff" "#87cefa" "#8b0000" "#ff0000" "#006400" "#32cd32" "#9400d3" "#9370db" "#8b4513" "#cd853f" \
    --markers o o o o o o o o o o \
    -o output/bigwig/multiBigwigSummary_53dN_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_53dN.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_53dN_heatmap.pdf





