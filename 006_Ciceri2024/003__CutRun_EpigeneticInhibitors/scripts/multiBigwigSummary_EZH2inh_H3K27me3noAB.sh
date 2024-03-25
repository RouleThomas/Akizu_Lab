#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/DMSO_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/DMSO_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/DMSO_noAB_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/EZH2inh_noAB_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_noAB_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_EZH2inh_H3K27me3noAB.npz




# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_EZH2inh_H3K27me3noAB.npz \
    --transpose \
    --ntop 0 \
    --labels DMSO_H3K27me3_R1 DMSO_H3K27me3_R2 DMSO_noAB_R1 EZH2inh_H3K27me3_R1 EZH2inh_H3K27me3_R2 EZH2inh_noAB_R1 EZH2inh_noAB_R2 \
    --colors "#1e90ff" "#1e90ff" "#8b4513" "#6b8e23" "#6b8e23" "#daa520" "#daa520" \
    --markers o o o o o o o \
    -o output/bigwig/multiBigwigSummary_EZH2inh_H3K27me3noAB_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_EZH2inh_H3K27me3noAB.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels DMSO_H3K27me3_R1 DMSO_H3K27me3_R2 DMSO_noAB_R1 EZH2inh_H3K27me3_R1 EZH2inh_H3K27me3_R2 EZH2inh_noAB_R1 EZH2inh_noAB_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_EZH2inh_H3K27me3noAB_heatmap.pdf





