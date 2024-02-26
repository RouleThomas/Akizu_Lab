#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3_005.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3_008.unique.dupmark.sorted.bw -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K27me3_DiffBindTMM.npz


# Plot
## PCA
plotPCA -in output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K27me3_DiffBindTMM.npz \
    --transpose \
    --ntop 0 \
    --labels WT_005 WT_008 KO_005 KO_008 \
    -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K27me3_DiffBindTMM_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K27me3_DiffBindTMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K27me3_DiffBindTMM_heatmap.pdf





