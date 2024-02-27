#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary BED-file -b output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3_005.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3_008.unique.dupmark.sorted.bw -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene.npz --BED /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_gene.bed


# Plot
## PCA
plotPCA -in output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene.npz \
    --transpose \
    --ntop 0 \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --colors black darkgrey red darkred \
    --markers o o o o \
    -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_005 WT_008 KO_005 KO_008 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_MG1655_DiffBind_TMM/multiBigwigSummary_H3K4me3_DiffBindTMM_BEDgene_heatmap.pdf





