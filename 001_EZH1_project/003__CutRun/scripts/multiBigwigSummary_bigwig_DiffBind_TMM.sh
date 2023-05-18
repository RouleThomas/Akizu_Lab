#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
-o output/bigwig_DiffBind_TMM/multiBigwigSummary_DiffBind_TMM.npz



plotCorrelation \
    -in output/bigwig_DiffBind_TMM/multiBigwigSummary_DiffBind_TMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig_DiffBind_TMM/multiBigwigSummary_DiffBind_TMM_heatmap.pdf




