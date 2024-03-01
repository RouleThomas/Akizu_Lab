#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig_hg19/GSM5860155_CnR_H3K27ac_NPC_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860156_CnR_H3K27ac_NPC_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860159_CnR_H3K27me3_NPC_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860160_CnR_H3K27me3_NPC_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860163_CnR_H3K4me3_NPC_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860164_CnR_H3K4me3_NPC_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860167_CnR_H3K9me3_NPC_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860168_CnR_H3K9me3_NPC_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860171_CnR_IgG_NPC_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860172_CnR_IgG_NPC_rep2.hg19.sorted.RmDup.10mNorm.bw -o output/bigwig_hg19/multiBigwigSummary_NPC_Ciceri.npz


# Plot
## PCA
plotPCA -in output/bigwig_hg19/multiBigwigSummary_NPC_Ciceri.npz \
    --transpose \
    --ntop 0 \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --colors "#00bfff" "#87cefa" "#8b0000" "#ff0000" "#006400" "#32cd32" "#9400d3" "#9370db" "#8b4513" "#cd853f" \
    --markers o o o o o o o o o o \
    -o output/bigwig_hg19/multiBigwigSummary_NPC_Ciceri_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_hg19/multiBigwigSummary_NPC_Ciceri.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_hg19/multiBigwigSummary_NPC_Ciceri_heatmap.pdf





