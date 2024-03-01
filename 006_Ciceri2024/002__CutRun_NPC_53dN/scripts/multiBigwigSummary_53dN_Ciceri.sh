#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig_hg19/GSM5860153_CnR_H3K27ac_D53_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860154_CnR_H3K27ac_D53_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860157_CnR_H3K27me3_D53_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860158_CnR_H3K27me3_D53_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860161_CnR_H3K4me3_D53_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860162_CnR_H3K4me3_D53_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860165_CnR_H3K9me3_D53_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860166_CnR_H3K9me3_D53_rep2.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860169_CnR_IgG_D53_rep1.hg19.sorted.RmDup.10mNorm.bw output/bigwig_hg19/GSM5860170_CnR_IgG_D53_rep2.hg19.sorted.RmDup.10mNorm.bw -o output/bigwig_hg19/multiBigwigSummary_53dN_Ciceri.npz



# Plot
## PCA
plotPCA -in output/bigwig_hg19/multiBigwigSummary_53dN_Ciceri.npz \
    --transpose \
    --ntop 0 \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --colors "#00bfff" "#87cefa" "#8b0000" "#ff0000" "#006400" "#32cd32" "#9400d3" "#9370db" "#8b4513" "#cd853f" \
    --markers o o o o o o o o o o \
    -o output/bigwig_hg19/multiBigwigSummary_53dN_Ciceri_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_hg19/multiBigwigSummary_53dN_Ciceri.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels H3K27ac_R1 H3K27ac_R2 H3K27me3_R1 H3K27me3_R2 H3K4me3_R1 H3K4me3_R2 H3K9me3_R1 H3K9me3_R2 IGG_R1 IGG_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_hg19/multiBigwigSummary_53dN_Ciceri_heatmap.pdf





