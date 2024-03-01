#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bw  ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_IGG_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_IGG_008.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_NPC_CutRun001008_Ciceri.npz


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_NPC_CutRun001008_Ciceri.npz \
    --transpose \
    --ntop 0 \
    --labels Akizu_H3K27me3_R1 Akizu_H3K27me3_R2 Akizu_H3K4me3_R1 Akizu_H3K4me3_R2 Akizu_IGG_R1 Akizu_IGG_R2 Ciceri_H3K27me3_R1 Ciceri_H3K27me3_R2 Ciceri_H3K4me3_R1 Ciceri_H3K4me3_R2 Ciceri_IGG_R1 Ciceri_IGG_R2 \
    --colors "#8b0000" "#ff0000" "#006400" "#32cd32" "#8b4513" "#cd853f" "#8b0000" "#ff0000" "#006400" "#32cd32" "#8b4513" "#cd853f" \
    --markers x x x x x x o o o o o o \
    -o output/bigwig/multiBigwigSummary_NPC_CutRun001008_Ciceri_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_NPC_CutRun001008_Ciceri.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels Akizu_H3K27me3_R1 Akizu_H3K27me3_R2 Akizu_H3K4me3_R1 Akizu_H3K4me3_R2 Akizu_IGG_R1 Akizu_IGG_R2 Ciceri_H3K27me3_R1 Ciceri_H3K27me3_R2 Ciceri_H3K4me3_R1 Ciceri_H3K4me3_R2 Ciceri_IGG_R1 Ciceri_IGG_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_NPC_CutRun001008_Ciceri_heatmap.pdf





