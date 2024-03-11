#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/DMSO_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/DMSO_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/DMSO_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/DMSO_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/DMSO_noAB_R1.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_noAB_R1.unique.dupmark.sorted.bw output/bigwig/DOT1Linh_noAB_R2.unique.dupmark.sorted.bw output/bigwig/EHMTinh_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/EHMTinh_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/EHMTinh_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/EHMTinh_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/EHMTinh_noAB_R1.unique.dupmark.sorted.bw output/bigwig/EHMTinh_noAB_R2.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/EZH2inh_noAB_R1.unique.dupmark.sorted.bw output/bigwig/EZH2inh_noAB_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_all.npz




# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels DMSO_H3K27me3_R1 DMSO_H3K27me3_R2 DMSO_H3K4me3_R1 DMSO_H3K4me3_R2 DMSO_noAB_R1 DOT1Linh_H3K27me3_R1 DOT1Linh_H3K27me3_R2 DOT1Linh_H3K4me3_R1 DOT1Linh_H3K4me3_R2 DOT1Linh_noAB_R1 DOT1Linh_noAB_R2 EHMTinh_H3K27me3_R1 EHMTinh_H3K27me3_R2 EHMTinh_H3K4me3_R1 EHMTinh_H3K4me3_R2 EHMTinh_noAB_R1 EHMTinh_noAB_R2 EZH2inh_H3K27me3_R1 EZH2inh_H3K27me3_R2 EZH2inh_H3K4me3_R1 EZH2inh_H3K4me3_R2 EZH2inh_noAB_R1 EZH2inh_noAB_R2 \
    --colors "#1e90ff" "#1e90ff" "#ff6347" "#ff6347" "#8b4513" "#32cd32" "#32cd32" "#db7093" "#db7093" "#cd853f" "#cd853f" "#3cb371" "#3cb371" "#b22222" "#b22222" "#a0522d" "#a0522d" "#6b8e23" "#6b8e23" "#800000" "#800000" "#daa520" "#daa520" \
    --markers o o o o o o o o o o o o o o o o o o o o o \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels DMSO_H3K27me3_R1 DMSO_H3K27me3_R2 DMSO_H3K4me3_R1 DMSO_H3K4me3_R2 DMSO_noAB_R1 DOT1Linh_H3K27me3_R1 DOT1Linh_H3K27me3_R2 DOT1Linh_H3K4me3_R1 DOT1Linh_H3K4me3_R2 DOT1Linh_noAB_R1 DOT1Linh_noAB_R2 EHMTinh_H3K27me3_R1 EHMTinh_H3K27me3_R2 EHMTinh_H3K4me3_R1 EHMTinh_H3K4me3_R2 EHMTinh_noAB_R1 EHMTinh_noAB_R2 EZH2inh_H3K27me3_R1 EZH2inh_H3K27me3_R2 EZH2inh_H3K4me3_R1 EZH2inh_H3K4me3_R2 EZH2inh_noAB_R1 EZH2inh_noAB_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf





