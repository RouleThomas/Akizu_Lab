#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary BED-file -b output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s1-rep0.bw output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s1-rep1.bw output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s1-rep2.bw output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s2-rep0.bw output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s2-rep1.bw output/THOR/THOR_PSC_WTvsKO_EZH2_TMM/PSCWTvsKOEZH2TMM-s2-rep2.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM/PSCWTvsKOEF1aEZH1EZH2TMM-s2-rep0.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM/PSCWTvsKOEF1aEZH1EZH2TMM-s2-rep1.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM/PSCWTvsKOEF1aEZH1EZH2TMM-s2-rep2.bw -o output/deeptools/multiBigwigSummary_EZH2_THOR_TMM_bed.npz --BED /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI_gene.bed






# Plot
## PCA
plotPCA -in output/deeptools/multiBigwigSummary_EZH2_THOR_TMM_bed.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --colors black black black red red red blue blue blue \
    -o output/deeptools/multiBigwigSummary_EZH2_THOR_TMM_bed_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/deeptools/multiBigwigSummary_EZH2_THOR_TMM_bed.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/deeptools/multiBigwigSummary_EZH2_THOR_TMM_bed_heatmap.pdf





