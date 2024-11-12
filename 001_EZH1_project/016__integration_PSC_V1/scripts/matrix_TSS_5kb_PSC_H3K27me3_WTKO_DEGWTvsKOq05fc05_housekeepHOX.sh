#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=100:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../015__RNAseq_PSC/meta/ENCFF159KBI_upregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf ../015__RNAseq_PSC/meta/ENCFF159KBI_downregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf \
    -S output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s1-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s1-rep1.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s1-rep2.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s2-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s2-rep1.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOX/PSCWTvsKOH3K27me3housekeepHOX-s2-rep2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_H3K27me3_WTKO_DEGWTvsKOq05fc05_housekeepHOX.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_H3K27me3_WTKO_DEGWTvsKOq05fc05_housekeepHOX.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_H3K27me3_WTKO_DEGWTvsKOq05fc05_housekeepHOX_heatmap_colorSmall.pdf \
    --samplesLabel "WT_H3K27me3_R1" "WT_H3K27me3_R2" "WT_H3K27me3_R3" "KO_H3K27me3_R1" "KO_H3K27me3_R2" "KO_H3K27me3_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


