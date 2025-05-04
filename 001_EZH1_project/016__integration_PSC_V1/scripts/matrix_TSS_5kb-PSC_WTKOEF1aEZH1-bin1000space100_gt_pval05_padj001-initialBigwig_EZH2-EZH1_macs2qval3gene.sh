#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_PSC_KOEF1aEZH1_EZH1_qval3_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive
plotHeatmap -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 20

plotHeatmap -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 0.5

plotHeatmap -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_heatmap_colorSmall3.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 5 5 10 10 150 150
        



plotProfile -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_plotProfile1.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black blue \
    --perGroup


plotProfile -m output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene.gz \
    -out output/deeptools/matrix_TSS_5kb-PSC_WTKOEF1aEZH1-bin1000space100_gt_pval05_padj001-initialBigwig_EZH2-EZH1_macs2qval3gene_plotProfile2.pdf \
    --samplesLabel "WT_EZH2" "KOEF1aEZH1_EZH2" \
    --colors black blue \
    --perGroup \
    --plotHeight 10 \
    --plotWidth 7

