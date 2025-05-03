#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6





computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/ENCFF159KBI_PSC_WTKOEF1aEZH1_H3K27me3_bin1000space100_gt_pval05_padj001__GainLost_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12



plotHeatmap -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 12 12

plotHeatmap -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_heatmap3.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 50 50

plotProfile -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colors black blue \
    --perGroup \
    --plotWidth 4 \
    --yMax 6




plotProfile -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_plotProfile4.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colors black blue \
    --perGroup \
    --plotWidth 8 

plotProfile -m output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene.gz \
    -out output/deeptools/matrix_TSSTES_250bp100bp-DIFFREPS_GainLostMix-PSC_WTKOEF1aEZH1_H3K27me3-bin1000space100_gt_pval05_padj001-initialBigwig_H3K27me3-gene_plotProfile5.pdf \
    --samplesLabel "WT_H3K27me3" "KOEF1aEZH1_H3K27me3" \
    --colors black blue \
    --perGroup \
    --plotWidth 4 \
    --yMax 20
