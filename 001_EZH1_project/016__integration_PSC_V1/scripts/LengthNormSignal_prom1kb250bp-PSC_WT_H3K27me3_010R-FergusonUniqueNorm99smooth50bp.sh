#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bw \
    -R ../../Master/meta/ENCFF159KBI_geneSymbol_prom1kb250bp.bed \
    --outFileName output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_prom1kb250bp-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99smooth50bp.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 1250 \
    --regionBodyLength 1250




