#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99.bw \
    -R ../../Master/meta/ENCFF159KBI_geneSymbol_prom2500bp2500bp.bed \
    --outFileName output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_KOEF1aEZH1_H3K27me3_005R-FergusonUniqueNorm99smooth50bp.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 5000 \
    --regionBodyLength 5000




