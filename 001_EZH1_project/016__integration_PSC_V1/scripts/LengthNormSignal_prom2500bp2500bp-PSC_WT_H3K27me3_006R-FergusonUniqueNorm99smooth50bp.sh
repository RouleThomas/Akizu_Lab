#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bw \
    -R ../../Master/meta/ENCFF159KBI_geneSymbol_prom2500bp2500bp.bed \
    --outFileName output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_prom2500bp2500bp-PSC_WT_H3K27me3_006R-FergusonUniqueNorm99smooth50bp.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 5000 \
    --regionBodyLength 5000




