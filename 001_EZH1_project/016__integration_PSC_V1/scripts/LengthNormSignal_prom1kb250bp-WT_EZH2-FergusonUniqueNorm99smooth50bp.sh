#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw \
    -R ../../Master/meta/ENCFF159KBI_geneSymbol_prom1kb250bp.bed \
    --outFileName output/edgeR/LengthNormSignal_prom1kb250bp-WT_EZH2-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_prom1kb250bp-WT_EZH2-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_prom1kb250bp-WT_EZH2-FergusonUniqueNorm99smooth50bp.bed













