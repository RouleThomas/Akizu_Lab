#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_NotGenotypeGroup/8wN_WT_H3K27me3_median.bw output/bigwig_histone_NotGenotypeGroup/8wN_HET_H3K27me3_median.bw output/bigwig_histone_NotGenotypeGroup/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --minThreshold 2 \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_corr_min2.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_corr_min2.bed


