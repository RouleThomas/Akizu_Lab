#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_NotGenotypeGroup_IggNorm_subtract/8wN_WT_H3K27me3_R1_subtract.bw output/bigwig_histone_NotGenotypeGroup_IggNorm_subtract/8wN_WT_H3K27me3_R2_subtract.bw output/bigwig_histone_NotGenotypeGroup_IggNorm_subtract/8wN_WT_H3K27me3_R3_subtract.bw output/bigwig_histone_NotGenotypeGroup_IggNorm_subtract/8wN_WT_H3K27me3_R4_subtract.bw \
    --skipZeros \
    --missingDataAsZero \
    --minThreshold 2 \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_WT_corr_IggNorm_subtract_min2.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_WT_corr_IggNorm_subtract_min2.bed


