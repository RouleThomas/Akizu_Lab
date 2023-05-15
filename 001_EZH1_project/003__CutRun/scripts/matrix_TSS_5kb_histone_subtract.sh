#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_IggNorm_subtract/8wN_WT_H3K27me3_median.bw output/bigwig_histone_IggNorm_subtract/8wN_HET_H3K27me3_median.bw output/bigwig_histone_IggNorm_subtract/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_histone_subtract.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_histone_subtract.bed



