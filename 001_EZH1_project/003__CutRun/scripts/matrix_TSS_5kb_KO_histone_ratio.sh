#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_IggNorm/8wN_KO_H3K27me3_R1_ratio.bw output/bigwig_histone_IggNorm/8wN_KO_H3K27me3_R2_ratio.bw output/bigwig_histone_IggNorm/8wN_KO_H3K27me3_R3_ratio.bw output/bigwig_histone_IggNorm/8wN_KO_H3K27me3_R4_ratio.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_KO_histone_ratio.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_KO_histone_ratio.bed


