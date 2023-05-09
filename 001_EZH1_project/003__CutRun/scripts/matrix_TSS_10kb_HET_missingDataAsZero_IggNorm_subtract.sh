#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_subtract/8wN_HET_H3K27me3_R1_subtract.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm/8wN_HET_H3K27me3_R2_subtract.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm/8wN_HET_H3K27me3_R3_subtract.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm/8wN_HET_H3K27me3_R4_subtract.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_subtract.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_subtract.bed




