#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_KO_H3K27me3_R1_log2ratio.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_KO_H3K27me3_R2_log2ratio.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_KO_H3K27me3_R3_log2ratio.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_KO_H3K27me3_R4_log2ratio.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_log2ratio.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_log2ratio.bed




