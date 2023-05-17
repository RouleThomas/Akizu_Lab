#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_DiffBind_TMM_subtract/8wN_HET_H3K27me3_R1_subtract.bw output/bigwig_DiffBind_TMM_subtract/8wN_HET_H3K27me3_R2_subtract.bw output/bigwig_DiffBind_TMM_subtract/8wN_HET_H3K27me3_R3_subtract.bw output/bigwig_DiffBind_TMM_subtract/8wN_HET_H3K27me3_R4_subtract.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_HET_DiffBind_TMM_subtract.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_HET_DiffBind_TMM_subtract.bed


