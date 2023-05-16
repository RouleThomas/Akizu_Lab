#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R3.dupmark.sorted.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_KO_DiffBind_TMM.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_KO_DiffBind_TMM.bed


