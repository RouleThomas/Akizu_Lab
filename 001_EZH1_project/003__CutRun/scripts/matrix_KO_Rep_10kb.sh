#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
-b 10000 -a 10000 \
-R ../../Master/meta/ENCFF159KBI.gtf \
-S output/bigwig_histone_NotGenotypeGroup_lib/8wN_KO_H3K27me3_R1.bw output/bigwig_histone_NotGenotypeGroup_lib/8wN_KO_H3K27me3_R2.bw output/bigwig_histone_NotGenotypeGroup_lib/8wN_KO_H3K27me3_R3.bw output/bigwig_histone_NotGenotypeGroup_lib/8wN_KO_H3K27me3_R4.bw \
--skipZeros \
--blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
-o output/deeptools/matrix_KO_Rep_10kb.gz \
-p 6 \
--outFileSortedRegions output/deeptools/matrix_KO_Rep_10kb.bed