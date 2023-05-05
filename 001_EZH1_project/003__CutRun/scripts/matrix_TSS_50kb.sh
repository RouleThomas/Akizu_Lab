#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
-b 50000 -a 50000 \
-R ../../Master/meta/ENCFF159KBI.gtf \
-S output/bigwig_histone_NotGenotypeGroup_lib/8wN_WT_H3K27me3_median.bw output/bigwig_histone_NotGenotypeGroup_lib/8wN_HET_H3K27me3_median.bw output/bigwig_histone_NotGenotypeGroup_lib/8wN_KO_H3K27me3_median.bw \
--skipZeros \
--blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
-o output/deeptools/matrix_TSS_50kb.gz \
-p 6 \
--outFileSortedRegions output/deeptools/matrix_TSS_50kb.bed