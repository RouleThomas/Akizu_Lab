#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_WT_H3K27me3_ratio_median.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_HET_H3K27me3_ratio_median.bw output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio/8wN_KO_H3K27me3_ratio_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio.bed


