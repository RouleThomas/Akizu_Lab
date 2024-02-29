#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb1kb_THOR_WT_allGenes.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb1kb_THOR_WT_allGenes.bed \
    --outFileNameMatrix output/deeptools/matrix_gene_1kb1kb_THOR_WT_allGenes.txt



