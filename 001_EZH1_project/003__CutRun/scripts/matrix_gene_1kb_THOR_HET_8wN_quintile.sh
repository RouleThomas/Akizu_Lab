#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_quintile_8wN_HET_notExpress.gtf meta/ENCFF159KBI_quintile_8wN_HET_quint1.gtf meta/ENCFF159KBI_quintile_8wN_HET_quint2.gtf meta/ENCFF159KBI_quintile_8wN_HET_quint3.gtf meta/ENCFF159KBI_quintile_8wN_HET_quint4.gtf \
    -S output/THOR/THOR_WTvsHET/WTvsHET-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_THOR_HET_8wN_quintile.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_THOR_HET_8wN_quintile.bed


plotProfile -m output/deeptools/matrix_gene_1kb_THOR_HET_8wN_quintile.gz \
    -out output/deeptools/matrix_gene_1kb_THOR_HET_8wN_quintile_profile.pdf




