#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_EZH2_H3K27me3_noIntergenic_ENCODE_sort.gtf \
    -S ../003__CutRun/output/THOR/THOR_WTvsHET/WTvsHET-s1_median.bw ../003__CutRun/output/THOR/THOR_WTvsHET/WTvsHET-s2_median.bw ../003__CutRun/output/THOR/THOR_WTvsKO/WTvsKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_THOR_EZH2_H3K27me3_noIntergenic_ENCODE.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_THOR_EZH2_H3K27me3_noIntergenic_ENCODE.bed


plotProfile -m output/deeptools/matrix_gene_1kb_THOR_EZH2_H3K27me3_noIntergenic_ENCODE.gz \
    -out output/deeptools/matrix_gene_1kb_THOR_EZH2_H3K27me3_noIntergenic_ENCODE_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

