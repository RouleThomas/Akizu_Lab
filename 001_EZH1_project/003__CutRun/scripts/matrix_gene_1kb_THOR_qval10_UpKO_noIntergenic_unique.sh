#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_WTvsKO_THOR_qval10_UpKO_noIntergenic.gtf \
    -S output/THOR/THOR_WTvsHET/WTvsHET-s1_median.bw output/THOR/THOR_WTvsHET/WTvsHET-s2_median.bw output/THOR/THOR_WTvsKO/WTvsKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_THOR_qval10_UpKO_noIntergenic.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_THOR_qval10_UpKO_noIntergenic.bed


plotProfile -m output/deeptools/matrix_gene_1kb_THOR_qval10_UpKO_noIntergenic.gz \
    -out output/deeptools/matrix_gene_1kb_THOR_qval10_UpKO_noIntergenic_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

