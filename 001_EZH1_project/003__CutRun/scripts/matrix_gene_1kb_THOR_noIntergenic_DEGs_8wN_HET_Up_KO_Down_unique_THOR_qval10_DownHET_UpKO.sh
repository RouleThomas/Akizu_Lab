#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00




computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_peak_noIntergenic_DEGs_8wN_HET_Up_KO_Down_unique_THOR_qval10_DownHET_UpKO_sort_bedsort.gtf \
    -S output/THOR/THOR_WTvsHET/WTvsHET-s1_median.bw output/THOR/THOR_WTvsHET/WTvsHET-s2_median.bw output/THOR/THOR_WTvsKO/WTvsKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_noIntergenic_DEGs_8wN_HET_Up_KO_Down_unique_THOR_qval10_DownHET_UpKO.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_noIntergenic_DEGs_8wN_HET_Up_KO_Down_unique_THOR_qval10_DownHET_UpKO.bed


plotProfile -m output/deeptools/matrix_gene_1kb_noIntergenic_DEGs_8wN_HET_Up_KO_Down_unique_THOR_qval10_DownHET_UpKO.gz \
    -out output/deeptools/matrix_gene_1kb_noIntergenic_DEGs_8wN_HET_Up_KO_Down_unique_THOR_qval10_DownHET_UpKO_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

