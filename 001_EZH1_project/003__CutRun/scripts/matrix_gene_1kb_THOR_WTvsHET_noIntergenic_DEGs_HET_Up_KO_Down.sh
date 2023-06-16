#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_sort.gtf \
    -S output/THOR/THOR_WTvsHET/WTvsHET-s1_median.bw output/THOR/THOR_WTvsHET/WTvsHET-s2_median.bw output/THOR/THOR_WTvsKO/WTvsKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_THOR_WTvsHET_noIntergenic_DEGs_HET_Up_KO_Down.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_THOR_WTvsHET_noIntergenic_DEGs_HET_Up_KO_Down.bed


plotProfile -m output/deeptools/matrix_gene_1kb_THOR_WTvsHET_noIntergenic_DEGs_HET_Up_KO_Down.gz \
    -out output/deeptools/matrix_gene_1kb_THOR_WTvsHET_noIntergenic_DEGs_HET_Up_KO_Down_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

