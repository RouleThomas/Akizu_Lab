#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_peak_noIntergenic_DEGs_8wN_KO_Down.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_KO_Down.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_KO_Down.bed


plotProfile -m output/deeptools/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_KO_Down.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_KO_Down_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

