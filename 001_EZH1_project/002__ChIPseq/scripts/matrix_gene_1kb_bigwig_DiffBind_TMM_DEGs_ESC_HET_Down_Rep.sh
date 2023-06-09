#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R /scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_HET_Down.gtf \
    -S output/bigwig_DiffBind_TMM/ESC_HET_H3K27me3_R1.bw output/bigwig_DiffBind_TMM/ESC_HET_H3K27me3_R2.bw output/bigwig_DiffBind_TMM/ESC_KO_H3K27me3_R1.bw output/bigwig_DiffBind_TMM/ESC_KO_H3K27me3_R2.bw output/bigwig_DiffBind_TMM/ESC_WT_H3K27me3_R1.bw output/bigwig_DiffBind_TMM/ESC_WT_H3K27me3_R2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_Down_Rep.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_Down_Rep.bed



plotProfile -m output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_Down_Rep.gz \
    -out output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_Down_Rep_profile.png \
    --perGroup \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



