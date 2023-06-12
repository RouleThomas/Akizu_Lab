#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00




computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf \
    -S output/bigwig_DiffBind_TMM/ESC_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/NPC_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/2dN_WT_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT.bed


plotProfile -m output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT.gz \
    -out output/deeptools/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT_profile.png \
    --perGroup \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


