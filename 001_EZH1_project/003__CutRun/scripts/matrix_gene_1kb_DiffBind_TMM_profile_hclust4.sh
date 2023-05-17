#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_gene_1kb_DiffBind_TMM.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_profile_hclust4.png \
    --perGroup \
    --hclust 4 \
    --colors black blue red \
    --plotTitle "" \
    -T "H3K27me3 read density" \
    -z ""



