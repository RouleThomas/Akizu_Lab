#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_gene_1kb_DiffBind_TMM_subtract.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_subtract_profile.png \
    --perGroup \
    --colors black blue red \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" \
    -T "H3K27me3 read density" \
    -z ""



