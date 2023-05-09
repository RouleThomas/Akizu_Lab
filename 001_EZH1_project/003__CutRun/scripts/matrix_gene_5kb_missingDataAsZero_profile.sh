#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_gene_5kb_missingDataAsZero.gz \
    -out output/deeptools/matrix_gene_5kb_missingDataAsZero_profile.png \
    --perGroup \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" \
    -T "H3K27me3 read density" \
    -z ""



