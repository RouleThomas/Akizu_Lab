#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00




computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R ../003__CutRun/meta/ENCFF159KBI_peak_noIntergenic.gtf \
    -S output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_KO_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_KO_H3K27me3_R2.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_KO_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_KO_H3K27me3_R2.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/2dN_KO_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/2dN_KO_H3K27me3_R2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_KO_Rep.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_KO_Rep.bed



plotProfile -m output/deeptools/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_KO_Rep.gz \
    -out output/deeptools/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_KO_Rep_profile.png \
    --perGroup \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



