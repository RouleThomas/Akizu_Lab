#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


 


computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs_8wN_bedtools.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.bed


plotProfile -m output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.gz \
    -out output/deeptools/mmatrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""




plotHeatmap -m output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs_heatmap_kmeans6.png \
    --perGroup \
    --kmeans 6 \
    --colorMap bwr \
    --samplesLabel "WT" "HET" "KO" \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs_heatmap_kmeans6.txt
