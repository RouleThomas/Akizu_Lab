#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00




computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R ../003__CutRun/meta/ENCFF159KBI_peak_noIntergenic.gtf \
    -S output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_R1.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_R2.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_R1.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_R2.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_R1.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_R2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_ESC_noIntergenic_Rep.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_ESC_noIntergenic_Rep.bed


plotProfile -m output/deeptools/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_ESC_noIntergenic_Rep.gz \
    -out output/deeptools/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_ESC_noIntergenic_Rep_profile.png \
    --perGroup \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
