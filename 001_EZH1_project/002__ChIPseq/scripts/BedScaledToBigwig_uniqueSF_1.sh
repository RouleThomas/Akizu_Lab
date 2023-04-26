#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00




samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1"	1.97
  "2dN_HET_H3K27me3_R2"	1.75
  "2dN_KO_H3K27me3_R1"	1.46
  "2dN_KO_H3K27me3_R2"	1
  "2dN_WT_H3K27me3_R1"	1.29
  "2dN_WT_H3K27me3_R2"	1.69
  "NPC_HET_H3K27me3_R1"	1.13
  "NPC_HET_H3K27me3_R2"	1.43
  "NPC_KO_H3K27me3_R1"	1.55
  "NPC_KO_H3K27me3_R2"	2.64
  "NPC_WT_H3K27me3_R1"	1.45
  "NPC_WT_H3K27me3_R2"	1.51
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"


libSize=$(cat output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bed | wc -l)
scale=$(echo "15000000/($libSize*$SF)" | bc -l)

    genomeCoverageBed -bg -scale $scale -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bed \
    -g ../../Master/meta/GRCh38_chrom_sizes.tab > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/${sample}.bedGraph

    bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/${sample}.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/${sample}.sorted.bedGraph

    bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/${sample}.sorted.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/${sample}.bw

done



