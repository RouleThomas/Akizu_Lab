#!/bin/bash
#SBATCH --mem=75G
#SBATCH --time=100:00:00




samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1"	2
  "2dN_HET_H3K27me3_R2"	1.79
  "2dN_KO_H3K27me3_R1"	1.49
  "2dN_KO_H3K27me3_R2"	1
  "2dN_WT_H3K27me3_R1"	1.23
  "2dN_WT_H3K27me3_R2"	1.62
  "NPC_HET_H3K27me3_R1"	1.09
  "NPC_HET_H3K27me3_R2"	1.42
  "NPC_KO_H3K27me3_R1"	1.5
  "NPC_KO_H3K27me3_R2"	2.76
  "NPC_WT_H3K27me3_R1"	1.49
  "NPC_WT_H3K27me3_R2"	1.46
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"


libSize=$(cat output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bed | wc -l)
scale=$(echo "15000000/($libSize*$SF)" | bc -l)

    bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.sorted.bedGraph

    bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.sorted.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bw

done



