#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=100:00:00




samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1" 2.27
  "2dN_HET_H3K27me3_R2" 2.03
  "2dN_KO_H3K27me3_R1" 1.69
  "2dN_KO_H3K27me3_R2" 1.14
  "2dN_WT_H3K27me3_R1" 1.4
  "2dN_WT_H3K27me3_R2" 1.84
  "ESC_HET_H3K27me3_R1" 10.12
  "ESC_HET_H3K27me3_R2" 23.95
  "ESC_KO_H3K27me3_R1" 9.7
  "ESC_KO_H3K27me3_R2" 15.53
  "ESC_WT_H3K27me3_R1" 1.98
  "ESC_WT_H3K27me3_R2" 4.65
  "ESC_WT_H3K27me3_R3" 1
  "NPC_HET_H3K27me3_R1" 1.24
  "NPC_HET_H3K27me3_R2" 1.61
  "NPC_KO_H3K27me3_R1" 1.7
  "NPC_KO_H3K27me3_R2" 3.14
  "NPC_WT_H3K27me3_R1" 1.7
  "NPC_WT_H3K27me3_R2" 1.66
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"


libSize=$(cat output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bed | wc -l)
scale=$(echo "15000000/($libSize*$SF)" | bc -l)

    genomeCoverageBed -bg -scale $scale -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bed \
    -g ../../Master/meta/GRCh38_chrom_sizes.tab > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bedGraph

    bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${sample}.bw

done



