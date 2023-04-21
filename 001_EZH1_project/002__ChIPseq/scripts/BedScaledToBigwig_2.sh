#!/bin/bash
#SBATCH --mem=75G
#SBATCH --time=100:00:00




samples_and_scaling_factors=(
  "ESC_HET_H3K27me3_R1"	8.91
  "ESC_HET_H3K27me3_R2"	21.08
  "ESC_KO_H3K27me3_R1"	8.54
  "ESC_KO_H3K27me3_R2"	13.67
  "ESC_WT_H3K27me3_R1"	1.74
  "ESC_WT_H3K27me3_R2"	4.1
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



