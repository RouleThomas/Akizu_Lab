#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00




samples_and_scaling_factors=(
  "ESC_HET_H3K27me3_R1"	10.51
  "ESC_HET_H3K27me3_R2"	23.35
  "ESC_KO_H3K27me3_R1"	10.06
  "ESC_KO_H3K27me3_R2"	15.78
  "ESC_WT_H3K27me3_R1"	7
  "ESC_WT_H3K27me3_R2"	4.31
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



