#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00


input_list=("8wN_HET_IGG_R1"
"8wN_HET_H3K27me3_R1"
"8wN_HET_IGG_R2"
"8wN_HET_H3K27me3_R2"
"8wN_KO_IGG_R1"
"8wN_KO_H3K27me3_R1"
"8wN_KO_IGG_R2"
"8wN_KO_H3K27me3_R2"
"8wN_iPSCpatient_IGG_R1"
"8wN_iPSCpatient_H3K27me3_R1"
"8wN_WT_IGG_R1"
"8wN_WT_H3K27me3_R1"
"8wN_WT_IGG_R2"
"8wN_WT_H3K27me3_R2"
"8wN_HET_IGG_R3"
"8wN_HET_H3K27me3_R3"
"8wN_HET_IGG_R4"
"8wN_HET_H3K27me3_R4"
"8wN_KO_IGG_R3"
"8wN_KO_H3K27me3_R3"
"8wN_KO_IGG_R4"
"8wN_KO_H3K27me3_R4"
"8wN_iPSCpatient_IGG_R2"
"8wN_iPSCpatient_H3K27me3_R2"
"8wN_WT_IGG_R3"
"8wN_WT_H3K27me3_R3"
"8wN_WT_IGG_R4"
"8wN_WT_H3K27me3_R4"
)

for x in "${input_list[@]}"; do
    bedtools bamtobed -i output/bowtie2/${x}.dupmark.sorted.bam > \
    output/bigwig_histone_NotGenotypeGroup_lib/${x}.bed
done



samples_and_scaling_factors=(
    "8wN_HET_H3K27me3_R1" 1.584793814
    "8wN_HET_H3K27me3_R2" 1.88492268
    "8wN_HET_H3K27me3_R3" 1.050902062
    "8wN_HET_H3K27me3_R4" 2.857345361
    "8wN_HET_IGG_R1" 1.384716733
    "8wN_HET_IGG_R2" 1.135704875
    "8wN_HET_IGG_R3" 1.052700922
    "8wN_HET_IGG_R4" 2.559947299
    "8wN_KO_H3K27me3_R1" 1.306572165
    "8wN_KO_H3K27me3_R2" 1
    "8wN_KO_H3K27me3_R3" 3.809536082
    "8wN_KO_H3K27me3_R4" 1.128994845
    "8wN_KO_IGG_R1" 1.376811594
    "8wN_KO_IGG_R2" 1
    "8wN_KO_IGG_R3" 1.425559947
    "8wN_KO_IGG_R4" 1.027667984
    "8wN_WT_H3K27me3_R1" 1.690850515
    "8wN_WT_H3K27me3_R2" 2.925515464
    "8wN_WT_H3K27me3_R3" 2.051546392
    "8wN_WT_H3K27me3_R4" 3.093685567
    "8wN_WT_IGG_R1" 2.187088274
    "8wN_WT_IGG_R2" 2.802371542
    "8wN_WT_IGG_R3" 1.200263505
    "8wN_WT_IGG_R4" 5.216073781
    "8wN_iPSCpatient_H3K27me3_R1" 2.268170103
    "8wN_iPSCpatient_H3K27me3_R2" 1.47564433
    "8wN_iPSCpatient_IGG_R1" 1.606060606
    "8wN_iPSCpatient_IGG_R2" 1.180500659
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"


libSize=$(cat output/bigwig_histone_NotGenotypeGroup_lib/${sample}.bed | wc -l)
scale=$(echo "15000000/($libSize*$SF)" | bc -l)

    genomeCoverageBed -bg -scale $scale -i output/bigwig_histone_NotGenotypeGroup_lib/${sample}.bed \
    -g ../../Master/meta/GRCh38_chrom_sizes.tab > \
    output/bigwig_histone_NotGenotypeGroup_lib/${sample}.bedGraph

    bedtools sort -i output/bigwig_histone_NotGenotypeGroup_lib/${sample}.bedGraph > \
    output/bigwig_histone_NotGenotypeGroup_lib/${sample}.sorted.bedGraph

    bedGraphToBigWig output/bigwig_histone_NotGenotypeGroup_lib/${sample}.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_NotGenotypeGroup_lib/${sample}.bw

done

