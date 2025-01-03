#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "PSC_KOEF1aEZH1_IGG_005R"
    "PSC_KOEF1aEZH1_IGG_006R"
    "PSC_KOEF1aEZH1_IGG_013R1"
    "PSC_KOEF1aEZH1_IGG_014R1"
    "PSC_KO_IGG_006R"
    "PSC_KO_IGG_013R1"
    "PSC_KO_IGG_014R1"
    "PSC_KO_IGG_014R2"
    "PSC_WT_IGG_006R"
    "PSC_WT_IGG_010R"
    "PSC_WT_IGG_013R1"
    "PSC_WT_IGG_014R1"
)




for sample in "${input_list[@]}"; do
    bedtools sort -i output/bigwig_Ferguson/${sample}_norm.bedGraph > output/bigwig_Ferguson/${sample}_norm.sorted.bedGraph

    bedGraphToBigWig \
        output/bigwig_Ferguson/${sample}_norm.sorted.bedGraph \
        ../../Master/meta/GRCh38_chrom_sizes.tab \
        output/bigwig_Ferguson/${sample}_norm.bw
done






