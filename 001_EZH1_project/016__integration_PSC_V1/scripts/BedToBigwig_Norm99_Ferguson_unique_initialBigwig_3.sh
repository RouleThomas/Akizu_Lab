#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "PSC_WT_EZH2_006R"
    "PSC_WT_EZH2_010R"
    "PSC_WT_EZH2_014R1"
    "PSC_WT_H3K27me3_006R"
    "PSC_WT_H3K27me3_010R"
    "PSC_WT_H3K27me3_013R1"
    "PSC_WT_SUZ12_006R"
    "PSC_WT_SUZ12_013R1"
    "PSC_WT_SUZ12_014R1"

)




for sample in "${input_list[@]}"; do
    bedtools sort -i output/bigwig_Ferguson/${sample}_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/${sample}_unique_norm99_initialBigwig.sorted.bedGraph

    bedGraphToBigWig \
        output/bigwig_Ferguson/${sample}_unique_norm99_initialBigwig.sorted.bedGraph \
        ../../Master/meta/GRCh38_chrom_sizes.tab \
        output/bigwig_Ferguson/${sample}_unique_norm99_initialBigwig.bw
done






