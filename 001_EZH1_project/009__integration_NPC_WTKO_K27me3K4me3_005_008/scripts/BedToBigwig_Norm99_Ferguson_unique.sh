#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "NPC_WT_H3K4me3_008"
    "NPC_WT_H3K4me3_005"
    "NPC_KO_H3K4me3_008"
    "NPC_KO_H3K4me3_005"
    "NPC_WT_H3K27me3_008"
    "NPC_WT_H3K27me3_005"
    "NPC_KO_H3K27me3_008"
    "NPC_KO_H3K27me3_005"
)




for sample in "${input_list[@]}"; do
    bedtools sort -i output/bigwig_Ferguson/${sample}_unique_norm99.bedGraph > output/bigwig_Ferguson/${sample}_unique_norm99.sorted.bedGraph

    bedGraphToBigWig \
        output/bigwig_Ferguson/${sample}_unique_norm99.sorted.bedGraph \
        ../../Master/meta/GRCh38_chrom_sizes.tab \
        output/bigwig_Ferguson/${sample}_unique_norm99.bw
done






