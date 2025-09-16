#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
  "ESC_WT_H3K27me3_R1" "ESC_WT_H3K27me3_R2" "ESC_WT_H3K27me3_R3"
  "ESC_KO_H3K27me3_R1" "ESC_KO_H3K27me3_R2" "ESC_KO_H3K27me3_R3"
  "ESC_OEKO_H3K27me3_R1" "ESC_OEKO_H3K27me3_R2" "ESC_OEKO_H3K27me3_R3"
)




for sample in "${input_list[@]}"; do
    bedtools sort -i output/bigwig_Ferguson/${sample}_noXchr_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/${sample}_noXchr_unique_norm99_initialBigwig.sorted.bedGraph

    bedGraphToBigWig \
        output/bigwig_Ferguson/${sample}_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
        ../../Master/meta/GRCh38_chrom_sizes.tab \
        output/bigwig_Ferguson/${sample}_noXchr_unique_norm99_initialBigwig.bw
done






