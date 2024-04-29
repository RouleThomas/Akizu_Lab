#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00




x=(
    "hESC_WT_DNMT3A_R1"
    "hESC_WT_DNMT3A_R2"
    "hESC_WT_DNMT3B_R1"
    "hESC_WT_DNMT3B_R2"
    "hESC_WT_H3K27me3_R1"
    "hESC_WT_H3K27me3_R2"
    "hESC_WT_H3K4me3_R1"
    "hESC_WT_H3K4me3_R2"
    "hESC_WT_input_R1"
    "hESC_WT_input_R2"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done