#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00



x=(
    "DMSO_H3K27me3_R2"
    "DMSO_H3K27me3_R1"
    "EZH2inh_H3K27me3_R2"
    "EZH2inh_H3K4me3_R2"
    "EZH2inh_H3K4me3_R1"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done