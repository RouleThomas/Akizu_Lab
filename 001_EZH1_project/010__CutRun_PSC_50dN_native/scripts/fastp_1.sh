#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"PSC_WT_EZH2"
"PSC_WT_H3K27me1"
"PSC_WT_H3K27me3"
"PSC_WT_IGG"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




