#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00




x=(
"50dN_WT_EZH1"
"50dN_WT_EZH2"
"50dN_WT_H3K27ac"
"50dN_WT_H3K27me1AM"
"50dN_WT_H3K27me1OR"
"50dN_WT_H3K27me3"
"50dN_WT_IGG"
"50dN_WT_SUZ12"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




