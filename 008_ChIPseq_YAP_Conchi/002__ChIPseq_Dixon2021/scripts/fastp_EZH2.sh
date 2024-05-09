#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00




x=(
    "hESC_WT_EZH2_R1"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done