#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



# Indicate here sample names
x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




