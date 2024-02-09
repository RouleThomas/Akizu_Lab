#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dNFA_KOEF1aEZH1_EZH1cs"
"50dNnative_KOEF1aEZH1_EZH1cs"
"50dNFA_KOEF1aEZH1_EZH2"
"50dNnative_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_H3K27me3"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




