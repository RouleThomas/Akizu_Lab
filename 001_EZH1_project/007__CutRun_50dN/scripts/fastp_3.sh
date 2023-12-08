#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dN_WTQ731E_EZH1cs_R1"
"50dN_WTQ731E_EZH1cs_R2"
"50dN_WTQ731E_EZH2_R1"
"50dN_WTQ731E_EZH2_R2"
"50dN_WTQ731E_H3K27me3_R3"
"50dN_WTQ731E_H3K27me3_R1"
"50dN_WTQ731E_H3K27me3_R2"
"50dN_WTQ731E_IGG_R1"
"50dN_WTQ731E_IGG_R2"
"50dN_WTQ731E_SUZ12_R1"
"50dN_WTQ731E_SUZ12_R2"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




