#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dN_KO_input_R1"
"50dN_KOEF1aEZH1_input_R1"
"50dN_WT_input_R1"
"50dN_WT_EZH1_R1"
"50dN_WT_EZH1_R2"
"50dN_KO_EZH2_R1"
"50dN_KO_EZH2_R2"
"50dN_KOEF1aEZH1_EZH2_R1"
"50dN_KOEF1aEZH1_EZH2_R2"
"50dN_WT_EZH2_R1"
"50dN_WT_EZH2_R2"
"50dN_WT_H3K27me1_R1"
"50dN_WT_H3K27me1_R2"
)

for x in "${x[@]}"; do
    fastp -i input_raw/${x}.fq.gz  \
    -o output/fastp/${x}.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




