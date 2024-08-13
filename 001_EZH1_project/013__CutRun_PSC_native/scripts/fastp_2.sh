#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"PSC_WT_IGG_R1"
"PSC_WT_SUZ12_R1"
"PSC_KO_EZH1_R2"
"PSC_KO_EZH2_R2"
"PSC_KO_H3K27me3_R2"
"PSC_KO_IGG_R2"
"PSC_KO_SUZ12_R2"
"PSC_KOEF1aEZH1_EZH1_R2"
"PSC_KOEF1aEZH1_EZH2_R2"
"PSC_KOEF1aEZH1_H3K27me3_R2"
"PSC_KOEF1aEZH1_IGG_R2"
"PSC_KOEF1aEZH1_SUZ12_R2"
"PSC_WTEF1aEZH1_EZH1_R2"
"PSC_WTEF1aEZH1_EZH2_R2"
"PSC_WTEF1aEZH1_H3K27me3_R2"
"PSC_WTEF1aEZH1_IGG_R2"
"PSC_WTEF1aEZH1_SUZ12_R2"
"PSC_WT_EZH1_R2"
"PSC_WT_EZH2_R2"
"PSC_WT_H3K27me3_R2"
"PSC_WT_IGG_R2"
"PSC_WT_SUZ12_R2"
)

for x in "${x[@]}"; do
    fastp -i input_raw/${x}_1.fq.gz -I input_raw/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




