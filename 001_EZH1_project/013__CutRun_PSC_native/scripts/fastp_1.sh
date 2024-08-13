#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"NEU_WT_H3K27me3_R1"
"NEU_WT_IGG_R1"
"NEU_WT_SUZ12_R1"
"PSC_KO_EZH1_R1"
"PSC_KO_EZH2_R1"
"PSC_KO_H3K27me3_R1"
"PSC_KO_IGG_R1"
"PSC_KO_SUZ12_R1"
"PSC_KOEF1aEZH1_EZH1_R1"
"PSC_KOEF1aEZH1_EZH2_R1"
"PSC_KOEF1aEZH1_H3K27me3_R1"
"PSC_KOEF1aEZH1_IGG_R1"
"PSC_KOEF1aEZH1_SUZ12_R1"
"PSC_WTEF1aEZH1_EZH1_R1"
"PSC_WTEF1aEZH1_EZH2_R1"
"PSC_WTEF1aEZH1_H3K27me3_R1"
"PSC_WTEF1aEZH1_IGG_R1"
"PSC_WTEF1aEZH1_SUZ12_R1"
"PSC_WT_EZH1_R1"
"PSC_WT_EZH2_R1"
"PSC_WT_H3K27me3_R1"
)

for x in "${x[@]}"; do
    fastp -i input_raw/${x}_1.fq.gz -I input_raw/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




