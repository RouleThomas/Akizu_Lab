#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00




x=(
    "PSC_WT_R1"
    "PSC_WT_R2"
    "PSC_WT_R3"
    "PSC_KO_R1"
    "PSC_KO_R2"
    "PSC_KO_R3"
    "PSC_KOEF1aEZH1_R1"
    "PSC_KOEF1aEZH1_R2"
    "PSC_KOEF1aEZH1_R3"
    )

for x in "${x[@]}"; do
    fastp -i input_raw/${x}_1.fq.gz -I input_raw/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done