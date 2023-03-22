#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL



x=("2dN_HET_H3K27me3_R1"
"2dN_HET_H3K27me3_R2"
"2dN_HET_input_R1"
"2dN_HET_input_R2"
"2dN_KO_H3K27me3_R1"
"2dN_KO_H3K27me3_R2"
"2dN_KO_input_R1"
"2dN_KO_input_R2"
"2dN_WT_H3K27me3_R1"
"2dN_WT_H3K27me3_R2"
"2dN_WT_input_R1"
"2dN_WT_input_R2")

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




