#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00




x=(
    "EZH2inh_H3K27me3_R1"
    "EHMTinh_noAB_R2"
    "EHMTinh_noAB_R1"
    "EHMTinh_H3K4me3_R2"
    "EHMTinh_H3K4me3_R1"
    "EHMTinh_H3K27me3_R2"
    "EHMTinh_H3K27me3_R1"
    "DOT1Linh_noAB_R2"
    "DOT1Linh_noAB_R1"
    "DOT1Linh_H3K4me3_R2"
    "DOT1Linh_H3K4me3_R1"
    "DOT1Linh_H3K27me3_R2"
    "DOT1Linh_H3K27me3_R1"
    "DMSO_noAB_R1"
    "DMSO_H3K4me3_R2"
    "DMSO_H3K4me3_R1"
    "EZH2inh_noAB_R2"
    "EZH2inh_noAB_R1"
    "EZH2inh_H3K4me3_R2"
    "EZH2inh_H3K4me3_R1"
    "EZH2inh_H3K27me3_R2"
    "DMSO_H3K27me3_R2"
    "DMSO_H3K27me3_R1"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done