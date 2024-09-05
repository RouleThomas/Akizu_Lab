#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00

nthreads=5

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15


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
    # count read
    samtools view -F 4 -c output/spikein/${x}_MG1655.unique.dupmark.sorted.bam > output/spikein/${x}_MG1655_count_unique.txt
done



