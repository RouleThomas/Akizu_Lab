#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00

nthreads=5

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"PSC_KOEF1aEZH1_H3K27me3"
"PSC_KOEF1aEZH1_HA"
"PSC_KOEF1aEZH1_EZH1cs"
"PSC_KOEF1aEZH1_EZH2"
"PSC_KOEF1aEZH1_SUZ12"
"PSC_KOEF1aEZH1_IGG"
"PSC_KO_H3K27me3"
"PSC_KO_EZH2"
"PSC_KO_SUZ12"
"PSC_KO_IGG"
"PSC_KO_EZH1cs"
"PSC_KO_HA"
"PSC_WT_EZH2"
"PSC_WT_IGG"
"PSC_WT_H3K27me3"
"PSC_WT_EZH1cs"
"PSC_WT_HA"
"PSC_WT_SUZ12"
)


for x in "${x[@]}"; do
    # count the uniq map reads
    samtools view -S -F 4 -c output/spikein/${x}_MG1655.unique.dupmark.sorted.bam > output/spikein/${x}_MG1655.unique.dupmark.sorted-count.txt
done



