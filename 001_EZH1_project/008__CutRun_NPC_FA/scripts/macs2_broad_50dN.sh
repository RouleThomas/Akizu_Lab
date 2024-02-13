#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

x=(
"50dNFA_KOEF1aEZH1_EZH1cs"
"50dNnative_KOEF1aEZH1_EZH1cs"
"50dNFA_KOEF1aEZH1_EZH2"
"50dNnative_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_H3K27me3"
)


for x in "${x[@]}"; do
    macs2 callpeak -t output/bowtie2/${x}.unique.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${x} --broad 
done





