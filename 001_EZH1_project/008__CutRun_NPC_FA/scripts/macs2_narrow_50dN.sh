#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "50dNFA_KOEF1aEZH1_EZH1cs" "50dNFA_KOEF1aEZH1_IGG"
    "50dNnative_KOEF1aEZH1_EZH1cs" "50dNnative_KOEF1aEZH1_IGG"
    "50dNFA_KOEF1aEZH1_EZH2" "50dNFA_KOEF1aEZH1_IGG"
    "50dNnative_KOEF1aEZH1_EZH2" "50dNnative_KOEF1aEZH1_IGG"
    "50dNFA_KOEF1aEZH1_H3K27me3" "50dNFA_KOEF1aEZH1_IGG"
    "50dNnative_KOEF1aEZH1_H3K27me3" "50dNnative_KOEF1aEZH1_IGG"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2/${sample}.unique.dupmark.sorted.bam \
        -c output/bowtie2/${control}.unique.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample} 
done


