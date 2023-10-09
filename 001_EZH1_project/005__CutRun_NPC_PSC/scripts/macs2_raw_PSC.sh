#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "PSC_KOEF1aEZH1_EZH1cs" "PSC_KOEF1aEZH1_IGG"
    "PSC_KOEF1aEZH1_EZH1pt" "PSC_KOEF1aEZH1_IGG"
    "PSC_KOEF1aEZH1_H3K27me3" "PSC_KOEF1aEZH1_IGG"
    "PSC_KOEF1aEZH1_HA" "PSC_KOEF1aEZH1_IGG"
    "PSC_KOEF1aEZH1_SUZ12" "PSC_KOEF1aEZH1_IGG"

    "PSC_KOsynEZH1_EZH1cs" "PSC_KOsynEZH1_IGG"
    "PSC_KOsynEZH1_EZH1pt" "PSC_KOsynEZH1_IGG"
    "PSC_KOsynEZH1_H3K27me3" "PSC_KOsynEZH1_IGG"
    "PSC_KOsynEZH1_HA" "PSC_KOsynEZH1_IGG"
    "PSC_KOsynEZH1_SUZ12" "PSC_KOsynEZH1_IGG"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2/${sample}.unique.dupmark.sorted.bam \
        -c output/bowtie2/${control}.unique.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample} --broad 
done


