#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "NPC_WT_EZH1cs" "NPC_WT_IGG"
    "NPC_WT_EZH1pt" "NPC_WT_IGG"
    "NPC_WT_EZH2" "NPC_WT_IGG"
    "NPC_WT_H3K27me1" "NPC_WT_IGG"
    "NPC_WT_H3K27me3" "NPC_WT_IGG"
    "NPC_WT_SUZ12" "NPC_WT_IGG"

    "NPC_KO_EZH1cs" "NPC_KO_IGG"
    "NPC_KO_EZH1pt" "NPC_KO_IGG"
    "NPC_KO_EZH2" "NPC_KO_IGG"
    "NPC_KO_H3K27me1" "NPC_KO_IGG"
    "NPC_KO_H3K27me3" "NPC_KO_IGG"
    "NPC_KO_SUZ12" "NPC_KO_IGG"
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





samples_and_controls=(
    "NPC_WT_H3K4me3" "NPC_WT_IGG"

    "NPC_KO_H3K4me3" "NPC_KO_IGG"
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





