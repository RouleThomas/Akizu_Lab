#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "PSC_KO_EZH1_R1" "PSC_KO_IGG_R1"
    "PSC_KO_EZH1_R2" "PSC_KO_IGG_R2"
    "PSC_KO_EZH2_R1" "PSC_KO_IGG_R1"
    "PSC_KO_EZH2_R2" "PSC_KO_IGG_R2"
    "PSC_KO_H3K27me3_R1" "PSC_KO_IGG_R1"
    "PSC_KO_H3K27me3_R2" "PSC_KO_IGG_R2"
    "PSC_KO_SUZ12_R1" "PSC_KO_IGG_R1"
    "PSC_KO_SUZ12_R2" "PSC_KO_IGG_R2"

    "PSC_WT_EZH1_R1" "PSC_WT_IGG_R1"
    "PSC_WT_EZH1_R2" "PSC_WT_IGG_R2"
    "PSC_WT_EZH2_R1" "PSC_WT_IGG_R1"
    "PSC_WT_EZH2_R2" "PSC_WT_IGG_R2"
    "PSC_WT_H3K27me3_R1" "PSC_WT_IGG_R1"
    "PSC_WT_H3K27me3_R2" "PSC_WT_IGG_R2"
    "PSC_WT_SUZ12_R1" "PSC_WT_IGG_R1"
    "PSC_WT_SUZ12_R2" "PSC_WT_IGG_R2"


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







