#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2_unique"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "8wN_WT_H3K27me3_R1" "8wN_WT_IGG_R1"
    "8wN_WT_H3K27me3_R2" "8wN_WT_IGG_R2"
    "8wN_WT_H3K27me3_R3" "8wN_WT_IGG_R3"
    "8wN_WT_H3K27me3_R4" "8wN_WT_IGG_R4"
    "8wN_KO_H3K27me3_R1" "8wN_KO_IGG_R1"
    "8wN_KO_H3K27me3_R2" "8wN_KO_IGG_R2"
    "8wN_KO_H3K27me3_R3" "8wN_KO_IGG_R3"
    "8wN_KO_H3K27me3_R4" "8wN_KO_IGG_R4"
    "8wN_HET_H3K27me3_R1" "8wN_HET_IGG_R1" 
    "8wN_HET_H3K27me3_R2" "8wN_HET_IGG_R2"
    "8wN_HET_H3K27me3_R3" "8wN_HET_IGG_R3"
    "8wN_HET_H3K27me3_R4" "8wN_HET_IGG_R4"
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
