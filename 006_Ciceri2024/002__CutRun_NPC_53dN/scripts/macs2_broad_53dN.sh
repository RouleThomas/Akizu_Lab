#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "53dN_WT_H3K27me3_R1" "53dN_WT_IGG_R1"
    "53dN_WT_H3K27me3_R2" "53dN_WT_IGG_R2"
    "53dN_WT_H3K4me3_R1" "53dN_WT_IGG_R1"
    "53dN_WT_H3K4me3_R2" "53dN_WT_IGG_R2"
    "53dN_WT_H3K9me3_R1" "53dN_WT_IGG_R1"
    "53dN_WT_H3K9me3_R2" "53dN_WT_IGG_R2"
    "53dN_WT_H3K27ac_R1" "53dN_WT_IGG_R1"
    "53dN_WT_H3K27ac_R2" "53dN_WT_IGG_R2"
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


