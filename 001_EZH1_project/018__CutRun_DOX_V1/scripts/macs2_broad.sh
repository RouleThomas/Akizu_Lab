#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"




samples_and_controls=(
    # WT samples
    "ESC_WT_EZH1_R1" "ESC_WT_IGG_R1"
    "ESC_WT_EZH1_R2" "ESC_WT_IGG_R2"
    "ESC_WT_EZH1_R3" "ESC_WT_IGG_R3"
    "ESC_WT_EZH2_R1" "ESC_WT_IGG_R1"
    "ESC_WT_EZH2_R2" "ESC_WT_IGG_R2"
    "ESC_WT_EZH2_R3" "ESC_WT_IGG_R3"
    "ESC_WT_H3K27me3_R1" "ESC_WT_IGG_R1"
    "ESC_WT_H3K27me3_R2" "ESC_WT_IGG_R2"
    "ESC_WT_H3K27me3_R3" "ESC_WT_IGG_R3"

    # KO samples (control = WT IGG same replicate)
    "ESC_KO_EZH1_R1" "ESC_WT_IGG_R1"
    "ESC_KO_EZH1_R2" "ESC_WT_IGG_R2"
    "ESC_KO_EZH1_R3" "ESC_WT_IGG_R3"
    "ESC_KO_EZH2_R1" "ESC_WT_IGG_R1"
    "ESC_KO_EZH2_R2" "ESC_WT_IGG_R2"
    "ESC_KO_EZH2_R3" "ESC_WT_IGG_R3"
    "ESC_KO_H3K27me3_R1" "ESC_WT_IGG_R1"
    "ESC_KO_H3K27me3_R2" "ESC_WT_IGG_R2"
    "ESC_KO_H3K27me3_R3" "ESC_WT_IGG_R3"

    # OEKO samples (control = WT IGG same replicate)
    "ESC_OEKO_EZH1_R1" "ESC_WT_IGG_R1"
    "ESC_OEKO_EZH1_R2" "ESC_WT_IGG_R2"
    "ESC_OEKO_EZH1_R3" "ESC_WT_IGG_R3"
    "ESC_OEKO_EZH2_R1" "ESC_WT_IGG_R1"
    "ESC_OEKO_EZH2_R2" "ESC_WT_IGG_R2"
    "ESC_OEKO_EZH2_R3" "ESC_WT_IGG_R3"
    "ESC_OEKO_H3K27me3_R1" "ESC_WT_IGG_R1"
    "ESC_OEKO_H3K27me3_R2" "ESC_WT_IGG_R2"
    "ESC_OEKO_H3K27me3_R3" "ESC_WT_IGG_R3"
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







