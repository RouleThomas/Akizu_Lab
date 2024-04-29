#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "hESC_WT_DVL2_R1" "hESC_WT_input_R1"
    "hESC_YAPKO_DVL2_R1" "hESC_WT_input_R1"
    "hESC_WT_EZH2_R1" "hESC_WT_input_R1"
    "hESC_YAPKO_EZH2_R1" "hESC_WT_input_R1"
    "hESC_WT_EZH2_R2" "hESC_WT_input_R1"
    "hESC_YAPKO_EZH2_R2" "hESC_WT_input_R1"
    "hESC_WT_QSER1_R1" "hESC_WT_input_R1"
    "hESC_YAPKO_QSER1_R1" "hESC_WT_input_R1"
    "hESC_WT_QSER1_R2" "hESC_WT_input_R1"
    "hESC_YAPKO_QSER1_R2" "hESC_WT_input_R1"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2/${sample}.unique.dupmark.sorted.bam \
        -c output/bowtie2/${control}.unique.dupmark.sorted.bam \
        --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample}
done







