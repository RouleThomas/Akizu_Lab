#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "CPC_untreated_YAP1_R1" "CPC_untreated_input_R3"
    "CPC_untreated_YAP1_R2" "CPC_untreated_input_R3"
    "CPC_untreated_TEAD4_R1" "CPC_untreated_input_R3"
    "CPC_untreated_TEAD4_R2" "CPC_untreated_input_R3"
    "CPC_untreated_NR2F2_R1" "CPC_untreated_input_R3"
    "CPC_untreated_NR2F2_R2" "CPC_untreated_input_R3"
    "CPC_untreated_NR2F2_R3" "CPC_untreated_input_R3"
    "CPC_untreated_YAP1_R3" "CPC_untreated_input_R3"
    "CPC_untreated_TEAD4_R3" "CPC_untreated_input_R3"

    "CPC_RA_YAP1_R1" "CPC_RA_input_R3"
    "CPC_RA_YAP1_R2" "CPC_RA_input_R3"
    "CPC_RA_TEAD4_R1" "CPC_RA_input_R3"
    "CPC_RA_TEAD4_R2" "CPC_RA_input_R3"
    "CPC_RA_NR2F2_R1" "CPC_RA_input_R3"
    "CPC_RA_NR2F2_R2" "CPC_RA_input_R3"
    "CPC_RA_NR2F2_R3" "CPC_RA_input_R3"
    "CPC_RA_YAP1_R3" "CPC_RA_input_R3"
    "CPC_RA_TEAD4_R3" "CPC_RA_input_R3"
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


## run in interactive

  
macs2 callpeak -t output/bowtie2/CPC_RA_YAP1_R3.unique.dupmark.sorted.bam \
    -c output/bowtie2/CPC_RA_input_R3.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir output/macs2/narrow -n CPC_RA_YAP1_R3

macs2 callpeak -t output/bowtie2/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bam \
    -c output/bowtie2/CPC_RA_input_R3.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir output/macs2/narrow -n CPC_RA_TEAD4_R3

