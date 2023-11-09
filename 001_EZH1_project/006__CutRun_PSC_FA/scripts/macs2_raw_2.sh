#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "PSC_KO_H3K27me3" "PSC_KO_IGG"
    "PSC_KO_EZH2" "PSC_KO_IGG"
    "PSC_KO_SUZ12" "PSC_KO_IGG"
    "PSC_KO_EZH1cs" "PSC_KO_IGG"
    "PSC_KO_HA" "PSC_KO_IGG"
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







