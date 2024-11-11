#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "PSC_KO_EZH2_013R1"	"PSC_KO_IGG_013R1"
    "PSC_KO_EZH2_014R1"	"PSC_KO_IGG_014R1"
    "PSC_KO_EZH2_014R2"	"PSC_KO_IGG_014R2"
    "PSC_KO_EZH1_006R"	"PSC_KO_IGG_006R"
    "PSC_KO_EZH1_013R1"	"PSC_KO_IGG_013R1"
    "PSC_KO_EZH1_014R2"	"PSC_KO_IGG_014R2"
    "PSC_KO_H3K27me3_006R"	"PSC_KO_IGG_006R"
    "PSC_KO_H3K27me3_013R1"	"PSC_KO_IGG_013R1"
    "PSC_KO_H3K27me3_014R2"	"PSC_KO_IGG_014R2"
    "PSC_KO_SUZ12_013R1"	"PSC_KO_IGG_013R1"
    "PSC_KO_SUZ12_014R1"	"PSC_KO_IGG_014R1"
    "PSC_KO_SUZ12_014R2"	"PSC_KO_IGG_014R2"
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







