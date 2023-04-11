#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
    "8wN_iPSCpatient_H3K27me3_R1" "8wN_iPSCpatient_IGG_R1"
    "8wN_iPSCpatient_H3K27me3_R2" "8wN_iPSCpatient_IGG_R2"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2/${sample}.dupmark.sorted.bam \
        -c output/bowtie2/${control}.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample} --broad 
done
