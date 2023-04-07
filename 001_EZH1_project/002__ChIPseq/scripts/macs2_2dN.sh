#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
  "2dN_HET_H3K27me3_R1" "2dN_HET_input_R1"
  "2dN_HET_H3K27me3_R2" "2dN_HET_input_R2"
  "2dN_KO_H3K27me3_R1" "2dN_KO_input_R1"
  "2dN_KO_H3K27me3_R2" "2dN_KO_input_R2"
  "2dN_WT_H3K27me3_R1" "2dN_WT_input_R1"
  "2dN_WT_H3K27me3_R2" "2dN_WT_input_R2"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
        -c output/bowtie2_endtoend/${control}.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample} --broad 
done
