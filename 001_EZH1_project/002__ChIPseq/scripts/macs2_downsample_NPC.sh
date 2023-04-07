#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G


macs2_out="output/macs2_downsample"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


samples_and_controls=(
  "NPC_HET_H3K27me3_R1" "NPC_HET_input_R1"
  "NPC_HET_H3K27me3_R2" "NPC_HET_input_R2"
  "NPC_KO_H3K27me3_R1" "NPC_KO_input_R1"
  "NPC_KO_H3K27me3_R2" "NPC_KO_input_R2"
  "NPC_WT_H3K27me3_R1" "NPC_WT_input_R1"
  "NPC_WT_H3K27me3_R2" "NPC_WT_input_R2"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    macs2 callpeak -t output/bowtie2_endtoend/downsample/${sample}.dupmark.sorted.bam \
        -c output/bowtie2_endtoend/downsample/${control}.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto \
        --nomodel -g hs \
        --outdir ${macs2_out} -n ${sample} --broad 
done
