#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



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
  
    bamCompare -b1 output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
        -b2 output/bowtie2_endtoend/${control}.dupmark.sorted.bam \
        -o output/bigwig_inputNorm/${sample}_ratio.bw
done


