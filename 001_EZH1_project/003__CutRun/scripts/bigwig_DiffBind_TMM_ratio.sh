#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


samples_and_controls=(
"8wN_WT_H3K27me3_R1" "8wN_WT_IGG_R1"
"8wN_WT_H3K27me3_R2" "8wN_WT_IGG_R2"
"8wN_WT_H3K27me3_R3" "8wN_WT_IGG_R3"
"8wN_WT_H3K27me3_R4" "8wN_WT_IGG_R4"

"8wN_KO_H3K27me3_R1" "8wN_KO_IGG_R1"
"8wN_KO_H3K27me3_R2" "8wN_KO_IGG_R2"
"8wN_KO_H3K27me3_R3" "8wN_KO_IGG_R3"
"8wN_KO_H3K27me3_R4" "8wN_KO_IGG_R4"

"8wN_HET_H3K27me3_R1" "8wN_HET_IGG_R1"
"8wN_HET_H3K27me3_R2" "8wN_HET_IGG_R2"
"8wN_HET_H3K27me3_R3" "8wN_HET_IGG_R3"
"8wN_HET_H3K27me3_R4" "8wN_HET_IGG_R4"

"8wN_iPSCpatient_H3K27me3_R1" "8wN_iPSCpatient_IGG_R1"
"8wN_iPSCpatient_H3K27me3_R2" "8wN_iPSCpatient_IGG_R2"
)


for ((i=0; i<${#samples_and_controls[@]}; i+=2)); do
  sample="${samples_and_controls[i]}"
  control="${samples_and_controls[i+1]}"
  
    bigwigCompare -b1 output/bigwig_DiffBind_TMM/${sample}.dupmark.sorted.bw \
        -b2 output/bigwig_DiffBind_TMM/${control}.dupmark.sorted.bw \
        --operation ratio \
        -o output/bigwig_DiffBind_TMM_ratio/${sample}_ratio.bw
done


