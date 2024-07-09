#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"hESC_WT_input_R4"
"hESC_WT_input_R5"
"hESC_YAPKO_input_R4"
"hESC_WT_QSER1_R4"
"hESC_WT_QSER1_R5"
"hESC_YAPKO_QSER1_R4"
"hESC_YAPKO_QSER1_R5"
"hESC_WT_TEAD4_R4"
"hESC_WT_TEAD4_R5"
"hESC_WT_YAP1_R4"
"hESC_WT_YAP1_R5"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input_raw/${x}.fq.gz
done