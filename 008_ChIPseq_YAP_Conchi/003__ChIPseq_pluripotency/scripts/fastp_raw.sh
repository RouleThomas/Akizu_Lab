#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"hESC_WT_TEAD4_R1"
"hESC_WT_YAP1_R1"
"hESC_WT_input_R1"
)

for x in "${x[@]}"; do
    fastp -i input/${x}.fq.gz  \
    -o output/fastp/${x}.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




