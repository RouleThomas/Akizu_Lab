#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"hESC_WT_DVL2_R1"
"hESC_YAPKO_DVL2_R1"
"hESC_WT_EZH2_R1"
"hESC_YAPKO_EZH2_R1"
"hESC_WT_EZH2_R2"
"hESC_YAPKO_EZH2_R2"
"hESC_WT_QSER1_R1"
"hESC_YAPKO_QSER1_R1"
"hESC_WT_QSER1_R2"
"hESC_YAPKO_QSER1_R2"
"hESC_WT_input_R1"
)

for x in "${x[@]}"; do
    fastp -i input/${x}.fq.gz  \
    -o output/fastp/${x}.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




