#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


x=(
    "PSC_WT_R1"
    "PSC_WT_R2"
    "PSC_WT_R3"
    )
        
for x in "${x[@]}"; do
echo "Processing sample ${x}"
salmon quant -i ../../Master/meta/salmon/Homo_sapiens -l A \
         -1 output/fastp/${x}_1.fq.gz \
         -2 output/fastp/${x}_2.fq.gz \
         -p 8 --validateMappings -o salmon/${x}_quant
done







