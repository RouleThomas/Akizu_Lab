#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp${x}_2.fq.gz
done



