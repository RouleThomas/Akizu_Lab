#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


x=("HP14_Het" "HP42_Het" "HP43_Het" "HP20_KO" "HP38_KO" "HP41_KO")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done



