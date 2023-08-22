#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


x=("CT14_Het" "CT42_Het" "CT43_Het" "CT20_KO" "CT38_KO" "CT41_KO")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done



