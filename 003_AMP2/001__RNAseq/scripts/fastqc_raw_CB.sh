#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


x=("CB14_Het" "CB42_Het" "CB43_Het" "CB20_KO" "CB38_KO" "CB41_KO")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done



