#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


x=(
    "PSC_KO_R1"
    "PSC_KO_R2"
    "PSC_KO_R3"
    )
        
for x in "${x[@]}"; do
echo "Processing sample ${x}"
kallisto quant -i ../../Master/meta/kallisto/transcripts.idx \
    -o output/kallisto/${x}_quant \
    -b 100 output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
    -g ../../Master/meta/gencode.v47.annotation.gtf -t 8
done







