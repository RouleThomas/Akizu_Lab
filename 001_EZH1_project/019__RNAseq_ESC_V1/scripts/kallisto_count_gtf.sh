#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


x=(
    "ESC_WT_R1"
    "ESC_KO_R1"
    "ESC_OEKO_R1"
    "ESC_WT_R2"
    "ESC_KO_R2"
    "ESC_OEKO_R2"
    "ESC_WT_R3"
    "ESC_KO_R3"
    "ESC_OEKO_R3"
    )
        
for x in "${x[@]}"; do
echo "Processing sample ${x}"
kallisto quant -i ../../Master/meta/kallisto/transcripts.idx \
    -o output/kallisto/${x}_quant \
    -b 100 output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
    -g ../../Master/meta/gencode.v47.annotation.gtf -t 8 --rf-stranded --genomebam --chromosomes ../../Master/meta/GRCh38_chrom_sizes.tab
done







