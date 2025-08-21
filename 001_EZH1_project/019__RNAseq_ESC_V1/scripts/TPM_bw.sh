#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00

x=(
    "PSC_WT_R1"
    "PSC_WT_R2"
    "PSC_WT_R3"
    "PSC_KO_R1"
    "PSC_KO_R2"
    "PSC_KO_R3"
    "PSC_KOEF1aEZH1_R1"
    "PSC_KOEF1aEZH1_R2"
    "PSC_KOEF1aEZH1_R3"
    )
        
for x in "${x[@]}"; do
    bamCoverage --bam output/kallisto/${x}_quant/pseudoalignments.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 1
done


