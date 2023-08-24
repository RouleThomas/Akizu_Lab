#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

x=("HP14_Het" "HP42_Het" "HP43_Het" "HP20_KO" "HP38_KO" "HP41_KO")


for x in "${x[@]}"; do
    bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done


