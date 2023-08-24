#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

x=("CT14_Het" "CT42_Het" "CT43_Het" "CT20_KO" "CT38_KO" "CT41_KO")


for x in "${x[@]}"; do
    bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done


