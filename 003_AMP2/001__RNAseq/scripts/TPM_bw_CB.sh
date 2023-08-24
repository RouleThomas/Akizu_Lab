#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

x=("CB14_Het" "CB42_Het" "CB43_Het" "CB20_KO" "CB38_KO" "CB41_KO")


for x in "${x[@]}"; do
    bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done


