#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=72:00:00

x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")



for x in "${x[@]}"; do
    bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done

