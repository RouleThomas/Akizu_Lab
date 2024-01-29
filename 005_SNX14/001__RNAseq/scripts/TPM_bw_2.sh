#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00

x=("171HetCB" "174MTCB" "175HetCB" "177MTCB" "474WTCB" "171HetCX" "174MTCX" "175HetCX" "177MTCX" "474WTCX")


for x in "${x[@]}"; do
    bamCoverage --bam output/bam/1year_AL1804271_R2_new_analysis/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done


