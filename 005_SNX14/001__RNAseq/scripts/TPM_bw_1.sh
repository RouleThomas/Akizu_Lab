#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00

x=("S_CB_KO1" "S_CB_KO2" "S_CB_KO3" "S_CB_WT1" "S_CB_WT2" "S_CB_WT3" "S_CX_KO1" "S_CX_KO2" "S_CX_KO3" "S_CX_WT1" "S_CX_WT2" "S_CX_WT3")


for x in "${x[@]}"; do
    bamCoverage --bam output/bam/1month_rnaseqyz072420/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 10  
done


