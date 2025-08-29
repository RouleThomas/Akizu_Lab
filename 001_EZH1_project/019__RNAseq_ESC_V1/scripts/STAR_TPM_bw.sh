#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00

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
    bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig_STAR/${x}.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 1
done


