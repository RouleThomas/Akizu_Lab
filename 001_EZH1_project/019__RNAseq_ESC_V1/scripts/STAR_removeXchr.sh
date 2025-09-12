#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6   # Number of CPU cores per task





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
    samtools view -b output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam $(seq -f "chr%g" 1 22) > output/STAR/fastp/${x}_AlignednoXchr.sortedByCoord.out.bam
    samtools index output/STAR/fastp/${x}_AlignednoXchr.sortedByCoord.out.bam
done



