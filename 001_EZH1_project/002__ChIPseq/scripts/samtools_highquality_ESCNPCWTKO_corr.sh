#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("ESC_KO_H3K27me3_R1"
"ESC_KO_H3K27me3_R2"

"ESC_WT_H3K27me3_R1"
"ESC_WT_H3K27me3_R2"
"ESC_WT_H3K27me3_R3"

"NPC_KO_H3K27me3_R1"
"NPC_KO_H3K27me3_R2"

"NPC_WT_H3K27me3_R1"
"NPC_WT_H3K27me3_R2"
)


for x in "${x[@]}"; do
    # sort reads after marking the duplicates
    samtools sort -o output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam \
        output/bowtie2_endtoend/${x}.MAPQ30.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam
done



