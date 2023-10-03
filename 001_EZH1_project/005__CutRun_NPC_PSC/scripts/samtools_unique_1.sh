#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("NPC_KO_EZH1cs"
"NPC_KO_EZH1pt"
"NPC_KO_EZH2"
"NPC_KO_H3K27me1"
"NPC_KO_H3K27me3"
"NPC_KO_H3K4me3"
"NPC_KO_IGG"
"NPC_KO_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH1pt")


for x in "${x[@]}"; do
    # Remove duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/bowtie2/${x}.filter.bam \
        -O output/bowtie2/${x}.unique.dupmark.bam \
        -M output/bowtie2/${x}.unique.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES true \
	-ASSUME_SORTED true
    # sort reads after removing the duplicates
    samtools sort -o output/bowtie2/${x}.unique.dupmark.sorted.bam \
        output/bowtie2/${x}.unique.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2/${x}.unique.dupmark.sorted.bam
done



