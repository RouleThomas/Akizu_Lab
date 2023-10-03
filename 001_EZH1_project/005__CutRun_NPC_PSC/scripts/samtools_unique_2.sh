#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("NPC_WT_EZH2"
"NPC_WT_H3K27me1"
"NPC_WT_H3K27me3"
"NPC_WT_H3K4me3"
"NPC_WT_IGG"
"NPC_WT_SUZ12"
"PSC_KOEF1aEZH1_EZH1cs"
"PSC_KOEF1aEZH1_EZH1pt"
"PSC_KOEF1aEZH1_H3K27me3"
"PSC_KOEF1aEZH1_HA")


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



