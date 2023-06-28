#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.14-GCC-11.2.0
module load picard/2.26.10-Java-15

x=("8wN_WT_IGG_R1" "8wN_WT_H3K27me3_R1"
   "8wN_WT_IGG_R2" "8wN_WT_H3K27me3_R2" 
   "8wN_WT_IGG_R3" "8wN_WT_H3K27me3_R3"
   "8wN_WT_IGG_R4" "8wN_WT_H3K27me3_R4"
   "8wN_HET_IGG_R1" "8wN_HET_H3K27me3_R1"
   "8wN_HET_IGG_R2" "8wN_HET_H3K27me3_R2")


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



