#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("PSC_KOEF1aEZH1_IGG"
"PSC_KOEF1aEZH1_SUZ12"
"PSC_KOsynEZH1_EZH1cs"
"PSC_KOsynEZH1_EZH1pt"
"PSC_KOsynEZH1_H3K27me3"
"PSC_KOsynEZH1_HA"
"PSC_KOsynEZH1_IGG"
"PSC_KOsynEZH1_SUZ12")


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



