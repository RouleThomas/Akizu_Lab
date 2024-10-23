#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"PSC_KOEF1aEZH1_EZH2_R1"
"PSC_KOEF1aEZH1_IGG_R1"
"PSC_WT_EZH1_R1"
"PSC_WT_EZH2_R1"
"PSC_WT_IGG_R1"
"PSC_WT_SUZ12_R1"
"PSC_KO_EZH2_R1"
"PSC_KO_IGG_R1"
"PSC_KO_SUZ12_R1"
"PSC_KO_EZH1_R2"
"PSC_KO_EZH2_R2"
"PSC_KO_H3K27me3_R2"
"PSC_KO_IGG_R2"
"PSC_KO_SUZ12_R2"
"PSC_KO_IGG_R3"
"PSC_KO_SUZ12_R3"
)


for x in "${x[@]}"; do
    # sort read
    samtools sort -o output/bowtie2/${x}.bam \
        output/bowtie2/${x}.sam
    # index the bam file
    samtools index output/bowtie2/${x}.bam
    # remove reads without MAPQ>=20 and unmapped reads, secondary alignments, and reads failing quality check
    samtools view -@ ${nthreads} -F 772 -q 20 \
        -b output/bowtie2/${x}.bam \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | \
        samtools sort -o output/bowtie2/${x}.filter.bam
    # index filtered reads
    samtools index output/bowtie2/${x}.filter.bam
    
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



