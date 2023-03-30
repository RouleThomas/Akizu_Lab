#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G
#SBATCH --time=100:00:00

nthreads=5

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("8wN_KO_IGG_R1" "8wN_KO_H3K27me3_R1"
   "8wN_KO_IGG_R2" "8wN_KO_H3K27me3_R2"
   "8wN_KO_IGG_R3" "8wN_KO_H3K27me3_R3"
   "8wN_KO_IGG_R4" "8wN_KO_H3K27me3_R4")


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
    # mark duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/bowtie2/${x}.filter.bam \
        -O output/bowtie2/${x}.dupmark.bam \
        -M output/bowtie2/${x}.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES false \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o output/bowtie2/${x}.dupmark.sorted.bam \
        output/bowtie2/${x}.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2/${x}.dupmark.sorted.bam
done



