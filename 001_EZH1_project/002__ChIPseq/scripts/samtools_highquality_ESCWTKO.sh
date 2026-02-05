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
    # sort read
    samtools sort -o output/bowtie2_endtoend/${x}.bam \
        output/bowtie2_endtoend/${x}.sam
    # index the bam file
    samtools index output/bowtie2_endtoend/${x}.bam
    
    # remove reads without MAPQ>=30 and unmapped reads, secondary alignments, and reads failing quality check
   samtools view -@ ${nthreads} -F 772 -q 30 \
        -b output/bowtie2_endtoend/${x}.bam \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | \
        samtools sort -o output/bowtie2_endtoend/${x}.MAPQ30.filter.bam
    # index filtered reads
    samtools index output/bowtie2_endtoend/${x}.MAPQ30.filter.bam
    # mark duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/bowtie2_endtoend/${x}.MAPQ30.filter.bam \
        -O output/bowtie2_endtoend/${x}.MAPQ30.dupmark.bam \
        -M output/bowtie2_endtoend/${x}.MAPQ30.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES true \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam \
        output/bowtie2_endtoend/${x}.unique.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam
done



