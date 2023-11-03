#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00

nthreads=5

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
    # sort read
    samtools sort -o output/spikein/${x}_MG1655.bam \
        output/spikein/${x}_MG1655.sam
    # index the bam file
    samtools index output/spikein/${x}_MG1655.bam
    # remove reads without MAPQ>=20 and unmapped reads, secondary alignments, and reads failing quality check
    samtools view -@ ${nthreads} -F 772 -q 20 \
        -b output/spikein/${x}_MG1655.bam \
	U00096.3 | \
        samtools sort -o output/spikein/${x}_MG1655.filter.bam
    # index filtered reads
    samtools index output/spikein/${x}_MG1655.filter.bam
    # Remove duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/spikein/${x}_MG1655.filter.bam \
        -O output/spikein/${x}_MG1655.unique.dupmark.bam \
        -M output/spikein/${x}_MG1655.unique.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES true \
	-ASSUME_SORTED true
    # sort reads after removing the duplicates
    samtools sort -o output/spikein/${x}_MG1655.unique.dupmark.sorted.bam \
        output/spikein/${x}_MG1655.unique.dupmark.bam
    # index the sorted reads
    samtools index output/spikein/${x}_MG1655.unique.dupmark.sorted.bam
done



