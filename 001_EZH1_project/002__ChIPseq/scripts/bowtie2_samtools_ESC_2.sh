#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=150G
#SBATCH --time=100:00:00

nthreads=10

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("ESC_KO_input_R2"
"ESC_WT_H3K27me3_R1"
"ESC_WT_H3K27me3_R2"
"ESC_WT_H3K27me3_R3"
"ESC_WT_input_R1"
"ESC_WT_input_R2"
"ESC_WT_input_R3")


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2_endtoend/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
    samtools sort -o output/bowtie2_endtoend/${x}.bam \
        output/bowtie2_endtoend/${x}.sam
    # index the bam file
    samtools index output/bowtie2_endtoend/${x}.bam
    # remove reads without MAPQ>=20 and unmapped reads, secondary alignments, and reads failing quality check
    samtools view -@ ${nthreads} -F 772 -q 20 \
        -b output/bowtie2_endtoend/${x}.bam \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | \
        samtools sort -o output/bowtie2_endtoend/${x}.filter.bam
    # index filtered reads
    samtools index output/bowtie2_endtoend/${x}.filter.bam
    # mark duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/bowtie2_endtoend/${x}.filter.bam \
        -O output/bowtie2_endtoend/${x}.dupmark.bam \
        -M output/bowtie2_endtoend/${x}.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES false \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o output/bowtie2_endtoend/${x}.dupmark.sorted.bam \
        output/bowtie2_endtoend/${x}.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2_endtoend/${x}.dupmark.sorted.bam
done




