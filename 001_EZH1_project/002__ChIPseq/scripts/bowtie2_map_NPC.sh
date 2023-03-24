#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G

nthreads=5

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("NPC_HET_H3K27me3_R1"
"NPC_HET_H3K27me3_R2"
"NPC_HET_input_R1"
"NPC_HET_input_R2"
"NPC_KO_H3K27me3_R1"
"NPC_KO_H3K27me3_R2"
"NPC_KO_input_R1"
"NPC_KO_input_R2"
"NPC_WT_H3K27me3_R1"
"NPC_WT_H3K27me3_R2"
"NPC_WT_input_R1"
"NPC_WT_input_R2")



for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --local --no-mixed --no-unal --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz | \ # here we stream the sam output directly to samtools to avoid creating the file
    samtools view -Sb - | \
    samtools sort -o output/bowtie2/${x}.bam
    # index the bam file
    samtools index output/bowtie2/${x}.bam
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -f 0x2 -q 20 \
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


