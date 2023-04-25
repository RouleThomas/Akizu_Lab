#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00

nthreads=10

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("ESC_KO_H3K27me3_R1"
"ESC_KO_H3K27me3_R2"
"ESC_KO_input_R1"
"ESC_KO_input_R2"
"ESC_WT_H3K27me3_R1"
"ESC_WT_H3K27me3_R2"
"ESC_WT_H3K27me3_R3"
"ESC_WT_input_R1"
"ESC_WT_input_R2"
"ESC_WT_input_R3"
"NPC_HET_H3K27me3_R1"
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
    # Remove duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/bowtie2_endtoend/${x}.filter.bam \
        -O output/bowtie2_endtoend/${x}.unique.dupmark.bam \
        -M output/bowtie2_endtoend/${x}.unique.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES true \
	-ASSUME_SORTED true
    # sort reads after removing the duplicates
    samtools sort -o output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam \
        output/bowtie2_endtoend/${x}.unique.dupmark.bam
    # index the sorted reads
    samtools index output/bowtie2_endtoend/${x}.unique.dupmark.sorted.bam
done



