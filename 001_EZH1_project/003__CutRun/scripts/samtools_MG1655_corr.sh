#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G
#SBATCH --time=100:00:00

nthreads=5

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("8wN_KO_IGG_R1_MG1655" "8wN_KO_H3K27me3_R1_MG1655"
   "8wN_KO_IGG_R2_MG1655" "8wN_KO_H3K27me3_R2_MG1655"
   "8wN_KO_IGG_R3_MG1655" "8wN_KO_H3K27me3_R3_MG1655"
   "8wN_KO_IGG_R4_MG1655" "8wN_KO_H3K27me3_R4_MG1655"
   "8wN_WT_IGG_R1_MG1655" "8wN_WT_H3K27me3_R1_MG1655"
   "8wN_WT_IGG_R2_MG1655" "8wN_WT_H3K27me3_R2_MG1655" 
   "8wN_WT_IGG_R3_MG1655" "8wN_WT_H3K27me3_R3_MG1655"
   "8wN_WT_IGG_R4_MG1655" "8wN_WT_H3K27me3_R4_MG1655"
   "8wN_HET_IGG_R1_MG1655" "8wN_HET_H3K27me3_R1_MG1655"
   "8wN_HET_IGG_R2_MG1655" "8wN_HET_H3K27me3_R2_MG1655"
   "8wN_HET_IGG_R3_MG1655" "8wN_HET_H3K27me3_R3_MG1655"
   "8wN_HET_IGG_R4_MG1655" "8wN_HET_H3K27me3_R4_MG1655"
   "8wN_iPSCpatient_IGG_R1_MG1655" "8wN_iPSCpatient_H3K27me3_R1_MG1655"
   "8wN_iPSCpatient_IGG_R2_MG1655" "8wN_iPSCpatient_H3K27me3_R2_MG1655")


for x in "${x[@]}"; do
    # sort read
    samtools sort -o output/spikein/${x}.bam \
        output/spikein/${x}.sam
    # index the bam file
    samtools index output/spikein/${x}.bam
    # remove reads without MAPQ>=20 and unmapped reads, secondary alignments, and reads failing quality check
    samtools view -@ ${nthreads} -F 772 -q 20 \
        -b output/spikein/${x}.bam \
	U00096.3 | \
        samtools sort -o output/spikein/${x}.filter.bam
    # index filtered reads
    samtools index output/spikein/${x}.filter.bam
    # mark duplicates with picard
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -I output/spikein/${x}.filter.bam \
        -O output/spikein/${x}.dupmark.bam \
        -M output/spikein/${x}.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES false \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o output/spikein/${x}.dupmark.sorted.bam \
        output/spikein/${x}.dupmark.bam
    # index the sorted reads
    samtools index output/spikein/${x}.dupmark.sorted.bam
done

