#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


macs2 callpeak -t output/bowtie2/NPC_WT_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_H3K27me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/NPC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_WT_H3K27me3_pool --broad 

macs2 callpeak -t output/bowtie2/NPC_WT_H3K4me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_H3K4me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/NPC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_WT_H3K4me3_pool --broad 

macs2 callpeak -t output/bowtie2/NPC_WT_H3K9me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_H3K9me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/NPC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_WT_H3K9me3_pool --broad 

macs2 callpeak -t output/bowtie2/NPC_WT_H3K27ac_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_H3K27ac_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/NPC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/NPC_WT_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_WT_H3K27ac_pool --broad 



