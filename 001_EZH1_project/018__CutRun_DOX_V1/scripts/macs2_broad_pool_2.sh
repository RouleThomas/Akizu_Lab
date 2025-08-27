#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/ESC_KO_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_H3K27me3_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_H3K27me3_R3.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_KO_H3K27me3_pool --broad 



macs2 callpeak -t output/bowtie2/ESC_KO_EZH2_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_EZH2_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_EZH2_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_KO_EZH2_pool --broad 


macs2 callpeak -t output/bowtie2/ESC_KO_EZH1_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_EZH1_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_KO_EZH1_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_KO_EZH1_pool --broad 


