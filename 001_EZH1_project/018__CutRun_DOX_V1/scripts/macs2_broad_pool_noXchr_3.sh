#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/ESC_OEKO_H3K27me3_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_H3K27me3_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_H3K27me3_R3_noXchr.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3_noXchr.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_OEKO_H3K27me3_noXchr_pool --broad 



macs2 callpeak -t output/bowtie2/ESC_OEKO_EZH2_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_EZH2_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_EZH2_R3_noXchr.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3_noXchr.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_OEKO_EZH2_noXchr_pool --broad 


macs2 callpeak -t output/bowtie2/ESC_OEKO_EZH1_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_EZH1_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_OEKO_EZH1_R3_noXchr.unique.dupmark.sorted.bam \
    -c output/bowtie2/ESC_WT_IGG_R1_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R2_noXchr.unique.dupmark.sorted.bam \
    output/bowtie2/ESC_WT_IGG_R3_noXchr.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_OEKO_EZH1_noXchr_pool --broad 


