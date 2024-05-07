#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/hESC_WT_DNMT3A_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_DNMT3A_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_DNMT3A_pool --broad 

macs2 callpeak -t output/bowtie2/hESC_WT_DNMT3B_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_DNMT3B_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_DNMT3B_pool --broad 

macs2 callpeak -t output/bowtie2/hESC_WT_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_H3K27me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2/hESC_WT_H3K4me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_H3K4me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_H3K4me3_pool --broad 


macs2 callpeak -t output/bowtie2/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_QSER1FLAG_pool --broad 




