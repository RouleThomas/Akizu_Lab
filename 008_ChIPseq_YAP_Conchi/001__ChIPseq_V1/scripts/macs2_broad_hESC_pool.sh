#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/hESC_WT_QSER1_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_QSER1_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_QSER1_pool --broad 

macs2 callpeak -t output/bowtie2/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_YAPKO_QSER1_pool --broad 


macs2 callpeak -t output/bowtie2/hESC_WT_EZH2_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_WT_EZH2_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_WT_EZH2_pool --broad 

macs2 callpeak -t output/bowtie2/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bam \
    output/bowtie2/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n hESC_YAPKO_EZH2_pool --broad 









