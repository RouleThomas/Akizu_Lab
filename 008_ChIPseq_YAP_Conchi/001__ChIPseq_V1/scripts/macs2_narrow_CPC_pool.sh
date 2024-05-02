#!/bin/bash
#SBATCH --mem=250G


macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bam \
    output/bowtie2/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/CPC_untreated_input_R3.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n CPC_untreated_NR2F2_pool

macs2 callpeak -t output/bowtie2/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bam \
    output/bowtie2/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/CPC_RA_input_R3.unique.dupmark.sorted.bam \
    --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n CPC_RA_NR2F2_pool
