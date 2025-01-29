#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"




macs2 callpeak -t output/bowtie2/PSC_WT_SUZ12_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_SUZ12_013R1.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_SUZ12_014R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_WT_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_013R1.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_014R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_WT_SUZ12_pool --broad 



macs2 callpeak -t output/bowtie2/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_WT_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_010R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_013R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_WT_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2/PSC_WT_EZH2_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_EZH2_010R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_EZH2_014R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_WT_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_010R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_WT_IGG_014R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_WT_EZH2_pool --broad 






