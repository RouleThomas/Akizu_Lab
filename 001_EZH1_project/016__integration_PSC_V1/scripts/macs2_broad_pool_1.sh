#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2/PSC_KOEF1aEZH1_SUZ12_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_SUZ12_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_SUZ12_013R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_KOEF1aEZH1_IGG_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_013R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_KOEF1aEZH1_SUZ12_pool --broad 




macs2 callpeak -t output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_KOEF1aEZH1_IGG_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_013R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_KOEF1aEZH1_H3K27me3_pool --broad 



macs2 callpeak -t output/bowtie2/PSC_KOEF1aEZH1_EZH2_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_EZH2_013R1.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_EZH2_014R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_KOEF1aEZH1_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_013R1.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_014R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_KOEF1aEZH1_EZH2_pool --broad 

macs2 callpeak -t output/bowtie2/PSC_KOEF1aEZH1_EZH1_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_EZH1_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_EZH1_013R1.unique.dupmark.sorted.bam \
    -c output/bowtie2/PSC_KOEF1aEZH1_IGG_005R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_006R.unique.dupmark.sorted.bam \
    output/bowtie2/PSC_KOEF1aEZH1_IGG_013R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n PSC_KOEF1aEZH1_EZH1_pool --broad 



