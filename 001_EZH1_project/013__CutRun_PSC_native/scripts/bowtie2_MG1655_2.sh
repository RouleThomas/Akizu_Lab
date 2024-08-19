#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=5

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15


x=(
"PSC_WT_IGG_R1"
"PSC_WT_SUZ12_R1"
"PSC_KO_EZH1_R2"
"PSC_KO_EZH2_R2"
"PSC_KO_H3K27me3_R2"
"PSC_KO_IGG_R2"
"PSC_KO_SUZ12_R2"
"PSC_KOEF1aEZH1_EZH1_R2"
"PSC_KOEF1aEZH1_EZH2_R2"
"PSC_KOEF1aEZH1_H3K27me3_R2"
"PSC_KOEF1aEZH1_IGG_R2"
"PSC_KOEF1aEZH1_SUZ12_R2"
"PSC_WTEF1aEZH1_EZH1_R2"
"PSC_WTEF1aEZH1_EZH2_R2"
"PSC_WTEF1aEZH1_H3K27me3_R2"
"PSC_WTEF1aEZH1_IGG_R2"
"PSC_WTEF1aEZH1_SUZ12_R2"
"PSC_WT_EZH1_R2"
"PSC_WT_EZH2_R2"
"PSC_WT_H3K27me3_R2"
"PSC_WT_IGG_R2"
"PSC_WT_SUZ12_R2"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../003__CutRun/meta/MG1655 \
            -S output/spikein/${x}_MG1655.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
    samtools view -S -F 4 -c output/spikein/${x}_MG1655.sam > output/spikein/${x}_MG1655_count.txt
done


