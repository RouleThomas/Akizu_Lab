#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=5

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("NPC_KO_EZH1cs"
"NPC_KO_EZH1pt"
"NPC_KO_EZH2"
"NPC_KO_H3K27me1"
"NPC_KO_H3K27me3"
"NPC_KO_H3K4me3"
"NPC_KO_IGG"
"NPC_KO_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH1pt")


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../003__CutRun/meta/MG1655 \
            -S output/spikein/${x}_MG1655.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
    samtools view -S -F 4 -c output/spikein/${x}_MG1655.sam > output/spikein/${x}_MG1655_count.txt
done


