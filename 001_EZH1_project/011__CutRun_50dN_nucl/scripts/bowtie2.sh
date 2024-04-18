#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"50dN_WT_EZH1"
"50dN_WT_EZH2"
"50dN_WT_H3K27ac"
"50dN_WT_H3K27me1AM"
"50dN_WT_H3K27me1OR"
"50dN_WT_H3K27me3"
"50dN_WT_IGG"
"50dN_WT_SUZ12"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
done


