#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"50dN_KO_input_R1"
"50dN_KOEF1aEZH1_input_R1"
"50dN_WT_input_R1"
"50dN_WT_EZH1_R1"
"50dN_WT_EZH1_R2"
"50dN_KO_EZH2_R1"
"50dN_KO_EZH2_R2"
"50dN_KOEF1aEZH1_EZH2_R1"
"50dN_KOEF1aEZH1_EZH2_R2"
"50dN_WT_EZH2_R1"
"50dN_WT_EZH2_R2"
"50dN_WT_H3K27me1_R1"
"50dN_WT_H3K27me1_R2"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -U output/fastp/${x}.fq.gz
done


