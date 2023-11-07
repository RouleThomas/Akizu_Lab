#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=("PSC_KO_H3K27me3"
"PSC_KO_EZH2"
"PSC_KO_SUZ12"
"PSC_KO_IGG"
"PSC_KO_EZH1cs"
"PSC_KO_HA")


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
done


