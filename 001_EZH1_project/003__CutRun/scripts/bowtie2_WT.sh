#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=150G
#SBATCH --time=100:00:00

nthreads=10

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("8wN_WT_IGG_R1" "8wN_WT_H3K27me3_R1"
   "8wN_WT_IGG_R2" "8wN_WT_H3K27me3_R2" 
   "8wN_WT_IGG_R3" "8wN_WT_H3K27me3_R3"
   "8wN_WT_IGG_R4" "8wN_WT_H3K27me3_R4")


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
done


