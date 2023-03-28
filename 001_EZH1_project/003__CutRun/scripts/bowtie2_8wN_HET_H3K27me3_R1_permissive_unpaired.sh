#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --time=100:00:00

nthreads=7

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

x=("8wN_HET_H3K27me3_R1")


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --local --no-unal --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/tmp/${x}_permissive_unpaired.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
done


