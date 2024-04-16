#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"hESC_WT_DVL2_R1"
"hESC_YAPKO_DVL2_R1"
"hESC_WT_EZH2_R1"
"hESC_YAPKO_EZH2_R1"
"hESC_WT_EZH2_R2"
"hESC_YAPKO_EZH2_R2"
"hESC_WT_QSER1_R1"
"hESC_YAPKO_QSER1_R1"
"hESC_WT_QSER1_R2"
"hESC_YAPKO_QSER1_R2"
"hESC_WT_input_R1"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -U output/fastp/${x}.fq.gz
done


