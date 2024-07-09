#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"hESC_WT_input_R4"
"hESC_WT_input_R5"
"hESC_YAPKO_input_R4"
"hESC_WT_QSER1_R4"
"hESC_WT_QSER1_R5"
"hESC_YAPKO_QSER1_R4"
"hESC_YAPKO_QSER1_R5"
"hESC_WT_TEAD4_R4"
"hESC_WT_TEAD4_R5"
"hESC_WT_YAP1_R4"
"hESC_WT_YAP1_R5"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -U output/fastp/${x}.fq.gz
done


