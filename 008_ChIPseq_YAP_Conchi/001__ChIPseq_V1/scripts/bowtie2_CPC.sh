#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"CPC_untreated_YAP1_R1"
"CPC_RA_YAP1_R1"
"CPC_untreated_YAP1_R2"
"CPC_RA_YAP1_R2"
"CPC_untreated_TEAD4_R1"
"CPC_RA_TEAD4_R1"
"CPC_untreated_TEAD4_R2"
"CPC_RA_TEAD4_R2"
"CPC_untreated_NR2F2_R1"
"CPC_RA_NR2F2_R1"
"CPC_untreated_NR2F2_R2"
"CPC_RA_NR2F2_R2"
"CPC_RA_NR2F2_R3"
"CPC_untreated_NR2F2_R3"
"CPC_untreated_input_R3"
"CPC_RA_input_R3"
"CPC_untreated_YAP1_R3"
"CPC_RA_YAP1_R3"
"CPC_untreated_TEAD4_R3"
"CPC_RA_TEAD4_R3"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -U output/fastp/${x}.fq.gz
done


