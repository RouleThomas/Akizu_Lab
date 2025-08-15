#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00

nthreads=10

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
"ESC_KO_EZH1_R1"
"ESC_KO_EZH2_R1"
"ESC_KO_H3K27me3_R1"
"ESC_OEKO_EZH1_R1"
"ESC_OEKO_EZH2_R1"
"ESC_OEKO_H3K27me3_R1"
"ESC_WT_EZH1_R1"
"ESC_WT_EZH2_R1"
"ESC_WT_H3K27me3_R1"
"ESC_WT_IGG_R1"
"ESC_KO_EZH1_R2"
"ESC_KO_EZH2_R2"
"ESC_KO_H3K27me3_R2"
"ESC_OEKO_EZH1_R2"
"ESC_OEKO_EZH2_R2"
"ESC_OEKO_H3K27me3_R2"
"ESC_WT_EZH1_R2"
"ESC_WT_EZH2_R2"
"ESC_WT_H3K27me3_R2"
"ESC_WT_IGG_R2"
"ESC_KO_EZH1_R3"
"ESC_KO_EZH2_R3"
"ESC_KO_H3K27me3_R3"
"ESC_OEKO_EZH1_R3"
"ESC_OEKO_EZH2_R3"
"ESC_OEKO_H3K27me3_R3"
"ESC_WT_EZH1_R3"
"ESC_WT_EZH2_R3"
"ESC_WT_H3K27me3_R3"
"ESC_WT_IGG_R3"
)


for x in "${x[@]}"; do
    # run bowtie2 and sort read
    bowtie2 --phred33 -q --no-unal --no-mixed --dovetail \
        -x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
            -S output/bowtie2/${x}.sam \
            -1 output/fastp/${x}_1.fq.gz  \
            -2 output/fastp/${x}_2.fq.gz 
done


