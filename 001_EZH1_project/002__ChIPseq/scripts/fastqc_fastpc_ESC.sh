#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL



x=("ESC_HET_H3K27me3_R1"
"ESC_HET_H3K27me3_R2"
"ESC_HET_input_R1"
"ESC_HET_input_R2"
"ESC_KO_H3K27me3_R1"
"ESC_KO_H3K27me3_R2"
"ESC_KO_input_R1"
"ESC_KO_input_R2"
"ESC_WT_H3K27me3_R1"
"ESC_WT_H3K27me3_R2"
"ESC_WT_H3K27me3_R3"
"ESC_WT_input_R1"
"ESC_WT_input_R2"
"ESC_WT_input_R3")


for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done