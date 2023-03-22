#!/bin/bash

fastqc -o output/fastqc input/2dN_HET_H3K27me3_R1_1.fq.gz
fastqc -o output/fastqc input/2dN_HET_H3K27me3_R1_2.fq.gz
fastqc -o output/fastqc input/2dN_HET_H3K27me3_R2_1.fq.gz
fastqc -o output/fastqc input/2dN_HET_H3K27me3_R2_2.fq.gz
fastqc -o output/fastqc input/2dN_HET_input_R1_1.fq.gz
fastqc -o output/fastqc input/2dN_HET_input_R1_2.fq.gz
fastqc -o output/fastqc input/2dN_HET_input_R2_1.fq.gz
fastqc -o output/fastqc input/2dN_HET_input_R2_2.fq.gz
fastqc -o output/fastqc input/2dN_KO_H3K27me3_R1_1.fq.gz
fastqc -o output/fastqc input/2dN_KO_H3K27me3_R1_2.fq.gz
fastqc -o output/fastqc input/2dN_KO_H3K27me3_R2_1.fq.gz
fastqc -o output/fastqc input/2dN_KO_H3K27me3_R2_2.fq.gz
fastqc -o output/fastqc input/2dN_KO_input_R1_1.fq.gz
fastqc -o output/fastqc input/2dN_KO_input_R1_2.fq.gz
fastqc -o output/fastqc input/2dN_KO_input_R2_1.fq.gz
fastqc -o output/fastqc input/2dN_KO_input_R2_2.fq.gz
fastqc -o output/fastqc input/2dN_WT_H3K27me3_R1_1.fq.gz
fastqc -o output/fastqc input/2dN_WT_H3K27me3_R1_2.fq.gz
fastqc -o output/fastqc input/2dN_WT_H3K27me3_R2_1.fq.gz

fastqc -o output/fastqc input/NPC_KO_input_R2_2.fq.gz
fastqc -o output/fastqc input/NPC_WT_H3K27me3_R1_1.fq.gz
fastqc -o output/fastqc input/NPC_WT_H3K27me3_R1_2.fq.gz
fastqc -o output/fastqc input/NPC_WT_H3K27me3_R2_1.fq.gz
fastqc -o output/fastqc input/NPC_WT_H3K27me3_R2_2.fq.gz
