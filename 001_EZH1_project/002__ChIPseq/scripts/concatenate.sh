#!/bin/bash


cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz > input/ESC_KO_input_R2_1.fq.gz
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz > input/ESC_KO_input_R2_2.fq.gz
cat input/10-H3K27me3-Neurons--DID2--Het2A-2_S20_L001_R1_001.fastq.gz input/10-H3K27me3-Neurons--DID2--Het2A-2_S20_L002_R1_001.fastq.gz > input/2dN_HET_H3K27me3_R2_1.fq.gz
cat input/10-H3K27me3-Neurons--DID2--Het2A-2_S20_L001_R2_001.fastq.gz input/10-H3K27me3-Neurons--DID2--Het2A-2_S20_L002_R2_001.fastq.gz > input/2dN_HET_H3K27me3_R2_2.fq.gz

cat input/11-ESCs-Pan-Het-2A-1-INPUT_S9_L001_R1_001.fastq.gz input/11-ESCs-Pan-Het-2A-1-INPUT_S9_L002_R1_001.fastq.gz > input/ESC_HET_input_R1_1.fq.gz
cat input/11-ESCs-Pan-Het-2A-1-INPUT_S9_L001_R2_001.fastq.gz input/11-ESCs-Pan-Het-2A-1-INPUT_S9_L002_R2_001.fastq.gz > input/ESC_HET_input_R1_2.fq.gz

cat input/11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L001_R1_001.fastq.gz input/11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L002_R1_001.fastq.gz > input/2dN_KO_H3K27me3_R1_1.fq.gz
cat input/11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L001_R2_001.fastq.gz input/11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L002_R2_001.fastq.gz > input/2dN_KO_H3K27me3_R1_2.fq.gz

cat input/12-ESCs-Pan-Het-2A-2-INPUT_S10_L001_R1_001.fastq.gz input/12-ESCs-Pan-Het-2A-2-INPUT_S10_L002_R1_001.fastq.gz > input/ESC_HET_input_R2_1.fq.gz
cat input/12-ESCs-Pan-Het-2A-2-INPUT_S10_L001_R2_001.fastq.gz input/12-ESCs-Pan-Het-2A-2-INPUT_S10_L002_R2_001.fastq.gz > input/ESC_HET_input_R2_2.fq.gz

cat input/12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L001_R1_001.fastq.gz input/12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L002_R1_001.fastq.gz > input/2dN_KO_H3K27me3_R2_1.fq.gz
cat input/12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L001_R2_001.fastq.gz input/12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L002_R2_001.fastq.gz > input/2dN_KO_H3K27me3_R2_2.fq.gz

cat input/13-INPUT-NPCs-WTS2-1_S23_L001_R1_001.fastq.gz input/13-INPUT-NPCs-WTS2-1_S23_L002_R1_001.fastq.gz > input/NPC_WT_input_R1_1.fq.gz
cat input/13-INPUT-NPCs-WTS2-1_S23_L001_R2_001.fastq.gz input/13-INPUT-NPCs-WTS2-1_S23_L002_R2_001.fastq.gz > input/NPC_WT_input_R1_2.fq.gz

cat input/14-INPUT-NPCs-WTS2-2_S24_L001_R1_001.fastq.gz input/14-INPUT-NPCs-WTS2-2_S24_L002_R1_001.fastq.gz > input/NPC_WT_input_R2_1.fq.gz
cat input/14-INPUT-NPCs-WTS2-2_S24_L001_R2_001.fastq.gz input/14-INPUT-NPCs-WTS2-2_S24_L002_R2_001.fastq.gz > input/NPC_WT_input_R2_2.fq.gz

cat input/15-INPUT-NPCs-Het2A-1_S25_L001_R1_001.fastq.gz input/15-INPUT-NPCs-Het2A-1_S25_L002_R1_001.fastq.gz > input/NPC_HET_input_R1_1.fq.gz
cat input/15-INPUT-NPCs-Het2A-1_S25_L001_R2_001.fastq.gz input/15-INPUT-NPCs-Het2A-1_S25_L002_R2_001.fastq.gz > input/NPC_HET_input_R1_2.fq.gz

cat input/16-INPUT-NPCs-Het2A-2_S26_L001_R1_001.fastq.gz input/16-INPUT-NPCs-Het2A-2_S26_L002_R1_001.fastq.gz > input/NPC_HET_input_R2_1.fq.gz
cat input/16-INPUT-NPCs-Het2A-2_S26_L001_R2_001.fastq.gz input/16-INPUT-NPCs-Het2A-2_S26_L002_R2_001.fastq.gz > input/NPC_HET_input_R2_2.fq.gz

cat input/17-INPUT-NPCs-EZH1-KO-1_S27_L001_R1_001.fastq.gz input/17-INPUT-NPCs-EZH1-KO-1_S27_L002_R1_001.fastq.gz > input/NPC_KO_input_R1_1.fq.gz
cat input/17-INPUT-NPCs-EZH1-KO-1_S27_L001_R2_001.fastq.gz input/17-INPUT-NPCs-EZH1-KO-1_S27_L002_R2_001.fastq.gz > input/NPC_KO_input_R1_2.fq.gz

cat input/18-INPUT-NPCs-EZH1-KO-2_S28_L001_R1_001.fastq.gz input/18-INPUT-NPCs-EZH1-KO-2_S28_L002_R1_001.fastq.gz > input/NPC_KO_input_R2_1.fq.gz
cat input/18-INPUT-NPCs-EZH1-KO-2_S28_L001_R2_001.fastq.gz input/18-INPUT-NPCs-EZH1-KO-2_S28_L002_R2_001.fastq.gz > input/NPC_KO_input_R2_2.fq.gz

# _2

cat input/19-INPUT-Neurons--DID2--WTS2-1_S29_L001_R1_001.fastq.gz input/19-INPUT-Neurons--DID2--WTS2-1_S29_L002_R1_001.fastq.gz > input/2dN_WT_input_R1_1.fq.gz
cat input/19-INPUT-Neurons--DID2--WTS2-1_S29_L001_R2_001.fastq.gz input/19-INPUT-Neurons--DID2--WTS2-1_S29_L002_R2_001.fastq.gz > input/2dN_WT_input_R1_2.fq.gz

cat input/1-ESCs-WT-S2-1-H3k27me3_S1_L001_R1_001.fastq.gz input/1-ESCs-WT-S2-1-H3k27me3_S1_L002_R1_001.fastq.gz > input/ESC_WT_H3K27me3_R3_1.fq.gz
cat input/1-ESCs-WT-S2-1-H3k27me3_S1_L001_R2_001.fastq.gz input/1-ESCs-WT-S2-1-H3k27me3_S1_L002_R2_001.fastq.gz > input/ESC_WT_H3K27me3_R3_2.fq.gz

cat input/1-H3K27me3-ESCs-H9-WT-1_S35_L001_R1_001.fastq.gz input/1-H3K27me3-ESCs-H9-WT-1_S35_L002_R1_001.fastq.gz > input/ESC_WT_H3K27me3_R1_1.fq.gz
cat input/1-H3K27me3-ESCs-H9-WT-1_S35_L001_R2_001.fastq.gz input/1-H3K27me3-ESCs-H9-WT-1_S35_L002_R2_001.fastq.gz > input/ESC_WT_H3K27me3_R1_2.fq.gz

cat input/1-H3K27me3-NPCs-WTS2-1_S11_L001_R1_001.fastq.gz input/1-H3K27me3-NPCs-WTS2-1_S11_L002_R1_001.fastq.gz > input/NPC_WT_H3K27me3_R1_1.fq.gz
cat input/1-H3K27me3-NPCs-WTS2-1_S11_L001_R2_001.fastq.gz input/1-H3K27me3-NPCs-WTS2-1_S11_L002_R2_001.fastq.gz > input/NPC_WT_H3K27me3_R1_2.fq.gz

cat input/20-INPUT-Neurons--DID2--WTS2-2_S30_L001_R1_001.fastq.gz input/20-INPUT-Neurons--DID2--WTS2-2_S30_L002_R1_001.fastq.gz > input/2dN_WT_input_R2_1.fq.gz
cat input/20-INPUT-Neurons--DID2--WTS2-2_S30_L001_R2_001.fastq.gz input/20-INPUT-Neurons--DID2--WTS2-2_S30_L002_R2_001.fastq.gz > input/2dN_WT_input_R2_2.fq.gz

cat input/21-INPUT-Neurons--DID2--Het2A-1_S31_L001_R1_001.fastq.gz input/21-INPUT-Neurons--DID2--Het2A-1_S31_L002_R1_001.fastq.gz > input/2dN_HET_input_R1_1.fq.gz
cat input/21-INPUT-Neurons--DID2--Het2A-1_S31_L001_R2_001.fastq.gz input/21-INPUT-Neurons--DID2--Het2A-1_S31_L002_R2_001.fastq.gz > input/2dN_HET_input_R1_2.fq.gz

cat input/22-INPUT-Neurons--DID2--Het2A-2_S32_L001_R1_001.fastq.gz input/22-INPUT-Neurons--DID2--Het2A-2_S32_L002_R1_001.fastq.gz > input/2dN_HET_input_R2_1.fq.gz
cat input/22-INPUT-Neurons--DID2--Het2A-2_S32_L001_R2_001.fastq.gz input/22-INPUT-Neurons--DID2--Het2A-2_S32_L002_R2_001.fastq.gz > input/2dN_HET_input_R2_2.fq.gz

cat input/23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L001_R1_001.fastq.gz input/23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L002_R1_001.fastq.gz > input/2dN_KO_input_R1_1.fq.gz
cat input/23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L001_R2_001.fastq.gz input/23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L002_R2_001.fastq.gz > input/2dN_KO_input_R1_2.fq.gz

cat input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R1_001.fastq.gz input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R1_001.fastq.gz > input/2dN_KO_input_R2_1.fq.gz
cat input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq.gz input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq.gz > input/2dN_KO_input_R2_2.fq.gz

cat input/2-H3K27me3-ESCs-H9-WT-2_S36_L001_R1_001.fastq.gz input/2-H3K27me3-ESCs-H9-WT-2_S36_L002_R1_001.fastq.gz > input/ESC_WT_H3K27me3_R2_1.fq.gz
cat input/2-H3K27me3-ESCs-H9-WT-2_S36_L001_R2_001.fastq.gz input/2-H3K27me3-ESCs-H9-WT-2_S36_L002_R2_001.fastq.gz > input/ESC_WT_H3K27me3_R2_2.fq.gz

cat input/2-H3K27me3-NPCs-WTS2-2_S12_L001_R1_001.fastq.gz input/2-H3K27me3-NPCs-WTS2-2_S12_L002_R1_001.fastq.gz > input/NPC_WT_H3K27me3_R2_1.fq.gz
cat input/2-H3K27me3-NPCs-WTS2-2_S12_L001_R2_001.fastq.gz input/2-H3K27me3-NPCs-WTS2-2_S12_L002_R2_001.fastq.gz > input/NPC_WT_H3K27me3_R2_2.fq.gz

cat input/3-ESCs-EZH1-KO-1-H3k27me3_S2_L001_R1_001.fastq.gz input/3-ESCs-EZH1-KO-1-H3k27me3_S2_L002_R1_001.fastq.gz > input/ESC_KO_H3K27me3_R1_1.fq.gz
cat input/3-ESCs-EZH1-KO-1-H3k27me3_S2_L001_R2_001.fastq.gz input/3-ESCs-EZH1-KO-1-H3k27me3_S2_L002_R2_001.fastq.gz > input/ESC_KO_H3K27me3_R1_2.fq.gz

cat input/3-H3K27me3-NPCs-Het2A-1_S13_L001_R1_001.fastq.gz input/3-H3K27me3-NPCs-Het2A-1_S13_L002_R1_001.fastq.gz > input/NPC_HET_H3K27me3_R1_1.fq.gz
cat input/3-H3K27me3-NPCs-Het2A-1_S13_L001_R2_001.fastq.gz input/3-H3K27me3-NPCs-Het2A-1_S13_L002_R2_001.fastq.gz > input/NPC_HET_H3K27me3_R1_2.fq.gz


# _3 

cat input/3-INPUT-ESCs-H9-WT-1_S37_L001_R1_001.fastq.gz input/3-INPUT-ESCs-H9-WT-1_S37_L002_R1_001.fastq.gz > input/ESC_WT_input_R1_1.fq.gz
cat input/3-INPUT-ESCs-H9-WT-1_S37_L001_R2_001.fastq.gz input/3-INPUT-ESCs-H9-WT-1_S37_L002_R2_001.fastq.gz > input/ESC_WT_input_R1_2.fq.gz

cat input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R1_001.fastq.gz input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R1_001.fastq.gz > input/ESC_KO_H3K27me3_R2_1.fq.gz
cat input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R2_001.fastq.gz input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R2_001.fastq.gz > input/ESC_KO_H3K27me3_R2_2.fq.gz

cat input/4-H3K27me3-NPCs-Het2A-2_S14_L001_R1_001.fastq.gz input/4-H3K27me3-NPCs-Het2A-2_S14_L002_R1_001.fastq.gz > input/NPC_HET_H3K27me3_R2_1.fq.gz
cat input/4-H3K27me3-NPCs-Het2A-2_S14_L001_R2_001.fastq.gz input/4-H3K27me3-NPCs-Het2A-2_S14_L002_R2_001.fastq.gz > input/NPC_HET_H3K27me3_R2_2.fq.gz

cat input/4-INPUT-ESCs-H9-WT-2_S38_L001_R1_001.fastq.gz input/4-INPUT-ESCs-H9-WT-2_S38_L002_R1_001.fastq.gz > input/ESC_WT_input_R2_1.fq.gz
cat input/4-INPUT-ESCs-H9-WT-2_S38_L001_R2_001.fastq.gz input/4-INPUT-ESCs-H9-WT-2_S38_L002_R2_001.fastq.gz > input/ESC_WT_input_R2_2.fq.gz

cat input/5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L001_R1_001.fastq.gz input/5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L002_R1_001.fastq.gz > input/ESC_HET_H3K27me3_R1_1.fq.gz
cat input/5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L001_R2_001.fastq.gz input/5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L002_R2_001.fastq.gz > input/ESC_HET_H3K27me3_R1_2.fq.gz

cat input/5-H3K27me3-NPCs-EZH1-KO-1_S15_L001_R1_001.fastq.gz input/5-H3K27me3-NPCs-EZH1-KO-1_S15_L002_R1_001.fastq.gz > input/NPC_KO_H3K27me3_R1_1.fq.gz
cat input/5-H3K27me3-NPCs-EZH1-KO-1_S15_L001_R2_001.fastq.gz input/5-H3K27me3-NPCs-EZH1-KO-1_S15_L002_R2_001.fastq.gz > input/NPC_KO_H3K27me3_R1_2.fq.gz

cat input/6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L001_R1_001.fastq.gz input/6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L002_R1_001.fastq.gz > input/ESC_HET_H3K27me3_R2_1.fq.gz
cat input/6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L001_R2_001.fastq.gz input/6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L002_R2_001.fastq.gz > input/ESC_HET_H3K27me3_R2_2.fq.gz

cat input/6-H3K27me3-NPCs-EZH1-KO-2_S16_L001_R1_001.fastq.gz input/6-H3K27me3-NPCs-EZH1-KO-2_S16_L002_R1_001.fastq.gz > input/NPC_KO_H3K27me3_R2_1.fq.gz
cat input/6-H3K27me3-NPCs-EZH1-KO-2_S16_L001_R2_001.fastq.gz input/6-H3K27me3-NPCs-EZH1-KO-2_S16_L002_R2_001.fastq.gz > input/NPC_KO_H3K27me3_R2_2.fq.gz

cat input/7-ESCs-WT-S2-1-INPUT_S6_L001_R1_001.fastq.gz input/7-ESCs-WT-S2-1-INPUT_S6_L002_R1_001.fastq.gz > input/ESC_WT_input_R3_1.fq.gz
cat input/7-ESCs-WT-S2-1-INPUT_S6_L001_R2_001.fastq.gz input/7-ESCs-WT-S2-1-INPUT_S6_L002_R2_001.fastq.gz > input/ESC_WT_input_R3_2.fq.gz

cat input/7-H3K27me3-Neurons--DID2--WTS2-1_S17_L001_R1_001.fastq.gz input/7-H3K27me3-Neurons--DID2--WTS2-1_S17_L002_R1_001.fastq.gz > input/2dN_WT_H3K27me3_R1_1.fq.gz
cat input/7-H3K27me3-Neurons--DID2--WTS2-1_S17_L001_R2_001.fastq.gz input/7-H3K27me3-Neurons--DID2--WTS2-1_S17_L002_R2_001.fastq.gz > input/2dN_WT_H3K27me3_R1_2.fq.gz

cat input/8-H3K27me3-Neurons--DID2--WTS2-2_S18_L001_R1_001.fastq.gz input/8-H3K27me3-Neurons--DID2--WTS2-2_S18_L002_R1_001.fastq.gz > input/2dN_WT_H3K27me3_R2_1.fq.gz
cat input/8-H3K27me3-Neurons--DID2--WTS2-2_S18_L001_R2_001.fastq.gz input/8-H3K27me3-Neurons--DID2--WTS2-2_S18_L002_R2_001.fastq.gz > input/2dN_WT_H3K27me3_R2_2.fq.gz

cat input/9-ESCs-EZH1-KO-1-INPUT_S7_L001_R1_001.fastq.gz input/9-ESCs-EZH1-KO-1-INPUT_S7_L002_R1_001.fastq.gz > input/ESC_KO_input_R1_1.fq.gz
cat input/9-ESCs-EZH1-KO-1-INPUT_S7_L001_R2_001.fastq.gz input/9-ESCs-EZH1-KO-1-INPUT_S7_L002_R2_001.fastq.gz > input/ESC_KO_input_R1_2.fq.gz

cat input/9-H3K27me3-Neurons--DID2--Het2A-1_S19_L001_R1_001.fastq.gz input/9-H3K27me3-Neurons--DID2--Het2A-1_S19_L002_R1_001.fastq.gz > input/2dN_HET_H3K27me3_R1_1.fq.gz
cat input/9-H3K27me3-Neurons--DID2--Het2A-1_S19_L001_R2_001.fastq.gz input/9-H3K27me3-Neurons--DID2--Het2A-1_S19_L002_R2_001.fastq.gz > input/2dN_HET_H3K27me3_R1_2.fq.gz


