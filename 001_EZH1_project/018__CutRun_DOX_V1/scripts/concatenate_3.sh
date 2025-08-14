#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00






cat input_raw_Novogene/R1_ESC_OEKO_EZH1_CKDL250021971-1A_22W5G2LT4_L2_1.fq.gz \
input_raw_Novogene/R1_ESC_OEKO_EZH1_CKDL250021971-1A_22W5TMLT4_L3_1.fq.gz > \
input/ESC_OEKO_EZH1_R1_1.fq.gz



